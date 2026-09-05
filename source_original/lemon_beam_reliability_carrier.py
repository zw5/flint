#!/usr/bin/env python3
"""LEMON beam reliability and carrier-alignment ceiling pass.

The current beam result appears to plateau near r=.40 in sensor space. This
script tests a lawful way to improve the measurement without adding arbitrary
feature knobs:

1. Estimate each subject's posterior/lateral theta carrier frequency.
2. Compare fixed 6 Hz phase lags against subject-carrier-aligned phase lags.
3. Estimate beam-score reliability from odd/even and first/second window splits.
4. Evaluate whether reliability/precision weighting raises the association
   without raising the WST vocabulary control association.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from collections import Counter
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
from scipy.signal import butter, filtfilt, hilbert, welch
from scipy.stats import rankdata, t as tdist

import lemon_intelligence_residual_search as stats
from lemon_beam_hypothesis_dossier import build_apertures
from lemon_frequency_sliced_beam_eval import concentration_from_spectrum, ffloat, static_plv, window_starts


warnings.filterwarnings("ignore")

ROOT = Path("/Users/ximon/Documents/interp")
EXP = ROOT / "experiments/eeg_coherence"
EXPORT = EXP / "rust_export"
OUT = EXP / "lemon_beam_reliability_carrier"

THETA = (4.0, 8.0)
FIXED_CARRIER_HZ = 6.0
PHASE_FRACTIONS = (0.25, 1.0 / 3.0, 0.5)
APERTURES = ("posterior_broad", "lateral_posterior")
MODES = ("fixed6", "carrier")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.nan_to_num((x - np.nanmean(x)) / (np.nanstd(x) + 1e-9))


def cohort(age: float) -> str:
    if age <= 35:
        return "young"
    if age >= 55:
        return "old"
    return "middle"


def load_manifest(limit: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_csv(EXPORT / "manifest.csv"):
        age = ffloat(row.get("age"))
        c = cohort(age)
        if c == "middle":
            continue
        rec = dict(row)
        rec["age"] = age
        rec["LPS"] = ffloat(row.get("LPS"))
        rec["WST"] = ffloat(row.get("WST"))
        rec["sfreq"] = ffloat(row.get("sfreq"))
        rec["n_ch"] = int(float(row["n_ch"]))
        rec["n_eo"] = int(float(row["n_eo"]))
        rec["cohort"] = c
        rec["cohort_code"] = 0.0 if c == "young" else 1.0
        if rec["n_eo"] > 0 and all(math.isfinite(float(rec[k])) for k in ("age", "LPS", "WST", "sfreq")):
            rows.append(rec)
    rows.sort(key=lambda r: str(r["sid"]))
    return rows[:limit] if limit else rows


def load_channels() -> list[str]:
    return json.loads((EXPORT / "channels.json").read_text())["channels"]


def load_eo(row: dict[str, object]) -> np.ndarray | None:
    sid = str(row["sid"])
    n = int(row["n_eo"])
    path = EXPORT / f"{sid}_EO.f32"
    if n <= 0 or not path.exists():
        return None
    n_ch = int(row["n_ch"])
    data = np.fromfile(path, dtype=np.float32)
    if data.size != n_ch * n:
        return None
    return data.reshape(n_ch, n).T.astype(np.float64, copy=False)


def bandpass_theta(X: np.ndarray, sfreq: float) -> np.ndarray:
    lo, hi = THETA
    hi = min(hi, sfreq / 2.0 - 1.0)
    b, a = butter(4, [lo / (sfreq / 2.0), hi / (sfreq / 2.0)], btype="band")
    return filtfilt(b, a, X, axis=0)


def analytic_phase(X: np.ndarray) -> np.ndarray:
    analytic = hilbert(X, axis=0)
    return analytic / (np.abs(analytic) + 1e-9)


def theta_peak_frequency(X: np.ndarray, sfreq: float, lo: float = 4.0, hi: float = 8.0) -> float:
    X = np.asarray(X, dtype=float)
    nperseg = min(len(X), max(256, int(round(4.0 * sfreq))))
    freqs, psd = welch(X, fs=sfreq, nperseg=nperseg, axis=0)
    band = (freqs >= lo) & (freqs <= hi)
    if not np.any(band):
        return 0.5 * (lo + hi)
    mean_psd = np.nanmean(psd[band], axis=1)
    if not np.any(np.isfinite(mean_psd)):
        return 0.5 * (lo + hi)
    return float(freqs[band][int(np.nanargmax(mean_psd))])


def phase_lag_samples(sfreq: float, carrier_hz: float, phase_fraction: float) -> int:
    carrier_hz = min(max(float(carrier_hz), THETA[0]), THETA[1])
    return max(1, int(round((phase_fraction / carrier_hz) * sfreq)))


def cross_concentration(Z: np.ndarray, lag: int) -> float:
    if lag <= 0:
        C = (Z.conj().T @ Z) / max(len(Z), 1)
    else:
        present = Z[:-lag]
        future = Z[lag:]
        C = (future.T @ present.conj()) / max(len(present), 1)
    s = np.linalg.svd(C, compute_uv=False)
    concentration, _effrank, _lock = concentration_from_spectrum(s)
    return concentration


def score_for_starts(phase: np.ndarray, lags: dict[str, int], starts: list[int], win: int) -> dict[str, float]:
    by_lag: dict[str, list[float]] = {name: [] for name in lags}
    window_scores: list[float] = []
    for start in starts:
        Z = phase[start : start + win]
        lag_scores = []
        for name, lag in lags.items():
            if len(Z) <= lag + 10:
                continue
            try:
                c = cross_concentration(Z, lag)
            except np.linalg.LinAlgError:
                continue
            if math.isfinite(c):
                by_lag[name].append(c)
                lag_scores.append(c)
        if lag_scores:
            window_scores.append(float(np.nanmean(lag_scores)))
    out = {f"{name}_concentration": float(np.nanmean(vals)) if vals else math.nan for name, vals in by_lag.items()}
    out["score"] = float(np.nanmean(window_scores)) if window_scores else math.nan
    out["window_sd"] = float(np.nanstd(window_scores)) if window_scores else math.nan
    out["n_score_windows"] = float(len(window_scores))
    return out


def split_start_sets(starts: list[int]) -> dict[str, list[int]]:
    half = len(starts) // 2
    return {
        "full": starts,
        "odd": starts[::2],
        "even": starts[1::2],
        "first": starts[:half],
        "second": starts[half:],
    }


def spearman_brown(r: float) -> float:
    if not math.isfinite(float(r)):
        return math.nan
    if abs(1.0 + r) <= 1e-12:
        return math.nan
    return float((2.0 * r) / (1.0 + r))


def weighted_rank_residual(x: np.ndarray, controls: list[np.ndarray], weights: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(x) & np.isfinite(weights) & (weights > 0)
    for c in controls:
        mask &= np.isfinite(c)
    out = np.full_like(x, np.nan, dtype=float)
    if int(mask.sum()) <= len(controls) + 3:
        return out
    xr = rankdata(x[mask])
    w = weights[mask]
    w = w / (np.nanmean(w) + 1e-12)
    if not controls:
        out[mask] = xr - np.average(xr, weights=w)
        return out
    C = np.column_stack([rankdata(np.asarray(c)[mask]) for c in controls])
    B = np.column_stack([np.ones(len(xr)), C])
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(B * sw[:, None], xr * sw, rcond=None)
    out[mask] = xr - B @ coef
    return out


def weighted_partial_spearman(
    x: np.ndarray,
    y: np.ndarray,
    controls: list[np.ndarray],
    weights: np.ndarray,
) -> tuple[float, float, int]:
    rx = weighted_rank_residual(x, controls, weights)
    ry = weighted_rank_residual(y, controls, weights)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(rx) & np.isfinite(ry) & np.isfinite(weights) & (weights > 0)
    n = int(mask.sum())
    if n <= len(controls) + 3:
        return math.nan, math.nan, n
    rx = rx[mask]
    ry = ry[mask]
    w = weights[mask]
    w = w / (np.nanmean(w) + 1e-12)
    mx = np.average(rx, weights=w)
    my = np.average(ry, weights=w)
    cov = np.average((rx - mx) * (ry - my), weights=w)
    vx = np.average((rx - mx) ** 2, weights=w)
    vy = np.average((ry - my) ** 2, weights=w)
    r = float(cov / math.sqrt(max(vx * vy, 1e-12)))
    dof = n - len(controls) - 2
    p = float(2 * tdist.sf(abs(r * np.sqrt(dof / max(1 - r * r, 1e-9))), dof))
    return r, p, n


def process_subject(job: tuple[dict[str, object], dict[str, object]]) -> dict[str, object] | None:
    row, cfg = job
    sid = str(row["sid"])
    X = load_eo(row)
    if X is None:
        return None
    sfreq = float(row["sfreq"])
    channels: list[str] = cfg["channels"]
    apertures = build_apertures(channels)
    win = int(round(float(cfg["window_sec"]) * sfreq))
    step = int(round(float(cfg["step_sec"]) * sfreq))
    starts = window_starts(len(X), win, step, int(cfg["max_windows"]))
    start_sets = split_start_sets(starts)
    try:
        Xb = bandpass_theta(X, sfreq)
        phase_all = analytic_phase(Xb)
    except Exception as exc:
        print(f"  {sid}: {exc}", flush=True)
        return None

    out: dict[str, object] = {
        "sid": sid,
        "LPS": float(row["LPS"]),
        "WST": float(row["WST"]),
        "age": float(row["age"]),
        "cohort": row["cohort"],
        "cohort_code": float(row["cohort_code"]),
    }

    for aperture in APERTURES:
        idx = apertures[aperture]
        if len(idx) < 3:
            continue
        Xsel = Xb[:, idx]
        phase = phase_all[:, idx]
        peak = theta_peak_frequency(Xsel, sfreq, *THETA)
        out[f"{aperture}__theta_peak_hz"] = peak
        out[f"{aperture}__theta_power"] = float(np.mean(Xsel**2))
        out[f"{aperture}__theta_static_plv"] = static_plv(phase)
        carriers = {"fixed6": FIXED_CARRIER_HZ, "carrier": peak}
        for mode, carrier_hz in carriers.items():
            lags = {f"p{int(round(frac * 1000)):03d}": phase_lag_samples(sfreq, carrier_hz, frac) for frac in PHASE_FRACTIONS}
            prefix = f"{aperture}_{mode}"
            for tag, lag in lags.items():
                out[f"{prefix}__{tag}_lag_samples"] = lag
                out[f"{prefix}__{tag}_lag_ms"] = float(1000.0 * lag / sfreq)
            for split_name, split_starts in start_sets.items():
                scores = score_for_starts(phase, lags, split_starts, win)
                for key, value in scores.items():
                    out[f"{prefix}__{split_name}_{key}"] = value
    return out


def eval_feature(rows: list[dict[str, object]], feature: str, controls: list[str]) -> dict[str, object]:
    x = np.asarray([ffloat(r.get(feature)) for r in rows])
    y = np.asarray([ffloat(r.get("LPS")) for r in rows])
    wst = np.asarray([ffloat(r.get("WST")) for r in rows])
    ctrl = [np.asarray([ffloat(r.get(c)) for r in rows]) for c in controls]
    r, p, n = stats.partial_spearman(x, y, ctrl)
    rw, pw, _ = stats.partial_spearman(x, wst, ctrl)
    rlw, plw, _ = stats.partial_spearman(x, y, ctrl + [wst])
    return {
        "feature": feature,
        "n": n,
        "controls": ";".join(controls),
        "r_lps": r,
        "p_lps": p,
        "r_lps_plus_wst": rlw,
        "p_lps_plus_wst": plw,
        "r_wst": rw,
        "p_wst": pw,
    }


def controls_for(aperture: str, include_peak: bool) -> list[str]:
    controls = [
        "cohort_code",
        "age",
        f"{aperture}__theta_power",
        f"{aperture}__theta_static_plv",
    ]
    if include_peak:
        controls.append(f"{aperture}__theta_peak_hz")
    return controls


def reliability_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for aperture in APERTURES:
        for mode in MODES:
            prefix = f"{aperture}_{mode}"
            odd = np.asarray([ffloat(r.get(f"{prefix}__odd_score")) for r in rows])
            even = np.asarray([ffloat(r.get(f"{prefix}__even_score")) for r in rows])
            first = np.asarray([ffloat(r.get(f"{prefix}__first_score")) for r in rows])
            second = np.asarray([ffloat(r.get(f"{prefix}__second_score")) for r in rows])
            for split_name, a, b in (("odd_even", odd, even), ("first_second", first, second)):
                mask = np.isfinite(a) & np.isfinite(b)
                r = float(np.corrcoef(rankdata(a[mask]), rankdata(b[mask]))[0, 1]) if int(mask.sum()) >= 5 else math.nan
                out.append(
                    {
                        "aperture": aperture,
                        "mode": mode,
                        "split": split_name,
                        "n": int(mask.sum()),
                        "split_r": r,
                        "spearman_brown": spearman_brown(r),
                    }
                )
    return out


def precision_weights(rows: list[dict[str, object]], prefix: str) -> np.ndarray:
    sd = np.asarray([ffloat(r.get(f"{prefix}__full_window_sd")) for r in rows], dtype=float)
    med = np.nanmedian(sd[np.isfinite(sd)]) if np.any(np.isfinite(sd)) else 1.0
    sd = np.where(np.isfinite(sd), sd, med)
    raw = 1.0 / (sd + 1e-4) ** 2
    lo, hi = np.nanquantile(raw, [0.05, 0.95])
    raw = np.clip(raw, lo, hi)
    return raw / (np.nanmean(raw) + 1e-12)


def evaluate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for aperture in APERTURES:
        for mode in MODES:
            prefix = f"{aperture}_{mode}"
            features = [f"{prefix}__full_score"] + [f"{prefix}__full_p{int(round(frac * 1000)):03d}_concentration" for frac in PHASE_FRACTIONS]
            for include_peak in (False, True):
                analysis = "strict_peak" if include_peak else "strict"
                controls = controls_for(aperture, include_peak)
                for feature in features:
                    out.append({"analysis": analysis, **eval_feature(rows, feature, controls)})

            weights = precision_weights(rows, prefix)
            x = np.asarray([ffloat(r.get(f"{prefix}__full_score")) for r in rows])
            y = np.asarray([ffloat(r.get("LPS")) for r in rows])
            wst = np.asarray([ffloat(r.get("WST")) for r in rows])
            ctrl = [np.asarray([ffloat(r.get(c)) for r in rows]) for c in controls_for(aperture, include_peak=False)]
            rw, pw, n = weighted_partial_spearman(x, y, ctrl, weights)
            rww, pww, _ = weighted_partial_spearman(x, wst, ctrl, weights)
            rlw, plw, _ = weighted_partial_spearman(x, y, ctrl + [wst], weights)
            out.append(
                {
                    "analysis": "weighted_precision",
                    "feature": f"{prefix}__full_score",
                    "n": n,
                    "controls": ";".join(controls_for(aperture, include_peak=False)),
                    "r_lps": rw,
                    "p_lps": pw,
                    "r_lps_plus_wst": rlw,
                    "p_lps_plus_wst": plw,
                    "r_wst": rww,
                    "p_wst": pww,
                }
            )
    stats.bh_q(out, "p_lps", "q_lps")
    return out


def feature_aperture(feature: str) -> str:
    return "lateral_posterior" if feature.startswith("lateral_posterior_") else "posterior_broad"


def split_stability(
    rows: list[dict[str, object]],
    features: list[str],
    n_splits: int,
    seed: int,
) -> list[dict[str, object]]:
    if n_splits <= 0:
        return []
    rng = np.random.default_rng(seed)
    idx = np.arange(len(rows))
    out: list[dict[str, object]] = []
    for feature in features:
        aperture = feature_aperture(feature)
        controls = controls_for(aperture, include_peak=False)
        train_rs = []
        held_rs = []
        held_ps = []
        held_wst_rs = []
        for _ in range(n_splits):
            perm = rng.permutation(idx)
            half = len(perm) // 2
            train = [rows[i] for i in perm[:half]]
            held = [rows[i] for i in perm[half:]]
            tr = eval_feature(train, feature, controls)
            he = eval_feature(held, feature, controls)
            train_rs.append(ffloat(tr["r_lps"]))
            held_rs.append(ffloat(he["r_lps"]))
            held_ps.append(ffloat(he["p_lps"]))
            held_wst_rs.append(ffloat(he["r_wst"]))
        train_arr = np.asarray(train_rs, dtype=float)
        held_arr = np.asarray(held_rs, dtype=float)
        held_p_arr = np.asarray(held_ps, dtype=float)
        out.append(
            {
                "feature": feature,
                "n_splits": n_splits,
                "train_r_median": float(np.nanmedian(train_arr)),
                "held_r_median": float(np.nanmedian(held_arr)),
                "held_r_q05": float(np.nanquantile(held_arr, 0.05)),
                "held_r_q95": float(np.nanquantile(held_arr, 0.95)),
                "same_direction_rate": float(np.nanmean(np.sign(train_arr) == np.sign(held_arr))),
                "held_p_lt_0_05_rate": float(np.nanmean(held_p_arr < 0.05)),
                "held_wst_r_median": float(np.nanmedian(held_wst_rs)),
            }
        )
    return out


def write_summary(
    path: Path,
    subjects: list[dict[str, object]],
    eval_rows: list[dict[str, object]],
    rel_rows: list[dict[str, object]],
    split_rows: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    def fmt(row: dict[str, object]) -> str:
        return (
            "- {analysis} {feature}: n={n}, LPS r={r_lps:+.3f} p={p_lps:.4g} q={q_lps:.4g}, "
            "LPS+WST r={r_lps_plus_wst:+.3f} p={p_lps_plus_wst:.4g}, WST r={r_wst:+.3f} p={p_wst:.4g}"
        ).format(
            analysis=row["analysis"],
            feature=row["feature"],
            n=row["n"],
            r_lps=ffloat(row["r_lps"]),
            p_lps=ffloat(row["p_lps"]),
            q_lps=ffloat(row.get("q_lps")),
            r_lps_plus_wst=ffloat(row["r_lps_plus_wst"]),
            p_lps_plus_wst=ffloat(row["p_lps_plus_wst"]),
            r_wst=ffloat(row["r_wst"]),
            p_wst=ffloat(row["p_wst"]),
        )

    primary_features = {
        f"{aperture}_{mode}__full_score"
        for aperture in APERTURES
        for mode in MODES
    }
    primary = [r for r in eval_rows if r["feature"] in primary_features]
    top = sorted(eval_rows, key=lambda r: ffloat(r["p_lps"]) if math.isfinite(ffloat(r["p_lps"])) else 1.0)[:24]
    lines = [
        "# LEMON Beam Reliability and Carrier Alignment",
        "",
        "Goal: test whether subject-specific theta carrier alignment improves the same low-entropy self-propagation beam object.",
        "",
        f"Subjects: {len(subjects)}",
        f"Cohorts: {dict(Counter(str(r['cohort']) for r in subjects))}",
        f"Window seconds: {args.window_sec}",
        f"Max windows: {args.max_windows}",
        "",
        "## Primary Fixed-vs-Carrier Rows",
        "",
    ]
    for row in sorted(primary, key=lambda r: (str(r["feature"]), str(r["analysis"]))):
        lines.append(fmt(row))
    lines.extend(["", "## Reliability", ""])
    for row in rel_rows:
        lines.append(
            "- {aperture} {mode} {split}: n={n}, split r={split_r:+.3f}, Spearman-Brown={spearman_brown:+.3f}".format(
                aperture=row["aperture"],
                mode=row["mode"],
                split=row["split"],
                n=row["n"],
                split_r=ffloat(row["split_r"]),
                spearman_brown=ffloat(row["spearman_brown"]),
            )
        )
    lines.extend(["", "## Split Stability", ""])
    for row in split_rows:
        lines.append(
            "- {feature}: held r median={held_r_median:+.3f}, 5-95%=[{held_r_q05:+.3f},{held_r_q95:+.3f}], "
            "same direction={same_direction_rate:.1%}, held p<0.05={held_p_lt_0_05_rate:.1%}, held WST r={held_wst_r_median:+.3f}".format(
                **row
            )
        )
    lines.extend(["", "## Best Rows", ""])
    for row in top:
        lines.append(fmt(row))
    lines.extend(
        [
            "",
            "## Reading Rule",
            "",
            "Carrier alignment is promoted only if it improves reliability and LPS association without lifting WST. Precision weighting is exploratory because window variance may include real dynamics as well as measurement noise.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=max(1, min(5, cpu_count() - 2)))
    parser.add_argument("--window-sec", type=float, default=8.0)
    parser.add_argument("--step-sec", type=float, default=4.0)
    parser.add_argument("--max-windows", type=int, default=24)
    parser.add_argument("--splits", type=int, default=300)
    parser.add_argument("--tag", default="full")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    suffix = "full" if not args.limit else f"limit{args.limit}"
    if args.tag:
        suffix = f"{suffix}_{args.tag}"
    OUT.mkdir(parents=True, exist_ok=True)
    subj_path = OUT / f"carrier_subjects_{suffix}.csv"
    eval_path = OUT / f"carrier_eval_{suffix}.csv"
    rel_path = OUT / f"carrier_reliability_{suffix}.csv"
    split_path = OUT / f"carrier_splits_{suffix}.csv"
    summary_path = OUT / f"SUMMARY_{suffix}.md"

    if subj_path.exists() and not args.force:
        subjects = read_csv(subj_path)
        for row in subjects:
            for key, value in list(row.items()):
                if key not in {"sid", "cohort"}:
                    row[key] = ffloat(value)
    else:
        manifest = load_manifest(args.limit)
        cfg = {
            "channels": load_channels(),
            "window_sec": args.window_sec,
            "step_sec": args.step_sec,
            "max_windows": args.max_windows,
        }
        print(f"Subjects queued: {len(manifest)}; workers={args.workers}; max_windows={args.max_windows}", flush=True)
        jobs = [(row, cfg) for row in manifest]
        subjects: list[dict[str, object]] = []
        if args.workers <= 1:
            for i, job in enumerate(jobs, start=1):
                print(f"  processing {i}/{len(jobs)} {job[0]['sid']}", flush=True)
                rec = process_subject(job)
                if rec is not None:
                    subjects.append(rec)
        else:
            with Pool(args.workers) as pool:
                for i, rec in enumerate(pool.imap_unordered(process_subject, jobs), start=1):
                    if rec is not None:
                        subjects.append(rec)
                    if i == 1 or i % 10 == 0:
                        print(f"  completed {i}/{len(jobs)} ({len(subjects)} subjects)", flush=True)
        subjects.sort(key=lambda r: str(r["sid"]))
        write_csv(subj_path, subjects)

    eval_rows = evaluate(subjects)
    rel_rows = reliability_rows(subjects)
    split_features = [
        f"{aperture}_{mode}__full_score"
        for aperture in APERTURES
        for mode in MODES
    ] + [
        "posterior_broad_carrier__full_p500_concentration",
        "lateral_posterior_carrier__full_p500_concentration",
    ]
    split_rows = split_stability(subjects, split_features, args.splits, 20260627)
    write_csv(eval_path, eval_rows)
    write_csv(rel_path, rel_rows)
    write_csv(split_path, split_rows)
    write_summary(summary_path, subjects, eval_rows, rel_rows, split_rows, args)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
