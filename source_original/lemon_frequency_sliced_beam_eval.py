#!/usr/bin/env python3
"""Frequency-sliced LEMON beam concentration evaluation.

This tests the concern that fixed millisecond lags such as 40 ms and 80 ms are
arbitrary for theta waves. Instead of using fixed time lags, this script slices
theta into narrow frequency bands and evaluates lags as fractions of each wave
cycle:

    tau = phase_fraction / band_center_frequency

The primary score is the cross-field spectral concentration that survived the
proper-operator run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
from scipy.signal import butter, filtfilt, hilbert

import lemon_intelligence_residual_search as stats


warnings.filterwarnings("ignore")

ROOT = Path("/Users/ximon/Documents/interp")
EXP = ROOT / "experiments/eeg_coherence"
EXPORT = EXP / "rust_export"
OUT = EXP / "lemon_frequency_sliced_beam"

BANDS = {
    "theta4_8": (4.0, 8.0),
    "theta4_5": (4.0, 5.0),
    "theta5_6": (5.0, 6.0),
    "theta6_7": (6.0, 7.0),
    "theta7_8": (7.0, 8.0),
}
PHASE_FRACTIONS = (0.25, 1.0 / 3.0, 0.5)
CHSETS = ("posterior", "all")


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


def ffloat(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.nan_to_num((x - np.nanmean(x)) / (np.nanstd(x) + 1e-9))


def cohort(age: float) -> str:
    if age <= 35:
        return "young"
    if age >= 55:
        return "old"
    return "middle"


def is_frontal(ch: str) -> bool:
    ch = ch.upper()
    return (
        ch.startswith("FP")
        or ch.startswith("AF")
        or (len(ch) >= 2 and ch[0] == "F" and ch[1] in "12345678Z")
    )


def channel_sets() -> dict[str, list[int]]:
    channels = json.loads((EXPORT / "channels.json").read_text())["channels"]
    posterior = [i for i, ch in enumerate(channels) if not is_frontal(ch)]
    return {"posterior": posterior, "all": list(range(len(channels)))}


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


def bandpass(X: np.ndarray, sfreq: float, lo: float, hi: float) -> np.ndarray:
    hi = min(hi, sfreq / 2.0 - 1.0)
    b, a = butter(4, [lo / (sfreq / 2.0), hi / (sfreq / 2.0)], btype="band")
    return filtfilt(b, a, X, axis=0)


def static_plv(phase: np.ndarray) -> float:
    plv = np.abs((phase.conj().T @ phase) / max(len(phase), 1))
    iu = np.triu_indices(plv.shape[0], 1)
    return float(plv[iu].mean()) if len(iu[0]) else math.nan


def lag_samples_for_phase_fraction(sfreq: float, lo: float, hi: float, phase_fraction: float) -> int:
    center = 0.5 * (lo + hi)
    tau_sec = phase_fraction / center
    return max(1, int(round(tau_sec * sfreq)))


def concentration_from_spectrum(s: np.ndarray) -> tuple[float, float, float]:
    energy = np.maximum(np.asarray(s, dtype=float), 0.0) ** 2
    total = float(energy.sum())
    if total <= 1e-12:
        return math.nan, math.nan, math.nan
    p = energy / total
    p_pos = p[p > 0]
    entropy = float(-(p_pos * np.log(p_pos)).sum())
    max_entropy = math.log(len(p)) if len(p) > 1 else 1.0
    concentration = 1.0 - entropy / max(max_entropy, 1e-12)
    effrank = float(np.exp(entropy))
    lock = float(p[0])
    return concentration, effrank, lock


def cross_spectrum(Z: np.ndarray, lag: int) -> np.ndarray:
    present = Z[:-lag]
    future = Z[lag:]
    return (future.T @ present.conj()) / max(len(present), 1)


def leading_mode_metrics(Z: np.ndarray, lag: int) -> tuple[dict[str, float], np.ndarray]:
    C = cross_spectrum(Z, lag)
    u, s, _vh = np.linalg.svd(C, full_matrices=False)
    concentration, effrank, lock = concentration_from_spectrum(s)
    mode = u[:, 0] / (np.linalg.norm(u[:, 0]) + 1e-12)
    return {"concentration": concentration, "effrank": effrank, "lock": lock}, mode


def mean_adjacent_alignment(modes: list[np.ndarray]) -> float:
    if len(modes) < 2:
        return math.nan
    return float(np.nanmean([abs(np.vdot(modes[i], modes[i + 1])) for i in range(len(modes) - 1)]))


def mean_pairwise_lock(modes: list[np.ndarray]) -> float:
    if len(modes) < 2:
        return math.nan
    vals = []
    for i in range(len(modes)):
        for j in range(i + 1, len(modes)):
            vals.append(abs(np.vdot(modes[i], modes[j])))
    return float(np.nanmean(vals))


def window_starts(n: int, win: int, step: int, max_windows: int) -> list[int]:
    starts = list(range(0, max(n - win + 1, 0), step))
    if not starts:
        return []
    if max_windows and len(starts) > max_windows:
        idx = np.linspace(0, len(starts) - 1, max_windows).round().astype(int)
        starts = [starts[i] for i in idx]
    return starts


def phase_lag_features(phase: np.ndarray, lag: int, starts: list[int], win: int) -> dict[str, float]:
    metrics: dict[str, list[float]] = defaultdict(list)
    modes: list[np.ndarray] = []
    for start in starts:
        Z = phase[start : start + win]
        if len(Z) <= lag + 10:
            continue
        try:
            observed, mode = leading_mode_metrics(Z, lag)
        except np.linalg.LinAlgError:
            continue
        if not all(math.isfinite(v) for v in observed.values()):
            continue
        for key, value in observed.items():
            metrics[key].append(value)
        modes.append(mode)
    if not metrics["concentration"] or not modes:
        return {}
    out = {key: float(np.nanmean(vals)) for key, vals in metrics.items()}
    out["local_align"] = mean_adjacent_alignment(modes)
    out["global_lock"] = mean_pairwise_lock(modes)
    out["bint_ca"] = out["concentration"] * out["local_align"]
    out["n_lag_windows"] = float(len(modes))
    return out


def process_subject(job: tuple[dict[str, object], dict[str, object]]) -> list[dict[str, object]]:
    row, cfg = job
    sid = str(row["sid"])
    X = load_eo(row)
    if X is None:
        return []
    sfreq = float(row["sfreq"])
    win = int(round(float(cfg["window_sec"]) * sfreq))
    step = int(round(float(cfg["step_sec"]) * sfreq))
    chidx = cfg["chidx"]

    out: list[dict[str, object]] = []
    for band_name, (lo, hi) in BANDS.items():
        try:
            Xb = bandpass(X, sfreq, lo, hi)
            analytic = hilbert(Xb, axis=0)
            phase_all = analytic / (np.abs(analytic) + 1e-9)
        except Exception as exc:
            print(f"  {sid} {band_name}: {exc}", flush=True)
            continue
        starts = window_starts(len(phase_all), win, step, int(cfg["max_windows"]))
        for chset in CHSETS:
            idx = chidx[chset]
            phase = phase_all[:, idx]
            Xsel = Xb[:, idx]
            base = {
                "sid": sid,
                "LPS": float(row["LPS"]),
                "WST": float(row["WST"]),
                "age": float(row["age"]),
                "cohort": row["cohort"],
                "cohort_code": float(row["cohort_code"]),
                "band": band_name,
                "band_lo": lo,
                "band_hi": hi,
                "band_center": 0.5 * (lo + hi),
                "channel_set": chset,
                "n_channels": len(idx),
                "n_times": len(phase),
                "n_windows": len(starts),
                "band_power": float(np.mean(Xsel**2)),
                "static_plv": static_plv(phase),
            }
            for phase_fraction in PHASE_FRACTIONS:
                lag = lag_samples_for_phase_fraction(sfreq, lo, hi, phase_fraction)
                feats = phase_lag_features(phase, lag, starts, win)
                if not feats:
                    continue
                rec = dict(base)
                rec["phase_fraction"] = phase_fraction
                rec["lag_samples"] = lag
                rec["lag_ms"] = float(1000.0 * lag / sfreq)
                rec.update(feats)
                out.append(rec)
    return out


def pivot_subjects(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_sid: dict[str, dict[str, object]] = {}
    for row in rows:
        sid = str(row["sid"])
        out = by_sid.setdefault(
            sid,
            {
                "sid": sid,
                "LPS": row["LPS"],
                "WST": row["WST"],
                "age": row["age"],
                "cohort": row["cohort"],
                "cohort_code": row["cohort_code"],
            },
        )
        frac_tag = f"p{int(round(float(row['phase_fraction']) * 1000)):03d}"
        prefix = f"{row['channel_set']}_{row['band']}_{frac_tag}"
        for key, value in row.items():
            if key in {
                "sid",
                "LPS",
                "WST",
                "age",
                "cohort",
                "cohort_code",
                "channel_set",
                "band",
                "phase_fraction",
            }:
                continue
            out[f"{prefix}__{key}"] = value
        chband = f"{row['channel_set']}_{row['band']}"
        out[f"{chband}__band_power"] = row["band_power"]
        out[f"{chband}__static_plv"] = row["static_plv"]
    return [by_sid[sid] for sid in sorted(by_sid)]


def add_derived_candidates(subjects: list[dict[str, object]]) -> list[str]:
    candidates: list[str] = []

    def add_col(name: str, values: np.ndarray) -> None:
        for row, value in zip(subjects, values):
            row[name] = float(value)
        candidates.append(name)

    for chset in CHSETS:
        for band in BANDS:
            frac_features = []
            for phase_fraction in PHASE_FRACTIONS:
                frac_tag = f"p{int(round(phase_fraction * 1000)):03d}"
                prefix = f"{chset}_{band}_{frac_tag}"
                c = np.array([ffloat(r.get(f"{prefix}__concentration")) for r in subjects])
                a = np.array([ffloat(r.get(f"{prefix}__local_align")) for r in subjects])
                g = np.array([ffloat(r.get(f"{prefix}__global_lock")) for r in subjects])
                b = np.array([ffloat(r.get(f"{prefix}__bint_ca")) for r in subjects])
                for metric, values in {
                    "concentration": c,
                    "local_align": a,
                    "global_lock": g,
                    "bint_ca": b,
                    "z_c_plus_a": zscore(c) + zscore(a),
                    "z_c_plus_a_minus_g": zscore(c) + zscore(a) - zscore(g),
                }.items():
                    name = f"{prefix}_{metric}"
                    add_col(name, values)
                    if metric in {"concentration", "bint_ca", "z_c_plus_a_minus_g"}:
                        frac_features.append(name)
            for metric in ("concentration", "bint_ca", "z_c_plus_a_minus_g"):
                cols = [
                    np.array([ffloat(r.get(f"{chset}_{band}_p{int(round(frac * 1000)):03d}_{metric}")) for r in subjects])
                    for frac in PHASE_FRACTIONS
                ]
                add_col(f"{chset}_{band}_phasefrac_{metric}_zmean", sum(zscore(c) for c in cols) / len(cols))
        # Narrow-band consensus across 4 one-Hz bins, separated from full theta.
        for phase_fraction in PHASE_FRACTIONS:
            frac_tag = f"p{int(round(phase_fraction * 1000)):03d}"
            for metric in ("concentration", "bint_ca", "z_c_plus_a_minus_g"):
                cols = [
                    np.array([ffloat(r.get(f"{chset}_{band}_{frac_tag}_{metric}")) for r in subjects])
                    for band in ("theta4_5", "theta5_6", "theta6_7", "theta7_8")
                ]
                add_col(f"{chset}_narrowtheta_{frac_tag}_{metric}_zmean", sum(zscore(c) for c in cols) / len(cols))
        for metric in ("concentration", "bint_ca", "z_c_plus_a_minus_g"):
            cols = [
                np.array([ffloat(r.get(f"{chset}_narrowtheta_p{int(round(frac * 1000)):03d}_{metric}_zmean")) for r in subjects])
                for frac in PHASE_FRACTIONS
            ]
            add_col(f"{chset}_narrowtheta_allphase_{metric}_zmean", sum(zscore(c) for c in cols) / len(cols))

    return candidates


def eval_feature(rows: list[dict[str, object]], feature: str, controls: list[str]) -> dict[str, object]:
    x = np.array([ffloat(r.get(feature)) for r in rows])
    y = np.array([ffloat(r.get("LPS")) for r in rows])
    wst = np.array([ffloat(r.get("WST")) for r in rows])
    ctrl = [np.array([ffloat(r.get(c)) for r in rows]) for c in controls]
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


def controls_for_feature(feature: str, analysis: str) -> list[str]:
    chset = "posterior" if feature.startswith("posterior_") else "all"
    band = next((b for b in BANDS if f"_{b}_" in feature), "theta4_8")
    local = [f"{chset}_{band}__band_power", f"{chset}_{band}__static_plv"]
    base = ["cohort_code", "age"] + local
    if analysis == "strict_both":
        other = "all" if chset == "posterior" else "posterior"
        return base + [f"{other}_{band}__band_power", f"{other}_{band}__static_plv"]
    return base


def evaluate_candidates(subjects: list[dict[str, object]], candidates: list[str]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for feature in candidates:
        for analysis in ("strict_local", "strict_both"):
            out.append({"analysis": analysis, **eval_feature(subjects, feature, controls_for_feature(feature, analysis))})
    stats.bh_q(out, "p_lps", "q_lps")
    stats.bh_q(out, "p_lps_plus_wst", "q_lps_plus_wst")
    return out


def split_stability(
    subjects: list[dict[str, object]],
    features: list[str],
    n_splits: int,
    seed: int,
) -> list[dict[str, object]]:
    if n_splits <= 0:
        return []
    rng = np.random.default_rng(seed)
    idx = np.arange(len(subjects))
    out: list[dict[str, object]] = []
    for feature in features:
        controls = controls_for_feature(feature, "strict_both")
        train_rs = []
        held_rs = []
        held_ps = []
        held_wst_rs = []
        for _ in range(n_splits):
            perm = rng.permutation(idx)
            half = len(perm) // 2
            train = [subjects[i] for i in perm[:half]]
            held = [subjects[i] for i in perm[half:]]
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
    raw_rows: list[dict[str, object]],
    subjects: list[dict[str, object]],
    eval_rows: list[dict[str, object]],
    split_rows: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    strict = [r for r in eval_rows if r["analysis"] == "strict_both"]
    top = sorted(
        strict,
        key=lambda r: ffloat(r["p_lps"]) if math.isfinite(ffloat(r["p_lps"])) else 1.0,
    )[:30]
    locked_names = {
        "posterior_theta4_8_p250_concentration",
        "posterior_theta4_8_p500_concentration",
        "posterior_theta4_8_phasefrac_concentration_zmean",
        "posterior_narrowtheta_p250_concentration_zmean",
        "posterior_narrowtheta_p500_concentration_zmean",
        "posterior_narrowtheta_allphase_concentration_zmean",
    }
    locked_rows = [r for r in strict if str(r["feature"]) in locked_names]

    def fmt(row: dict[str, object]) -> str:
        return (
            "- {feature}: n={n}, LPS r={r_lps:+.3f} p={p_lps:.4g} q={q_lps:.4g}, "
            "LPS+WST r={r_lps_plus_wst:+.3f} p={p_lps_plus_wst:.4g}, "
            "WST r={r_wst:+.3f} p={p_wst:.4g}"
        ).format(
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

    lines = [
        "# LEMON Frequency-Sliced Beam Evaluation",
        "",
        "Theory target: replace arbitrary millisecond lags with phase-relative wave slices.",
        "",
        f"Subjects: {len(subjects)}",
        f"Raw rows: {len(raw_rows)}",
        f"Cohorts: {dict(Counter(str(r['cohort']) for r in subjects))}",
        f"Bands: {', '.join(BANDS)}",
        f"Phase fractions: {', '.join(f'{x:.3f}' for x in PHASE_FRACTIONS)}",
        f"Window seconds: {args.window_sec}",
        f"Max windows: {args.max_windows}",
        "",
        "## Locked Rows",
        "",
    ]
    for row in sorted(locked_rows, key=lambda r: str(r["feature"])):
        lines.append(fmt(row))
    lines.extend(["", "## Best Strict-Both Rows", ""])
    for row in top:
        lines.append(fmt(row))
    lines.extend(["", "## Split Stability", ""])
    for row in split_rows:
        lines.append(
            "- {feature}: held r median={held_r_median:+.3f}, 5-95%=[{held_r_q05:+.3f},{held_r_q95:+.3f}], "
            "same direction={same_direction_rate:.1%}, held p<0.05={held_p_lt_0_05_rate:.1%}, held WST r={held_wst_r_median:+.3f}".format(
                **row
            )
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
    parser.add_argument("--tag", default="")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    suffix = "full" if not args.limit else f"limit{args.limit}"
    if args.tag:
        suffix = f"{suffix}_{args.tag}"
    OUT.mkdir(parents=True, exist_ok=True)
    raw_path = OUT / f"freq_sliced_raw_{suffix}.csv"
    subj_path = OUT / f"freq_sliced_subjects_{suffix}.csv"
    eval_path = OUT / f"freq_sliced_eval_{suffix}.csv"
    split_path = OUT / f"freq_sliced_splits_{suffix}.csv"
    summary_path = OUT / f"SUMMARY_{suffix}.md"

    if raw_path.exists() and not args.force:
        raw_rows = read_csv(raw_path)
        for row in raw_rows:
            for k, v in list(row.items()):
                if k not in {"sid", "cohort", "band", "channel_set"}:
                    row[k] = ffloat(v)
    else:
        manifest = load_manifest(args.limit)
        cfg = {
            "chidx": channel_sets(),
            "window_sec": args.window_sec,
            "step_sec": args.step_sec,
            "max_windows": args.max_windows,
        }
        print(f"Subjects queued: {len(manifest)}; workers={args.workers}; max_windows={args.max_windows}", flush=True)
        jobs = [(row, cfg) for row in manifest]
        raw_rows: list[dict[str, object]] = []
        if args.workers <= 1:
            for i, job in enumerate(jobs, start=1):
                print(f"  processing {i}/{len(jobs)} {job[0]['sid']}", flush=True)
                raw_rows.extend(process_subject(job))
        else:
            with Pool(args.workers) as pool:
                for i, rows in enumerate(pool.imap_unordered(process_subject, jobs), start=1):
                    raw_rows.extend(rows)
                    if i == 1 or i % 10 == 0:
                        print(f"  completed {i}/{len(jobs)} ({len(raw_rows)} rows)", flush=True)
        raw_rows.sort(
            key=lambda r: (
                str(r["sid"]),
                str(r["channel_set"]),
                str(r["band"]),
                float(r["phase_fraction"]),
            )
        )
        write_csv(raw_path, raw_rows)

    subjects = pivot_subjects(raw_rows)
    candidates = add_derived_candidates(subjects)
    write_csv(subj_path, subjects)
    eval_rows = evaluate_candidates(subjects, candidates)
    write_csv(eval_path, eval_rows)

    strict = [r for r in eval_rows if r["analysis"] == "strict_both"]
    top_features = [
        str(r["feature"])
        for r in sorted(strict, key=lambda r: ffloat(r["p_lps"]) if math.isfinite(ffloat(r["p_lps"])) else 1.0)[:8]
    ]
    locked = [
        "posterior_theta4_8_p250_concentration",
        "posterior_theta4_8_p500_concentration",
        "posterior_theta4_8_phasefrac_concentration_zmean",
        "posterior_narrowtheta_p250_concentration_zmean",
        "posterior_narrowtheta_p500_concentration_zmean",
        "posterior_narrowtheta_allphase_concentration_zmean",
    ]
    split_rows = split_stability(subjects, list(dict.fromkeys(locked + top_features)), args.splits, 20260627)
    write_csv(split_path, split_rows)
    write_summary(summary_path, raw_rows, subjects, eval_rows, split_rows, args)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
