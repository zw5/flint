#!/usr/bin/env python3
"""LEMON spatial-filter ceiling pass.

This is a conservative source-proxy test, not a claim of source localization.
It asks whether the carrier-aligned beam result improves when the sensor field
is transformed with fixed, outcome-blind spatial filters derived from standard
montage geometry.

Locked transforms:

    raw       no spatial transform
    smooth35  Gaussian spatial smoothing, sigma=.35 chord units
    smooth50  Gaussian spatial smoothing, sigma=.50 chord units
    deblur10  mild high-boost: X + .10 * (X - smooth35(X))
    deblur25  mild high-boost: X + .25 * (X - smooth35(X))
    deblur50  mild high-boost: X + .50 * (X - smooth35(X))
    deblur75  stronger high-boost: X + .75 * (X - smooth35(X))

The goal is to test whether preserving/deblurring the broad field mode improves
the same posterior/lateral theta self-propagation object.
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from collections import Counter
from multiprocessing import Pool, cpu_count
from pathlib import Path

import mne
import numpy as np

import lemon_intelligence_residual_search as stats
import lemon_beam_reliability_carrier as carrier
from lemon_beam_hypothesis_dossier import build_apertures
from lemon_frequency_sliced_beam_eval import ffloat, static_plv, window_starts


warnings.filterwarnings("ignore")
mne.set_log_level("ERROR")

ROOT = Path("/Users/ximon/Documents/interp")
EXP = ROOT / "experiments/eeg_coherence"
OUT = EXP / "lemon_spatial_filter_ceiling"

TRANSFORMS = ("raw", "smooth35", "smooth50", "deblur10", "deblur25", "deblur50", "deblur75")
APERTURES = ("posterior_broad", "lateral_posterior")


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


def channel_positions(channels: list[str]) -> np.ndarray:
    montage = mne.channels.make_standard_montage("standard_1020")
    pos_map = montage.get_positions()["ch_pos"]
    rows = []
    missing = []
    for ch in channels:
        key = ch
        if key not in pos_map and ch.upper() in pos_map:
            key = ch.upper()
        if key not in pos_map:
            missing.append(ch)
            continue
        rows.append(pos_map[key])
    if missing:
        raise ValueError(f"Missing standard montage positions: {missing}")
    pos = np.asarray(rows, dtype=float)
    pos = pos / (np.linalg.norm(pos, axis=1, keepdims=True) + 1e-12)
    return pos


def spatial_smoother(pos: np.ndarray, sigma: float) -> np.ndarray:
    pos = np.asarray(pos, dtype=float)
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=2)
    W = np.exp(-(d**2) / (2.0 * sigma * sigma))
    W = W / (W.sum(axis=1, keepdims=True) + 1e-12)
    return W


def apply_spatial_transform(X: np.ndarray, smoother: np.ndarray, kind: str, amount: float = 1.0) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    smooth = X @ smoother.T
    if kind == "raw":
        return X
    if kind == "smooth":
        return smooth
    if kind == "deblur":
        return X + amount * (X - smooth)
    raise ValueError(f"Unknown spatial transform: {kind}")


def transform_matrix(X: np.ndarray, transform: str, smooth35: np.ndarray, smooth50: np.ndarray) -> np.ndarray:
    if transform == "raw":
        return X
    if transform == "smooth35":
        return apply_spatial_transform(X, smooth35, "smooth")
    if transform == "smooth50":
        return apply_spatial_transform(X, smooth50, "smooth")
    if transform == "deblur10":
        return apply_spatial_transform(X, smooth35, "deblur", amount=0.10)
    if transform == "deblur25":
        return apply_spatial_transform(X, smooth35, "deblur", amount=0.25)
    if transform == "deblur50":
        return apply_spatial_transform(X, smooth35, "deblur", amount=0.50)
    if transform == "deblur75":
        return apply_spatial_transform(X, smooth35, "deblur", amount=0.75)
    raise ValueError(transform)


def process_subject(job: tuple[dict[str, object], dict[str, object]]) -> dict[str, object] | None:
    row, cfg = job
    sid = str(row["sid"])
    X = carrier.load_eo(row)
    if X is None:
        return None
    sfreq = float(row["sfreq"])
    channels: list[str] = cfg["channels"]
    apertures = build_apertures(channels)
    smooth35 = cfg["smooth35"]
    smooth50 = cfg["smooth50"]
    win = int(round(float(cfg["window_sec"]) * sfreq))
    step = int(round(float(cfg["step_sec"]) * sfreq))
    starts = window_starts(len(X), win, step, int(cfg["max_windows"]))
    start_sets = carrier.split_start_sets(starts)
    out: dict[str, object] = {
        "sid": sid,
        "LPS": float(row["LPS"]),
        "WST": float(row["WST"]),
        "age": float(row["age"]),
        "cohort": row["cohort"],
        "cohort_code": float(row["cohort_code"]),
    }

    for transform in TRANSFORMS:
        try:
            Xt = transform_matrix(X, transform, smooth35, smooth50)
            Xb = carrier.bandpass_theta(Xt, sfreq)
            phase_all = carrier.analytic_phase(Xb)
        except Exception as exc:
            print(f"  {sid} {transform}: {exc}", flush=True)
            continue
        for aperture in APERTURES:
            idx = apertures[aperture]
            if len(idx) < 3:
                continue
            Xsel = Xb[:, idx]
            phase = phase_all[:, idx]
            peak = carrier.theta_peak_frequency(Xsel, sfreq, *carrier.THETA)
            lags = {
                f"p{int(round(frac * 1000)):03d}": carrier.phase_lag_samples(sfreq, peak, frac)
                for frac in carrier.PHASE_FRACTIONS
            }
            prefix = f"{transform}_{aperture}"
            out[f"{prefix}__theta_peak_hz"] = peak
            out[f"{prefix}__theta_power"] = float(np.mean(Xsel**2))
            out[f"{prefix}__theta_static_plv"] = static_plv(phase)
            for split_name, split_starts in start_sets.items():
                scores = carrier.score_for_starts(phase, lags, split_starts, win)
                for key, value in scores.items():
                    out[f"{prefix}__{split_name}_{key}"] = value
    return out


def controls_for(feature: str, include_peak: bool = False) -> list[str]:
    prefix = feature.split("__", 1)[0]
    controls = [
        "cohort_code",
        "age",
        f"{prefix}__theta_power",
        f"{prefix}__theta_static_plv",
    ]
    if include_peak:
        controls.append(f"{prefix}__theta_peak_hz")
    return controls


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


def evaluate(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for transform in TRANSFORMS:
        for aperture in APERTURES:
            prefix = f"{transform}_{aperture}"
            features = [f"{prefix}__full_score"] + [
                f"{prefix}__full_p{int(round(frac * 1000)):03d}_concentration"
                for frac in carrier.PHASE_FRACTIONS
            ]
            for include_peak in (False, True):
                analysis = "strict_peak" if include_peak else "strict"
                for feature in features:
                    out.append({"analysis": analysis, **eval_feature(rows, feature, controls_for(feature, include_peak))})
    stats.bh_q(out, "p_lps", "q_lps")
    return out


def reliability_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for transform in TRANSFORMS:
        for aperture in APERTURES:
            prefix = f"{transform}_{aperture}"
            odd = np.asarray([ffloat(r.get(f"{prefix}__odd_score")) for r in rows])
            even = np.asarray([ffloat(r.get(f"{prefix}__even_score")) for r in rows])
            first = np.asarray([ffloat(r.get(f"{prefix}__first_score")) for r in rows])
            second = np.asarray([ffloat(r.get(f"{prefix}__second_score")) for r in rows])
            for split, a, b in (("odd_even", odd, even), ("first_second", first, second)):
                mask = np.isfinite(a) & np.isfinite(b)
                r = float(np.corrcoef(stats.rankdata(a[mask]), stats.rankdata(b[mask]))[0, 1]) if int(mask.sum()) >= 5 else math.nan
                out.append(
                    {
                        "transform": transform,
                        "aperture": aperture,
                        "split": split,
                        "n": int(mask.sum()),
                        "split_r": r,
                        "spearman_brown": carrier.spearman_brown(r),
                    }
                )
    return out


def split_stability(rows: list[dict[str, object]], features: list[str], n_splits: int, seed: int) -> list[dict[str, object]]:
    if n_splits <= 0:
        return []
    rng = np.random.default_rng(seed)
    idx = np.arange(len(rows))
    out: list[dict[str, object]] = []
    for feature in features:
        controls = controls_for(feature, include_peak=False)
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

    primary_features = {f"{transform}_{aperture}__full_score" for transform in TRANSFORMS for aperture in APERTURES}
    primary = [r for r in eval_rows if r["feature"] in primary_features]
    top = sorted(eval_rows, key=lambda r: ffloat(r["p_lps"]) if math.isfinite(ffloat(r["p_lps"])) else 1.0)[:30]
    lines = [
        "# LEMON Spatial-Filter Ceiling",
        "",
        "Goal: test whether fixed geometry-based smoothing/deblurring improves the carrier-aligned beam object.",
        "",
        f"Subjects: {len(subjects)}",
        f"Cohorts: {dict(Counter(str(r['cohort']) for r in subjects))}",
        f"Transforms: {', '.join(TRANSFORMS)}",
        f"Window seconds: {args.window_sec}",
        f"Max windows: {args.max_windows}",
        "",
        "## Primary Rows",
        "",
    ]
    for row in sorted(primary, key=lambda r: (str(r["feature"]), str(r["analysis"]))):
        lines.append(fmt(row))
    lines.extend(["", "## Reliability", ""])
    for row in rel_rows:
        lines.append(
            "- {transform} {aperture} {split}: n={n}, split r={split_r:+.3f}, Spearman-Brown={spearman_brown:+.3f}".format(
                transform=row["transform"],
                aperture=row["aperture"],
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
            "Spatial transforms are promoted only if they beat the raw carrier-aligned posterior score in held-out stability while preserving near-zero WST association. Smoothing/deblurring here is a sensor-space source proxy, not source localization.",
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
    subj_path = OUT / f"spatial_subjects_{suffix}.csv"
    eval_path = OUT / f"spatial_eval_{suffix}.csv"
    rel_path = OUT / f"spatial_reliability_{suffix}.csv"
    split_path = OUT / f"spatial_splits_{suffix}.csv"
    summary_path = OUT / f"SUMMARY_{suffix}.md"

    if subj_path.exists() and not args.force:
        subjects = read_csv(subj_path)
        for row in subjects:
            for key, value in list(row.items()):
                if key not in {"sid", "cohort"}:
                    row[key] = ffloat(value)
    else:
        manifest = carrier.load_manifest(args.limit)
        channels = carrier.load_channels()
        pos = channel_positions(channels)
        cfg = {
            "channels": channels,
            "smooth35": spatial_smoother(pos, 0.35),
            "smooth50": spatial_smoother(pos, 0.50),
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
    primary = [f"{transform}_{aperture}__full_score" for transform in TRANSFORMS for aperture in APERTURES]
    extras = [
        "raw_posterior_broad__full_p500_concentration",
        "smooth35_posterior_broad__full_p500_concentration",
        "deblur10_posterior_broad__full_p250_concentration",
        "deblur10_posterior_broad__full_p333_concentration",
        "deblur10_posterior_broad__full_score",
        "deblur25_posterior_broad__full_p250_concentration",
        "deblur25_posterior_broad__full_p333_concentration",
        "deblur25_posterior_broad__full_p500_concentration",
        "deblur50_lateral_posterior__full_p250_concentration",
        "deblur50_lateral_posterior__full_p333_concentration",
        "deblur75_posterior_broad__full_p250_concentration",
        "deblur75_posterior_broad__full_p333_concentration",
        "deblur75_lateral_posterior__full_score",
    ]
    split_rows = split_stability(subjects, list(dict.fromkeys(primary + extras)), args.splits, 20260627)
    write_csv(eval_path, eval_rows)
    write_csv(rel_path, rel_rows)
    write_csv(split_path, split_rows)
    write_summary(summary_path, subjects, eval_rows, rel_rows, split_rows, args)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
