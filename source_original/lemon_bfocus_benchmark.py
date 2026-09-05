#!/usr/bin/env python3
"""Locked LEMON predictive-focus benchmark.

This script tests the empirical beam-focus candidate discovered in the theta
predictive-beam pass:

    B_focus_ratio = lock_theta(tau) / effrank_theta(tau)
    B_focus_zdiff = z(lock_theta) - z(effrank_theta)

It deliberately benchmarks the score across conditions, frequency bands,
channel sets, and null families. The null families are diagnostic benchmarks;
the primary score is the observed lagged predictive operator.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
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
OUT = EXP / "lemon_bfocus_benchmark"

BANDS = {
    "delta": (1.0, 4.0),
    "theta": (4.0, 8.0),
    "alpha": (8.0, 13.0),
    "beta": (13.0, 30.0),
}
NULL_FAMILIES = ("sync_time_permute", "independent_circular_shift", "timewise_channel_permute")


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


def ffloat(x: object) -> float:
    try:
        out = float(x)
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


def load_channels() -> list[str]:
    return json.loads((EXPORT / "channels.json").read_text())["channels"]


def channel_sets(names: list[str]) -> dict[str, list[int]]:
    frontal = [i for i, ch in enumerate(names) if is_frontal(ch)]
    posterior = [i for i, ch in enumerate(names) if not is_frontal(ch)]
    return {
        "posterior": posterior,
        "frontal": frontal,
        "all": list(range(len(names))),
    }


def load_manifest(limit: int) -> list[dict[str, object]]:
    rows = []
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
        rec["n_ec"] = int(float(row["n_ec"]))
        rec["n_eo"] = int(float(row["n_eo"]))
        rec["cohort"] = c
        rec["cohort_code"] = 0.0 if c == "young" else 1.0
        if all(math.isfinite(float(rec[k])) for k in ("age", "LPS", "WST", "sfreq")):
            rows.append(rec)
    rows.sort(key=lambda r: str(r["sid"]))
    return rows[:limit] if limit else rows


def load_export(row: dict[str, object], condition: str) -> np.ndarray | None:
    sid = str(row["sid"])
    n_key = "n_eo" if condition == "EO" else "n_ec"
    n = int(row[n_key])
    if n <= 0:
        return None
    path = EXPORT / f"{sid}_{condition}.f32"
    if not path.exists():
        return None
    n_ch = int(row["n_ch"])
    data = np.fromfile(path, dtype=np.float32)
    expected = n_ch * n
    if data.size != expected:
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


def effective_rank(s: np.ndarray) -> float:
    e = np.maximum(s, 0.0) ** 2
    total = float(e.sum())
    if total <= 1e-12:
        return math.nan
    p = e / total
    p = p[p > 0]
    return float(np.exp(-(p * np.log(p)).sum()))


def lag_spectrum(Z: np.ndarray, lag: int) -> np.ndarray:
    present = Z[:-lag]
    future = Z[lag:]
    C = (future.T @ present.conj()) / max(len(present), 1)
    return np.linalg.svd(C, compute_uv=False)


def focus_from_spectrum(s: np.ndarray) -> tuple[float, float, float]:
    erank = effective_rank(s)
    lock = float((s[0] ** 2) / ((s**2).sum() + 1e-12))
    ratio = float(lock / (erank + 1e-12))
    return lock, erank, ratio


def window_starts(n: int, win: int, step: int, max_windows: int) -> list[int]:
    starts = list(range(0, max(n - win + 1, 0), step))
    if not starts:
        return []
    if max_windows and len(starts) > max_windows:
        idx = np.linspace(0, len(starts) - 1, max_windows).round().astype(int)
        starts = [starts[i] for i in idx]
    return starts


def null_versions(Z: np.ndarray, rng: np.random.Generator) -> dict[str, np.ndarray]:
    n, n_ch = Z.shape
    sync = Z[rng.permutation(n)]

    shifted = np.empty_like(Z)
    for j, shift in enumerate(rng.integers(0, n, size=n_ch)):
        shifted[:, j] = np.roll(Z[:, j], int(shift))

    perms = np.array([rng.permutation(n_ch) for _ in range(n)], dtype=int)
    timewise = Z[np.arange(n)[:, None], perms]

    return {
        "sync_time_permute": sync,
        "independent_circular_shift": shifted,
        "timewise_channel_permute": timewise,
    }


def beam_features_for_phase(
    phase: np.ndarray,
    sfreq: float,
    lags_ms: list[float],
    window_sec: float,
    step_sec: float,
    max_windows: int,
    n_null: int,
    seed: int,
) -> dict[str, float] | None:
    win = int(round(window_sec * sfreq))
    step = int(round(step_sec * sfreq))
    starts = window_starts(len(phase), win, step, max_windows)
    lags = [max(1, int(round(ms * sfreq / 1000.0))) for ms in lags_ms]
    if not starts:
        return None

    rng = np.random.default_rng(seed)
    locks = []
    ranks = []
    ratios = []
    null_ratios: dict[str, list[float]] = {name: [] for name in NULL_FAMILIES}

    for lag in lags:
        for start in starts:
            Z = phase[start : start + win]
            if len(Z) <= lag + 10:
                continue
            lock, erank, ratio = focus_from_spectrum(lag_spectrum(Z, lag))
            locks.append(lock)
            ranks.append(erank)
            ratios.append(ratio)
            for _ in range(n_null):
                for name, Zn in null_versions(Z, rng).items():
                    _, _, nr = focus_from_spectrum(lag_spectrum(Zn, lag))
                    null_ratios[name].append(nr)

    if not ratios:
        return None

    out: dict[str, float] = {
        "beam_lock": float(np.nanmean(locks)),
        "beam_effrank": float(np.nanmean(ranks)),
        "bfocus_ratio": float(np.nanmean(ratios)),
        "n_windows": float(len(starts)),
        "n_lag_windows": float(len(ratios)),
    }
    observed = out["bfocus_ratio"]
    for name, vals in null_ratios.items():
        vals_arr = np.asarray(vals, dtype=float)
        out[f"null_mean_{name}"] = float(np.nanmean(vals_arr))
        out[f"null_sd_{name}"] = float(np.nanstd(vals_arr))
        out[f"bfocus_nullz_{name}"] = float((observed - out[f"null_mean_{name}"]) / (out[f"null_sd_{name}"] + 1e-9))
    return out


def process_subject(job: tuple[dict[str, object], dict[str, object]]) -> list[dict[str, object]]:
    row, cfg = job
    conditions = cfg["conditions"]
    bands = cfg["bands"]
    chsets = cfg["chsets"]
    chidx = cfg["chidx"]
    lags_ms = cfg["lags_ms"]
    seed_base = cfg["seed_base"]

    sid = str(row["sid"])
    sfreq = float(row["sfreq"])
    out = []
    for condition in conditions:
        X = load_export(row, condition)
        if X is None:
            continue
        for band_name in bands:
            lo, hi = BANDS[band_name]
            try:
                Xb = bandpass(X, sfreq, lo, hi)
                analytic = hilbert(Xb, axis=0)
                phase_all = analytic / (np.abs(analytic) + 1e-9)
            except Exception as exc:
                print(f"  {sid} {condition} {band_name}: {exc}", flush=True)
                continue
            for chset in chsets:
                idx = chidx[chset]
                if len(idx) < 4:
                    continue
                phase = phase_all[:, idx]
                seed = seed_base + sum(ord(c) for c in f"{sid}:{condition}:{band_name}:{chset}")
                feats = beam_features_for_phase(
                    phase,
                    sfreq,
                    lags_ms,
                    float(cfg["window_sec"]),
                    float(cfg["step_sec"]),
                    int(cfg["max_windows"]),
                    int(cfg["nulls"]),
                    seed,
                )
                if feats is None:
                    continue
                Xsel = Xb[:, idx]
                rec = {
                    "sid": sid,
                    "LPS": float(row["LPS"]),
                    "WST": float(row["WST"]),
                    "age": float(row["age"]),
                    "cohort": row["cohort"],
                    "cohort_code": float(row["cohort_code"]),
                    "condition": condition,
                    "band": band_name,
                    "band_lo": lo,
                    "band_hi": hi,
                    "channel_set": chset,
                    "n_channels": len(idx),
                    "n_times": X.shape[0],
                    "band_power": float(np.mean(Xsel**2)),
                    "static_plv": static_plv(phase),
                }
                rec.update(feats)
                out.append(rec)
    return out


def add_config_scores(rows: list[dict[str, object]]) -> None:
    groups: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    for i, row in enumerate(rows):
        groups[(str(row["condition"]), str(row["band"]), str(row["channel_set"]))].append(i)
    for idxs in groups.values():
        lock = np.array([float(rows[i]["beam_lock"]) for i in idxs])
        erank = np.array([float(rows[i]["beam_effrank"]) for i in idxs])
        ratio = np.array([float(rows[i]["bfocus_ratio"]) for i in idxs])
        zdiff = zscore(lock) - zscore(erank)
        ratio_z = zscore(ratio)
        for local_i, row_i in enumerate(idxs):
            rows[row_i]["bfocus_zdiff"] = float(zdiff[local_i])
            rows[row_i]["bfocus_ratio_z"] = float(ratio_z[local_i])


def eval_feature(rows: list[dict[str, object]], feature: str, controls: list[str]) -> dict[str, object]:
    x = np.array([float(r[feature]) for r in rows])
    y = np.array([float(r["LPS"]) for r in rows])
    wst = np.array([float(r["WST"]) for r in rows])
    ctrl = [np.array([float(r[c]) for r in rows]) for c in controls]
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
    features = [
        "bfocus_ratio",
        "bfocus_zdiff",
        "bfocus_ratio_z",
        "beam_lock",
        "beam_effrank",
        "bfocus_nullz_sync_time_permute",
        "bfocus_nullz_independent_circular_shift",
        "bfocus_nullz_timewise_channel_permute",
    ]
    controls = {
        "cohort_age": ["cohort_code", "age"],
        "strict": ["cohort_code", "age", "band_power", "static_plv"],
    }
    out = []
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["condition"]), str(row["band"]), str(row["channel_set"]))].append(row)

    for (condition, band, chset), group in sorted(groups.items()):
        for analysis, ctrl in controls.items():
            for feature in features:
                out.append(
                    {
                        "condition": condition,
                        "band": band,
                        "channel_set": chset,
                        "analysis": analysis,
                        **eval_feature(group, feature, ctrl),
                    }
                )
        young = [r for r in group if r["cohort"] == "young"]
        for feature in features:
            out.append(
                {
                    "condition": condition,
                    "band": band,
                    "channel_set": chset,
                    "analysis": "young_age",
                    **eval_feature(young, feature, ["age"]),
                }
            )
            out.append(
                {
                    "condition": condition,
                    "band": band,
                    "channel_set": chset,
                    "analysis": "young_strict",
                    **eval_feature(young, feature, ["age", "band_power", "static_plv"]),
                }
            )

    stats.bh_q(out, "p_lps", "q_lps")
    stats.bh_q(out, "p_lps_plus_wst", "q_lps_plus_wst")
    return out


def split_stability(rows: list[dict[str, object]], n_splits: int, seed: int) -> list[dict[str, object]]:
    if n_splits <= 0:
        return []
    targets = [
        ("EO", "theta", "posterior", "bfocus_ratio", ["cohort_code", "age", "band_power", "static_plv"]),
        ("EO", "theta", "posterior", "bfocus_zdiff", ["cohort_code", "age", "band_power", "static_plv"]),
        ("EC", "theta", "posterior", "bfocus_ratio", ["cohort_code", "age", "band_power", "static_plv"]),
        ("EC", "theta", "posterior", "bfocus_zdiff", ["cohort_code", "age", "band_power", "static_plv"]),
    ]
    rng = np.random.default_rng(seed)
    out = []
    for condition, band, chset, feature, controls in targets:
        group = [
            r
            for r in rows
            if r["condition"] == condition and r["band"] == band and r["channel_set"] == chset
        ]
        if len(group) < 40:
            continue
        idx = np.arange(len(group))
        train_rs = []
        held_rs = []
        train_ps = []
        held_ps = []
        held_wst_rs = []
        for _ in range(n_splits):
            perm = rng.permutation(idx)
            half = len(perm) // 2
            train = [group[i] for i in perm[:half]]
            held = [group[i] for i in perm[half:]]
            tr = eval_feature(train, feature, controls)
            he = eval_feature(held, feature, controls)
            train_rs.append(float(tr["r_lps"]))
            held_rs.append(float(he["r_lps"]))
            train_ps.append(float(tr["p_lps"]))
            held_ps.append(float(he["p_lps"]))
            held_wst_rs.append(float(he["r_lps_plus_wst"]))
        train_rs_arr = np.asarray(train_rs)
        held_rs_arr = np.asarray(held_rs)
        train_ps_arr = np.asarray(train_ps)
        held_ps_arr = np.asarray(held_ps)
        out.append(
            {
                "condition": condition,
                "band": band,
                "channel_set": chset,
                "feature": feature,
                "analysis": "strict_split",
                "n": len(group),
                "held_r_median": float(np.nanmedian(held_rs_arr)),
                "held_r_q05": float(np.nanquantile(held_rs_arr, 0.05)),
                "held_r_q95": float(np.nanquantile(held_rs_arr, 0.95)),
                "same_direction_rate": float(np.nanmean(np.sign(train_rs_arr) == np.sign(held_rs_arr))),
                "held_p_lt_0_05_rate": float(np.nanmean(held_ps_arr < 0.05)),
                "both_p_lt_0_05_rate": float(np.nanmean((held_ps_arr < 0.05) & (train_ps_arr < 0.05))),
                "held_lps_plus_wst_r_median": float(np.nanmedian(held_wst_rs)),
                "train_r_median": float(np.nanmedian(train_rs_arr)),
            }
        )
    return out


def parse_csv_arg(value: str, valid: set[str]) -> list[str]:
    items = [x.strip() for x in value.split(",") if x.strip()]
    bad = [x for x in items if x not in valid]
    if bad:
        raise SystemExit(f"Unknown values {bad}; valid: {sorted(valid)}")
    return items


def write_summary(
    path: Path,
    rows: list[dict[str, object]],
    eval_rows: list[dict[str, object]],
    split_rows: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    def display_row(row: dict[str, object]) -> dict[str, object]:
        out = dict(row)
        for key in ("r_lps", "p_lps", "q_lps", "r_lps_plus_wst", "p_lps_plus_wst", "r_wst", "p_wst"):
            value = ffloat(out.get(key))
            out[key] = value if math.isfinite(value) else math.nan
        return out

    primary = [
        r
        for r in eval_rows
        if r["condition"] == "EO"
        and r["band"] == "theta"
        and r["channel_set"] == "posterior"
        and r["feature"] in {"bfocus_ratio", "bfocus_zdiff"}
    ]
    strict_primary = [r for r in primary if r["analysis"] in {"strict", "young_strict"}]
    top_strict = sorted(
        [r for r in eval_rows if r["analysis"] == "strict" and r["feature"] in {"bfocus_ratio", "bfocus_zdiff"}],
        key=lambda r: float(r["p_lps"]) if math.isfinite(float(r["p_lps"])) else 1.0,
    )[:16]
    null_top = sorted(
        [r for r in eval_rows if r["analysis"] == "strict" and str(r["feature"]).startswith("bfocus_nullz")],
        key=lambda r: float(r["p_lps"]) if math.isfinite(float(r["p_lps"])) else 1.0,
    )[:12]

    lines = [
        "# LEMON B_focus Benchmark",
        "",
        f"Usable feature rows: {len(rows)}",
        f"Subjects per cohort: {dict(Counter(str(r['cohort']) for r in rows if r['condition'] == args.conditions.split(',')[0] and r['band'] == args.bands.split(',')[0] and r['channel_set'] == args.channel_sets.split(',')[0]))}",
        f"Conditions: {args.conditions}",
        f"Bands: {args.bands}",
        f"Channel sets: {args.channel_sets}",
        f"Lags ms: {args.lags_ms}",
        f"Window seconds: {args.window_sec}",
        f"Max windows: {args.max_windows}",
        f"Nulls per window-lag-family: {args.nulls}",
        "",
        "## Locked Primary Rows",
        "",
    ]
    for r in strict_primary:
        lines.append(
            "- {analysis} / {feature}: n={n}, LPS r={r_lps:+.3f} p={p_lps:.4g} q={q_lps:.4g}, "
            "LPS+WST r={r_lps_plus_wst:+.3f} p={p_lps_plus_wst:.4g}, WST r={r_wst:+.3f} p={p_wst:.4g}".format(
                **display_row(r)
            )
        )
    lines.extend(["", "## Best Strict B_focus Benchmarks", ""])
    for r in top_strict:
        lines.append(
            "- {condition}/{band}/{channel_set} / {feature}: n={n}, LPS r={r_lps:+.3f} p={p_lps:.4g} q={q_lps:.4g}, "
            "LPS+WST r={r_lps_plus_wst:+.3f} p={p_lps_plus_wst:.4g}, WST r={r_wst:+.3f} p={p_wst:.4g}".format(
                **display_row(r)
            )
        )
    lines.extend(["", "## Null-Normalized Strict Benchmarks", ""])
    for r in null_top:
        lines.append(
            "- {condition}/{band}/{channel_set} / {feature}: n={n}, LPS r={r_lps:+.3f} p={p_lps:.4g} q={q_lps:.4g}, "
            "LPS+WST r={r_lps_plus_wst:+.3f} p={p_lps_plus_wst:.4g}, WST r={r_wst:+.3f} p={p_wst:.4g}".format(
                **display_row(r)
            )
        )
    lines.extend(["", "## Split-Half Stability", ""])
    for r in split_rows:
        lines.append(
            "- {condition}/{band}/{channel_set} / {feature}: held r median {held_r_median:+.3f} "
            "[{held_r_q05:+.3f}, {held_r_q95:+.3f}], same direction {same_direction_rate:.1%}, "
            "held p<0.05 {held_p_lt_0_05_rate:.1%}, both p<0.05 {both_p_lt_0_05_rate:.1%}".format(**r)
        )
    lines.extend(
        [
            "",
            "## Guardrails",
            "",
            "- Primary scores are observed lagged predictive operator scores, not fit to LPS.",
            "- B_focus_ratio is lock / effective-rank.",
            "- B_focus_zdiff is z(lock) - z(effective-rank), standardized within each condition/band/channel-set.",
            "- Null-normalized rows compare observed B_focus_ratio with null B_focus_ratio for three null families.",
            "- This is still sensor-space predictive beamforming, not source-space lead-field beamforming.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=max(1, min(6, cpu_count() - 2)))
    parser.add_argument("--conditions", default="EO,EC")
    parser.add_argument("--bands", default="theta,delta,alpha,beta")
    parser.add_argument("--channel-sets", default="posterior,frontal,all")
    parser.add_argument("--window-sec", type=float, default=8.0)
    parser.add_argument("--step-sec", type=float, default=4.0)
    parser.add_argument("--lags-ms", default="40,80,120,160")
    parser.add_argument("--max-windows", type=int, default=16)
    parser.add_argument("--nulls", type=int, default=8)
    parser.add_argument("--splits", type=int, default=300)
    parser.add_argument("--tag", default="")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    conditions = parse_csv_arg(args.conditions, {"EO", "EC"})
    bands = parse_csv_arg(args.bands, set(BANDS))
    names = load_channels()
    chidx = channel_sets(names)
    chsets = parse_csv_arg(args.channel_sets, set(chidx))
    lags_ms = [float(x.strip()) for x in args.lags_ms.split(",") if x.strip()]

    OUT.mkdir(parents=True, exist_ok=True)
    suffix = "full" if not args.limit else f"limit{args.limit}"
    if args.tag:
        suffix = f"{suffix}_{args.tag}"
    feature_path = OUT / f"bfocus_features_{suffix}.csv"
    eval_path = OUT / f"bfocus_eval_{suffix}.csv"
    split_path = OUT / f"bfocus_splits_{suffix}.csv"
    summary_path = OUT / f"SUMMARY_{suffix}.md"

    manifest = load_manifest(args.limit)
    if feature_path.exists() and not args.force:
        rows = read_csv(feature_path)
        for row in rows:
            for key in (
                "LPS",
                "WST",
                "age",
                "cohort_code",
                "band_lo",
                "band_hi",
                "n_channels",
                "n_times",
                "band_power",
                "static_plv",
                "beam_lock",
                "beam_effrank",
                "bfocus_ratio",
                "n_windows",
                "n_lag_windows",
                "bfocus_zdiff",
                "bfocus_ratio_z",
                *[f"null_mean_{n}" for n in NULL_FAMILIES],
                *[f"null_sd_{n}" for n in NULL_FAMILIES],
                *[f"bfocus_nullz_{n}" for n in NULL_FAMILIES],
            ):
                if key in row:
                    row[key] = ffloat(row[key])
    else:
        cfg = {
            "conditions": conditions,
            "bands": bands,
            "chsets": chsets,
            "chidx": chidx,
            "lags_ms": lags_ms,
            "window_sec": args.window_sec,
            "step_sec": args.step_sec,
            "max_windows": args.max_windows,
            "nulls": args.nulls,
            "seed_base": 20260627,
        }
        print(
            f"Subjects queued: {len(manifest)}; configs={len(conditions) * len(bands) * len(chsets)}; "
            f"workers={args.workers}; nulls={args.nulls}; max_windows={args.max_windows}",
            flush=True,
        )
        jobs = [(row, cfg) for row in manifest]
        rows = []
        if args.workers <= 1:
            for i, job in enumerate(jobs, start=1):
                print(f"  processing {i}/{len(jobs)} {job[0]['sid']}", flush=True)
                rows.extend(process_subject(job))
        else:
            with Pool(args.workers) as pool:
                for i, subject_rows in enumerate(pool.imap_unordered(process_subject, jobs), start=1):
                    rows.extend(subject_rows)
                    if i == 1 or i % 10 == 0:
                        print(f"  completed {i}/{len(jobs)} ({len(rows)} feature rows)", flush=True)
        rows.sort(key=lambda r: (str(r["sid"]), str(r["condition"]), str(r["band"]), str(r["channel_set"])))
        add_config_scores(rows)
        fields = [
            "sid",
            "LPS",
            "WST",
            "age",
            "cohort",
            "cohort_code",
            "condition",
            "band",
            "band_lo",
            "band_hi",
            "channel_set",
            "n_channels",
            "n_times",
            "band_power",
            "static_plv",
            "beam_lock",
            "beam_effrank",
            "bfocus_ratio",
            "bfocus_zdiff",
            "bfocus_ratio_z",
            "n_windows",
            "n_lag_windows",
            *[f"null_mean_{n}" for n in NULL_FAMILIES],
            *[f"null_sd_{n}" for n in NULL_FAMILIES],
            *[f"bfocus_nullz_{n}" for n in NULL_FAMILIES],
        ]
        write_csv(feature_path, rows, fields)

    if rows and "bfocus_zdiff" not in rows[0]:
        add_config_scores(rows)
    eval_rows = evaluate(rows)
    eval_fields = [
        "condition",
        "band",
        "channel_set",
        "analysis",
        "feature",
        "n",
        "controls",
        "r_lps",
        "p_lps",
        "q_lps",
        "r_lps_plus_wst",
        "p_lps_plus_wst",
        "q_lps_plus_wst",
        "r_wst",
        "p_wst",
    ]
    write_csv(eval_path, eval_rows, eval_fields)
    split_rows = split_stability(rows, args.splits, 20260627)
    if split_rows:
        write_csv(split_path, split_rows)
    write_summary(summary_path, rows, eval_rows, split_rows, args)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
