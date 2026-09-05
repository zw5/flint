#!/usr/bin/env python3
"""LEMON beam-theory hypothesis dossier.

This script tests a small, theory-motivated dossier of beam hypotheses rather
than an open-ended EEG feature search.

New hypotheses tested here:

1. Aperture geometry: the beam footprint may be stronger in a more specific
   posterior aperture than the broad "not frontal" set.
2. Directional propagation: the beam may be a rectangular cross-field operator
   from one anatomical aperture at t to another at t + tau.
3. Nested validity: a small beam-theory candidate family should beat the
   baseline on held-out splits if the extra knobs are real physics.

It reuses the current winning core: broad-theta analytic phase and
phase-natural lags at 1/4, 1/3, and 1/2 cycle.
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
from lemon_frequency_sliced_beam_eval import (
    PHASE_FRACTIONS,
    concentration_from_spectrum,
    ffloat,
    lag_samples_for_phase_fraction,
    mean_adjacent_alignment,
    mean_pairwise_lock,
    static_plv,
    window_starts,
)


warnings.filterwarnings("ignore")

ROOT = Path("/Users/ximon/Documents/interp")
EXP = ROOT / "experiments/eeg_coherence"
EXPORT = EXP / "rust_export"
OUT = EXP / "lemon_beam_hypothesis_dossier"

THETA = (4.0, 8.0)


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


def load_channels() -> list[str]:
    return json.loads((EXPORT / "channels.json").read_text())["channels"]


def is_frontal(ch: str) -> bool:
    ch = ch.upper()
    return (
        ch.startswith("FP")
        or ch.startswith("AF")
        or (len(ch) >= 2 and ch[0] == "F" and ch[1] in "12345678Z")
    )


def side(ch: str) -> str:
    digits = "".join(c for c in ch if c.isdigit())
    if not digits:
        return "midline"
    return "left" if int(digits[-1]) % 2 else "right"


def build_apertures(channels: list[str]) -> dict[str, list[int]]:
    def select(pred) -> list[int]:
        return [i for i, ch in enumerate(channels) if pred(ch.upper())]

    posterior = [i for i, ch in enumerate(channels) if not is_frontal(ch)]
    parietal = select(lambda ch: ch.startswith("P") and not ch.startswith("PO"))
    parieto_occipital = select(lambda ch: ch.startswith("P") or ch.startswith("O"))
    occipital_po = select(lambda ch: ch.startswith("O") or ch.startswith("PO"))
    cp_parietal = select(lambda ch: ch.startswith("CP") or (ch.startswith("P") and not ch.startswith("PO")))
    lateral_posterior = [i for i in posterior if side(channels[i]) in {"left", "right"}]
    left_posterior = [i for i in posterior if side(channels[i]) == "left"]
    right_posterior = [i for i in posterior if side(channels[i]) == "right"]
    midline_posterior = [i for i in posterior if side(channels[i]) == "midline"]

    return {
        "posterior_broad": posterior,
        "parietal": parietal,
        "parieto_occipital": parieto_occipital,
        "occipital_po": occipital_po,
        "cp_parietal": cp_parietal,
        "lateral_posterior": lateral_posterior,
        "left_posterior": left_posterior,
        "right_posterior": right_posterior,
        "midline_posterior": midline_posterior,
    }


def build_directions(apertures: dict[str, list[int]]) -> dict[str, tuple[list[int], list[int]]]:
    return {
        "occipitalpo_to_parietal": (apertures["occipital_po"], apertures["parietal"]),
        "parietal_to_occipitalpo": (apertures["parietal"], apertures["occipital_po"]),
        "occipitalpo_to_cpparietal": (apertures["occipital_po"], apertures["cp_parietal"]),
        "cpparietal_to_occipitalpo": (apertures["cp_parietal"], apertures["occipital_po"]),
        "leftpost_to_rightpost": (apertures["left_posterior"], apertures["right_posterior"]),
        "rightpost_to_leftpost": (apertures["right_posterior"], apertures["left_posterior"]),
        "midlinepost_to_lateralpost": (apertures["midline_posterior"], apertures["lateral_posterior"]),
        "lateralpost_to_midlinepost": (apertures["lateral_posterior"], apertures["midline_posterior"]),
    }


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


def bandpass_theta(X: np.ndarray, sfreq: float) -> np.ndarray:
    lo, hi = THETA
    hi = min(hi, sfreq / 2.0 - 1.0)
    b, a = butter(4, [lo / (sfreq / 2.0), hi / (sfreq / 2.0)], btype="band")
    return filtfilt(b, a, X, axis=0)


def analytic_phase(X: np.ndarray) -> np.ndarray:
    analytic = hilbert(X, axis=0)
    return analytic / (np.abs(analytic) + 1e-9)


def cross_operator(future: np.ndarray, present: np.ndarray, lag: int) -> np.ndarray:
    X = present[:-lag]
    Y = future[lag:]
    return (Y.T @ X.conj()) / max(len(X), 1)


def leading_mode_features(future: np.ndarray, present: np.ndarray, lag: int) -> tuple[dict[str, float], np.ndarray]:
    C = cross_operator(future, present, lag)
    u, s, _vh = np.linalg.svd(C, full_matrices=False)
    concentration, effrank, lock = concentration_from_spectrum(s)
    mode = u[:, 0] / (np.linalg.norm(u[:, 0]) + 1e-12)
    return {"concentration": concentration, "effrank": effrank, "lock": lock}, mode


def phase_lag_features(
    future_phase: np.ndarray,
    present_phase: np.ndarray,
    lag: int,
    starts: list[int],
    win: int,
) -> dict[str, float]:
    metrics: dict[str, list[float]] = defaultdict(list)
    modes: list[np.ndarray] = []
    for start in starts:
        future = future_phase[start : start + win]
        present = present_phase[start : start + win]
        if len(future) <= lag + 10:
            continue
        try:
            observed, mode = leading_mode_features(future, present, lag)
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
    apertures: dict[str, list[int]] = cfg["apertures"]
    directions: dict[str, tuple[list[int], list[int]]] = cfg["directions"]
    win = int(round(float(cfg["window_sec"]) * sfreq))
    step = int(round(float(cfg["step_sec"]) * sfreq))
    starts = window_starts(len(X), win, step, int(cfg["max_windows"]))
    lags = {
        f"p{int(round(frac * 1000)):03d}": lag_samples_for_phase_fraction(sfreq, *THETA, frac)
        for frac in PHASE_FRACTIONS
    }

    try:
        Xb = bandpass_theta(X, sfreq)
        phase_all = analytic_phase(Xb)
    except Exception as exc:
        print(f"  {sid}: {exc}", flush=True)
        return []

    rows: list[dict[str, object]] = []

    def base_row(kind: str, name: str, present_idx: list[int], future_idx: list[int]) -> dict[str, object]:
        union = sorted(set(present_idx) | set(future_idx))
        return {
            "sid": sid,
            "LPS": float(row["LPS"]),
            "WST": float(row["WST"]),
            "age": float(row["age"]),
            "cohort": row["cohort"],
            "cohort_code": float(row["cohort_code"]),
            "kind": kind,
            "candidate": name,
            "present_n": len(present_idx),
            "future_n": len(future_idx),
            "union_n": len(union),
            "theta_power": float(np.mean(Xb[:, union] ** 2)),
            "theta_static_plv": static_plv(phase_all[:, union]),
        }

    for name, idx in apertures.items():
        if len(idx) < 3:
            continue
        for lag_name, lag in lags.items():
            feats = phase_lag_features(phase_all[:, idx], phase_all[:, idx], lag, starts, win)
            if not feats:
                continue
            rec = base_row("aperture", name, idx, idx)
            rec["phase_lag"] = lag_name
            rec.update(feats)
            rows.append(rec)

    for name, (src, dst) in directions.items():
        if len(src) < 3 or len(dst) < 3:
            continue
        for lag_name, lag in lags.items():
            feats = phase_lag_features(phase_all[:, dst], phase_all[:, src], lag, starts, win)
            if not feats:
                continue
            rec = base_row("direction", name, src, dst)
            rec["phase_lag"] = lag_name
            rec.update(feats)
            rows.append(rec)

    return rows


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
        prefix = f"{row['kind']}__{row['candidate']}__{row['phase_lag']}"
        for key, value in row.items():
            if key in {"sid", "LPS", "WST", "age", "cohort", "cohort_code", "kind", "candidate", "phase_lag"}:
                continue
            out[f"{prefix}__{key}"] = value
        cprefix = f"{row['kind']}__{row['candidate']}"
        out[f"{cprefix}__theta_power"] = row["theta_power"]
        out[f"{cprefix}__theta_static_plv"] = row["theta_static_plv"]
    return [by_sid[sid] for sid in sorted(by_sid)]


def candidate_prefix_from_feature(feature: str) -> str:
    if "__phasefrac_" in feature:
        return feature.split("__phasefrac_", 1)[0]
    for lag in ("__p250_", "__p333_", "__p500_"):
        if lag in feature:
            return feature.split(lag, 1)[0]
    return feature


def add_derived_candidates(subjects: list[dict[str, object]]) -> list[str]:
    prefixes = sorted(
        {
            key.removesuffix("__concentration").rsplit("__p", 1)[0]
            for row in subjects
            for key in row
            if key.endswith("__concentration") and "__p" in key
        }
    )
    candidates: list[str] = []

    def add_col(name: str, values: np.ndarray) -> None:
        for row, value in zip(subjects, values):
            row[name] = float(value)
        candidates.append(name)

    for prefix in prefixes:
        for lag in ("p250", "p333", "p500"):
            base = f"{prefix}__{lag}"
            c = np.array([ffloat(r.get(f"{base}__concentration")) for r in subjects])
            a = np.array([ffloat(r.get(f"{base}__local_align")) for r in subjects])
            g = np.array([ffloat(r.get(f"{base}__global_lock")) for r in subjects])
            b = np.array([ffloat(r.get(f"{base}__bint_ca")) for r in subjects])
            add_col(f"{prefix}__{lag}_concentration", c)
            add_col(f"{prefix}__{lag}_bint_ca", b)
            add_col(f"{prefix}__{lag}_z_c_plus_a_minus_g", zscore(c) + zscore(a) - zscore(g))
        for metric in ("concentration", "bint_ca", "z_c_plus_a_minus_g"):
            cols = [np.array([ffloat(r.get(f"{prefix}__{lag}_{metric}")) for r in subjects]) for lag in ("p250", "p333", "p500")]
            add_col(f"{prefix}__phasefrac_{metric}_zmean", sum(zscore(c) for c in cols) / len(cols))
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


def controls_for_feature(feature: str) -> list[str]:
    prefix = candidate_prefix_from_feature(feature)
    return ["cohort_code", "age", f"{prefix}__theta_power", f"{prefix}__theta_static_plv"]


def evaluate_candidates(subjects: list[dict[str, object]], candidates: list[str]) -> list[dict[str, object]]:
    out = [{"analysis": "strict", **eval_feature(subjects, feature, controls_for_feature(feature))} for feature in candidates]
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
        controls = controls_for_feature(feature)
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


def nested_selection(
    subjects: list[dict[str, object]],
    candidates: list[str],
    baseline: str,
    n_splits: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    idx = np.arange(len(subjects))
    rows: list[dict[str, object]] = []
    for split_i in range(n_splits):
        perm = rng.permutation(idx)
        half = len(perm) // 2
        train = [subjects[i] for i in perm[:half]]
        held = [subjects[i] for i in perm[half:]]
        scored = []
        for feature in candidates:
            ev = eval_feature(train, feature, controls_for_feature(feature))
            scored.append((abs(ffloat(ev["r_lps"])), feature, ev))
        scored.sort(reverse=True, key=lambda x: x[0] if math.isfinite(x[0]) else -1.0)
        selected = scored[0][1]
        held_ev = eval_feature(held, selected, controls_for_feature(selected))
        base_ev = eval_feature(held, baseline, controls_for_feature(baseline))
        rows.append(
            {
                "split": split_i,
                "selected": selected,
                "train_r": scored[0][2]["r_lps"],
                "held_r": held_ev["r_lps"],
                "held_p": held_ev["p_lps"],
                "held_wst_r": held_ev["r_wst"],
                "baseline_held_r": base_ev["r_lps"],
                "baseline_held_p": base_ev["p_lps"],
                "beats_baseline": float(ffloat(held_ev["r_lps"]) > ffloat(base_ev["r_lps"])),
            }
        )
    return rows


def summarize_nested(rows: list[dict[str, object]]) -> dict[str, object]:
    if not rows:
        return {}
    held = np.array([ffloat(r["held_r"]) for r in rows])
    base = np.array([ffloat(r["baseline_held_r"]) for r in rows])
    held_p = np.array([ffloat(r["held_p"]) for r in rows])
    wst = np.array([ffloat(r["held_wst_r"]) for r in rows])
    selected = Counter(str(r["selected"]) for r in rows)
    top_selected = "; ".join(f"{k}={v}" for k, v in selected.most_common(8))
    return {
        "n_splits": len(rows),
        "held_r_median": float(np.nanmedian(held)),
        "held_r_q05": float(np.nanquantile(held, 0.05)),
        "held_r_q95": float(np.nanquantile(held, 0.95)),
        "baseline_held_r_median": float(np.nanmedian(base)),
        "delta_vs_baseline_median": float(np.nanmedian(held - base)),
        "beats_baseline_rate": float(np.nanmean(held > base)),
        "held_p_lt_0_05_rate": float(np.nanmean(held_p < 0.05)),
        "held_wst_r_median": float(np.nanmedian(wst)),
        "top_selected": top_selected,
    }


def write_summary(
    path: Path,
    raw_rows: list[dict[str, object]],
    subjects: list[dict[str, object]],
    eval_rows: list[dict[str, object]],
    split_rows: list[dict[str, object]],
    nested_rows: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    strict = [r for r in eval_rows if r["analysis"] == "strict"]
    top = sorted(strict, key=lambda r: ffloat(r["p_lps"]) if math.isfinite(ffloat(r["p_lps"])) else 1.0)[:50]
    primary_names = {
        "aperture__posterior_broad__phasefrac_concentration_zmean",
        "aperture__parieto_occipital__phasefrac_concentration_zmean",
        "aperture__occipital_po__phasefrac_concentration_zmean",
        "aperture__parietal__phasefrac_concentration_zmean",
        "direction__occipitalpo_to_parietal__phasefrac_concentration_zmean",
        "direction__parietal_to_occipitalpo__phasefrac_concentration_zmean",
        "direction__leftpost_to_rightpost__phasefrac_concentration_zmean",
        "direction__rightpost_to_leftpost__phasefrac_concentration_zmean",
    }
    primary = [r for r in strict if str(r["feature"]) in primary_names]
    nested = summarize_nested(nested_rows)

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
        "# LEMON Beam-Theory Hypothesis Dossier",
        "",
        "Goal: test lawful beam-theory reframes without turning this into an arbitrary feature hunt.",
        "",
        f"Subjects: {len(subjects)}",
        f"Raw rows: {len(raw_rows)}",
        f"Cohorts: {dict(Counter(str(r['cohort']) for r in subjects))}",
        f"Window seconds: {args.window_sec}",
        f"Max windows: {args.max_windows}",
        "",
        "## Primary Hypothesis Rows",
        "",
    ]
    for row in sorted(primary, key=lambda r: str(r["feature"])):
        lines.append(fmt(row))

    lines.extend(["", "## Best Strict Rows", ""])
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

    lines.extend(["", "## Nested Candidate Test", ""])
    if nested:
        lines.extend(
            [
                f"- splits: {nested['n_splits']}",
                f"- selected held r median: {nested['held_r_median']:+.3f}",
                f"- selected held r 5-95%: [{nested['held_r_q05']:+.3f}, {nested['held_r_q95']:+.3f}]",
                f"- baseline held r median: {nested['baseline_held_r_median']:+.3f}",
                f"- median delta vs baseline: {nested['delta_vs_baseline_median']:+.3f}",
                f"- beats baseline rate: {nested['beats_baseline_rate']:.1%}",
                f"- selected held p<0.05: {nested['held_p_lt_0_05_rate']:.1%}",
                f"- held WST r median: {nested['held_wst_r_median']:+.3f}",
                f"- top selected: {nested['top_selected']}",
            ]
        )

    lines.extend(
        [
            "",
            "## Guardrail",
            "",
            "A larger in-sample row is not a smoking gun unless the nested held-out",
            "candidate test beats the fixed posterior baseline while preserving WST specificity.",
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
    parser.add_argument("--nested-splits", type=int, default=300)
    parser.add_argument("--tag", default="")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    suffix = "full" if not args.limit else f"limit{args.limit}"
    if args.tag:
        suffix = f"{suffix}_{args.tag}"
    OUT.mkdir(parents=True, exist_ok=True)
    raw_path = OUT / f"dossier_raw_{suffix}.csv"
    subj_path = OUT / f"dossier_subjects_{suffix}.csv"
    eval_path = OUT / f"dossier_eval_{suffix}.csv"
    split_path = OUT / f"dossier_splits_{suffix}.csv"
    nested_path = OUT / f"dossier_nested_{suffix}.csv"
    summary_path = OUT / f"SUMMARY_{suffix}.md"

    channels = load_channels()
    apertures = build_apertures(channels)
    directions = build_directions(apertures)

    if raw_path.exists() and not args.force:
        raw_rows = read_csv(raw_path)
        for row in raw_rows:
            for k, v in list(row.items()):
                if k not in {"sid", "cohort", "kind", "candidate", "phase_lag"}:
                    row[k] = ffloat(v)
    else:
        manifest = load_manifest(args.limit)
        cfg = {
            "apertures": apertures,
            "directions": directions,
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
        raw_rows.sort(key=lambda r: (str(r["sid"]), str(r["kind"]), str(r["candidate"]), str(r["phase_lag"])))
        write_csv(raw_path, raw_rows)

    subjects = pivot_subjects(raw_rows)
    candidates = add_derived_candidates(subjects)
    write_csv(subj_path, subjects)
    eval_rows = evaluate_candidates(subjects, candidates)
    write_csv(eval_path, eval_rows)

    baseline = "aperture__posterior_broad__phasefrac_concentration_zmean"
    strict = [r for r in eval_rows if r["analysis"] == "strict"]
    top_features = [
        str(r["feature"])
        for r in sorted(strict, key=lambda r: ffloat(r["p_lps"]) if math.isfinite(ffloat(r["p_lps"])) else 1.0)[:10]
    ]
    primary = [
        baseline,
        "aperture__parieto_occipital__phasefrac_concentration_zmean",
        "aperture__occipital_po__phasefrac_concentration_zmean",
        "aperture__parietal__phasefrac_concentration_zmean",
        "direction__occipitalpo_to_parietal__phasefrac_concentration_zmean",
        "direction__parietal_to_occipitalpo__phasefrac_concentration_zmean",
        "direction__leftpost_to_rightpost__phasefrac_concentration_zmean",
        "direction__rightpost_to_leftpost__phasefrac_concentration_zmean",
    ]
    split_rows = split_stability(subjects, list(dict.fromkeys(primary + top_features)), args.splits, 20260627)
    write_csv(split_path, split_rows)

    nested_candidates = [
        c
        for c in candidates
        if c.endswith("_concentration")
        or c.endswith("_concentration_zmean")
        or c.endswith("_bint_ca")
        or c.endswith("_bint_ca_zmean")
    ]
    nested_rows = nested_selection(subjects, nested_candidates, baseline, args.nested_splits, 20260627)
    write_csv(nested_path, nested_rows)

    write_summary(summary_path, raw_rows, subjects, eval_rows, split_rows, nested_rows, args)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
