#!/usr/bin/env python3
"""LEMON intelligence residual correlate search.

This is the revised-theory LEMON pipeline: search for EEG dynamics that track
LPS fluid reasoning while testing the obvious failure modes:

  - LEMON young/old cohort structure
  - age-bin midpoint
  - WST vocabulary specificity
  - power/coherence nuisance axes
  - young-cohort-only survival

It uses existing extracted feature tables, so it does not recompute EEG.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, t as tdist


ROOT = Path("/Users/ximon/Documents/interp")
EXP = ROOT / "experiments/eeg_coherence"
OUT = EXP / "lemon_intelligence_residual_search"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
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


def cohort_name(age: float) -> str:
    if age <= 35:
        return "young"
    if age >= 55:
        return "old"
    return "middle"


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return np.nan_to_num((x - np.nanmean(x)) / (np.nanstd(x) + 1e-9))


def rank_residual(x: np.ndarray, controls: list[np.ndarray]) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    mask = np.isfinite(x)
    for c in controls:
        mask &= np.isfinite(c)
    out = np.full_like(x, np.nan, dtype=float)
    if int(mask.sum()) <= len(controls) + 3:
        return out
    xr = rankdata(x[mask])
    if not controls:
        out[mask] = xr - xr.mean()
        return out
    C = np.column_stack([rankdata(np.asarray(c)[mask]) for c in controls])
    B = np.column_stack([np.ones(len(xr)), C])
    coef, *_ = np.linalg.lstsq(B, xr, rcond=None)
    out[mask] = xr - B @ coef
    return out


def partial_spearman(x: np.ndarray, y: np.ndarray, controls: list[np.ndarray]) -> tuple[float, float, int]:
    rx = rank_residual(x, controls)
    ry = rank_residual(y, controls)
    mask = np.isfinite(rx) & np.isfinite(ry)
    n = int(mask.sum())
    if n <= len(controls) + 3:
        return math.nan, math.nan, n
    rx = rx[mask]
    ry = ry[mask]
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = n - len(controls) - 2
    p = float(2 * tdist.sf(abs(r * np.sqrt(dof / max(1 - r * r, 1e-9))), dof))
    return r, p, n


def perm_p(x: np.ndarray, y: np.ndarray, controls: list[np.ndarray], observed_abs_r: float, n_perm: int, seed: int) -> float:
    if n_perm <= 0 or not math.isfinite(observed_abs_r):
        return math.nan
    rng = np.random.default_rng(seed)
    ry = rank_residual(y, controls)
    rx = rank_residual(x, controls)
    mask = np.isfinite(rx) & np.isfinite(ry)
    if int(mask.sum()) < 10:
        return math.nan
    rx = rx[mask]
    ry = ry[mask]
    ge = 1
    for _ in range(n_perm):
        rp = float(np.corrcoef(rx, rng.permutation(ry))[0, 1])
        if abs(rp) >= observed_abs_r:
            ge += 1
    return float(ge / (n_perm + 1))


def bh_q(rows: list[dict[str, object]], p_key: str, q_key: str) -> None:
    finite = [(i, float(r[p_key])) for i, r in enumerate(rows) if str(r.get(p_key, "")) != "" and math.isfinite(float(r[p_key]))]
    finite.sort(key=lambda x: x[1])
    m = len(finite)
    prev = 1.0
    qvals = [1.0] * m
    for rank in range(m, 0, -1):
        idx, p = finite[rank - 1]
        prev = min(prev, p * m / rank)
        qvals[rank - 1] = prev
    for (idx, _), q in zip(finite, qvals):
        rows[idx][q_key] = q
    for r in rows:
        r.setdefault(q_key, "")


def add_static_features(subjects: dict[str, dict[str, object]]) -> None:
    for row in read_csv(EXP / "features.csv"):
        sid = row.get("subject")
        if not sid:
            continue
        rec = subjects.setdefault(sid, {})
        for k, v in row.items():
            if k in {"subject", "n_ch"}:
                continue
            if k in {"WST", "LPS", "age"}:
                rec[k] = ffloat(v)
            else:
                rec[f"static__{k}"] = ffloat(v)


def add_trajectory_features(subjects: dict[str, dict[str, object]]) -> None:
    for row in read_csv(EXP / "trajectory_features.csv"):
        sid = row.get("sid")
        if not sid:
            continue
        rec = subjects.setdefault(sid, {})
        for k, v in row.items():
            if k == "sid":
                continue
            if k in {"WST", "LPS", "age"}:
                rec[k] = ffloat(v)
            else:
                rec[f"traj__{k}"] = ffloat(v)


def add_compressed_features(subjects: dict[str, dict[str, object]]) -> None:
    for row in read_csv(EXP / "compressed_sensing/results.csv"):
        sid = row.get("sid")
        if not sid:
            continue
        rec = subjects.setdefault(sid, {})
        for k, v in row.items():
            if k == "sid":
                continue
            if k in {"WST", "LPS", "age"}:
                rec[k] = ffloat(v)
            else:
                rec[f"cs__{k}"] = ffloat(v)


def add_dynamic_features(subjects: dict[str, dict[str, object]]) -> None:
    for row in read_csv(EXP / "dynamic_features.csv"):
        sid = row.get("sid")
        band = row.get("band")
        cond = row.get("cond")
        if not sid or not band or not cond:
            continue
        rec = subjects.setdefault(sid, {})
        for k, v in row.items():
            if k in {"sid", "band", "cond"}:
                continue
            if k in {"WST", "LPS", "age"}:
                rec[k] = ffloat(v)
                continue
            val = ffloat(v)
            rec[f"dyn__{cond}__{band}__{k}"] = val
        for metric in (
            "h_mean",
            "h_range",
            "h_std",
            "beam_depth",
            "beam_max",
            "meta_r",
            "r_mean",
            "r_range",
            "traj_dim",
            "switch",
        ):
            obs = ffloat(row.get(f"{metric}_obs"))
            sn_mu = ffloat(row.get(f"{metric}_snull_mean"))
            sn_sd = ffloat(row.get(f"{metric}_snull_std"))
            if math.isfinite(obs) and math.isfinite(sn_mu):
                rec[f"dyn__{cond}__{band}__{metric}_excess"] = obs - sn_mu
            if math.isfinite(obs) and math.isfinite(sn_mu) and math.isfinite(sn_sd) and sn_sd > 0:
                rec[f"dyn__{cond}__{band}__{metric}_z"] = (obs - sn_mu) / sn_sd


def feature_columns(rows: list[dict[str, object]]) -> list[str]:
    meta = {"sid", "LPS", "WST", "age", "cohort", "cohort_code"}
    cols = []
    for col in sorted({k for r in rows for k in r} - meta):
        vals = np.array([ffloat(r.get(col)) for r in rows], dtype=float)
        finite = vals[np.isfinite(vals)]
        if len(finite) >= 30 and float(finite.max() - finite.min()) > 1e-12:
            cols.append(col)
    return cols


def array(rows: list[dict[str, object]], col: str) -> np.ndarray:
    return np.array([ffloat(r.get(col)) for r in rows], dtype=float)


def build_rows() -> list[dict[str, object]]:
    subjects: dict[str, dict[str, object]] = {}
    add_static_features(subjects)
    add_trajectory_features(subjects)
    add_compressed_features(subjects)
    add_dynamic_features(subjects)
    rows = []
    for sid, rec in subjects.items():
        if not all(math.isfinite(ffloat(rec.get(k))) for k in ("LPS", "WST", "age")):
            continue
        age = ffloat(rec["age"])
        coh = cohort_name(age)
        if coh == "middle":
            continue
        rec = {"sid": sid, **rec, "cohort": coh, "cohort_code": 0.0 if coh == "young" else 1.0}
        rows.append(rec)
    rows.sort(key=lambda r: str(r["sid"]))
    return rows


def nuisance_columns(cols: list[str]) -> list[str]:
    keep = []
    for c in cols:
        name = c.lower()
        if "pow_" in name or "rel_power" in name or "meanplv" in name or "meanwpli" in name:
            keep.append(c)
        elif "mean_abs_corr" in name or "theta_delta_ratio" in name or "alpha_delta_ratio" in name:
            keep.append(c)
        elif "hi_freq_noise" in name:
            keep.append(c)
    return keep


def score_feature_set(rows: list[dict[str, object]], cols: list[str], weights: dict[str, float], name: str) -> dict[str, object]:
    X = np.column_stack([zscore(array(rows, c)) for c in weights])
    w = np.array([weights[c] for c in weights], dtype=float)
    w = w / (np.sqrt(np.sum(w**2)) + 1e-12)
    score = X @ w
    return {"name": name, "score": score, "components": ";".join(f"{k}:{v:g}" for k, v in weights.items())}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = build_rows()
    cols = feature_columns(rows)
    nuisance = nuisance_columns(cols)
    LPS = array(rows, "LPS")
    WST = array(rows, "WST")
    age = array(rows, "age")
    cohort = array(rows, "cohort_code")
    young_mask = np.array([r["cohort"] == "young" for r in rows], dtype=bool)

    base_controls = [cohort, age]
    wst_controls = [cohort, age, WST]
    nuisance_controls = base_controls + [array(rows, c) for c in nuisance]
    nuisance_wst_controls = nuisance_controls + [WST]

    result_rows: list[dict[str, object]] = []
    for i, col in enumerate(cols):
        x = array(rows, col)
        r_base, p_base, n_base = partial_spearman(x, LPS, base_controls)
        r_wst, p_wst, n_wst = partial_spearman(x, LPS, wst_controls)
        r_resid, p_resid, n_resid = partial_spearman(x, LPS, nuisance_controls)
        r_resid_wst, p_resid_wst, _ = partial_spearman(x, LPS, nuisance_wst_controls)
        r_young, p_young, n_young = partial_spearman(x[young_mask], LPS[young_mask], [age[young_mask]])
        r_wst_target, p_wst_target, _ = partial_spearman(x, WST, base_controls)
        result_rows.append(
            {
                "feature": col,
                "family": col.split("__", 1)[0],
                "n_base": n_base,
                "r_lps_cohort_age": r_base,
                "p_lps_cohort_age": p_base,
                "r_lps_plus_wst": r_wst,
                "p_lps_plus_wst": p_wst,
                "r_lps_power_residual": r_resid,
                "p_lps_power_residual": p_resid,
                "r_lps_power_residual_plus_wst": r_resid_wst,
                "p_lps_power_residual_plus_wst": p_resid_wst,
                "perm_p_power_residual": perm_p(x, LPS, nuisance_controls, abs(r_resid), 2000, 17 + i),
                "r_lps_young_only": r_young,
                "p_lps_young_only": p_young,
                "n_young": n_young,
                "r_wst_cohort_age": r_wst_target,
                "p_wst_cohort_age": p_wst_target,
                "abs_r_resid": abs(r_resid) if math.isfinite(r_resid) else math.nan,
            }
        )
    bh_q(result_rows, "p_lps_power_residual", "q_lps_power_residual")
    result_rows.sort(key=lambda r: (-float(r["abs_r_resid"]) if math.isfinite(float(r["abs_r_resid"])) else 1, float(r["p_lps_power_residual"]) if math.isfinite(float(r["p_lps_power_residual"])) else 1))
    fields = list(result_rows[0])
    write_csv(OUT / "single_feature_results.csv", result_rows, fields)

    composites = {}
    for prefix in ("dyn__EO__theta", "dyn__EC__theta", "dyn__EO__beta", "dyn__EO__gamma"):
        candidates = {
            f"{prefix}__switch_obs": 1.0,
            f"{prefix}__switch_z": 1.0,
            f"{prefix}__h_std_obs": 1.0,
            f"{prefix}__h_range_obs": 1.0,
            f"{prefix}__traj_dim_obs": -1.0,
            f"{prefix}__beam_depth_obs": 1.0,
        }
        composites[f"{prefix.replace('__', '_')}_reconfig"] = {k: v for k, v in candidates.items() if k in cols}
    composites["theta_trajectory_capacity"] = {
        k: v
        for k, v in {
            "cs__cs_traj_dim_theta": -1.0,
            "cs__elec_traj_dim_theta": -1.0,
            "cs__cs_steering_theta": -1.0,
            "traj__theta_steer": -1.0,
            "traj__traj_dim": -1.0,
        }.items()
        if k in cols
    }

    comp_rows = []
    for name, weights in composites.items():
        if len(weights) < 2:
            continue
        comp = score_feature_set(rows, cols, weights, name)
        x = comp["score"]
        r_base, p_base, n_base = partial_spearman(x, LPS, base_controls)
        r_wst, p_wst, _ = partial_spearman(x, LPS, wst_controls)
        r_resid, p_resid, _ = partial_spearman(x, LPS, nuisance_controls)
        r_resid_wst, p_resid_wst, _ = partial_spearman(x, LPS, nuisance_wst_controls)
        r_young, p_young, n_young = partial_spearman(x[young_mask], LPS[young_mask], [age[young_mask]])
        r_wst_target, p_wst_target, _ = partial_spearman(x, WST, base_controls)
        comp_rows.append(
            {
                "composite": name,
                "n_features": len(weights),
                "n_base": n_base,
                "r_lps_cohort_age": r_base,
                "p_lps_cohort_age": p_base,
                "r_lps_plus_wst": r_wst,
                "p_lps_plus_wst": p_wst,
                "r_lps_power_residual": r_resid,
                "p_lps_power_residual": p_resid,
                "r_lps_power_residual_plus_wst": r_resid_wst,
                "p_lps_power_residual_plus_wst": p_resid_wst,
                "perm_p_power_residual": perm_p(x, LPS, nuisance_controls, abs(r_resid), 5000, 900 + len(comp_rows)),
                "r_lps_young_only": r_young,
                "p_lps_young_only": p_young,
                "n_young": n_young,
                "r_wst_cohort_age": r_wst_target,
                "p_wst_cohort_age": p_wst_target,
                "components": comp["components"],
                "abs_r_resid": abs(r_resid) if math.isfinite(r_resid) else math.nan,
            }
        )
    bh_q(comp_rows, "p_lps_power_residual", "q_lps_power_residual")
    comp_rows.sort(key=lambda r: -float(r["abs_r_resid"]))
    write_csv(OUT / "composite_results.csv", comp_rows, list(comp_rows[0]) if comp_rows else [])

    subject_fields = ["sid", "LPS", "WST", "age", "cohort", "cohort_code"] + sorted(
        {k for r in rows for k in r} - {"sid", "LPS", "WST", "age", "cohort", "cohort_code"}
    )
    write_csv(OUT / "analysis_subjects.csv", rows, subject_fields)

    summary = [
        "# LEMON Intelligence Residual Search",
        "",
        f"Subjects: {len(rows)}",
        f"Cohorts: {dict(Counter(r['cohort'] for r in rows))}",
        f"Features searched: {len(cols)}",
        f"Nuisance features residualized in strict model: {len(nuisance)}",
        "",
        "The strict column is LPS partial Spearman after controlling cohort, age-bin midpoint, and power/coherence nuisance features. The WST column checks whether the feature is fluid-specific rather than vocabulary/general age structure.",
        "",
        "## Best Strict Single Features",
        "",
    ]
    for row in result_rows[:20]:
        summary.append(
            "- {feature}: r_strict={r:+.3f}, p={p:.4g}, perm_p={pp:.4g}, "
            "q={q}, strict+WST={rsw:+.3f}, young_r={ry:+.3f}, WST_r={wst:+.3f}".format(
                feature=row["feature"],
                r=float(row["r_lps_power_residual"]),
                p=float(row["p_lps_power_residual"]),
                pp=float(row["perm_p_power_residual"]),
                q=("" if row["q_lps_power_residual"] == "" else f"{float(row['q_lps_power_residual']):.4g}"),
                rsw=float(row["r_lps_power_residual_plus_wst"]),
                ry=float(row["r_lps_young_only"]) if math.isfinite(float(row["r_lps_young_only"])) else math.nan,
                wst=float(row["r_wst_cohort_age"]),
            )
        )
    summary.extend(["", "## Composites", ""])
    for row in comp_rows:
        summary.append(
            "- {name}: r_strict={r:+.3f}, p={p:.4g}, perm_p={pp:.4g}, "
            "strict+WST={rsw:+.3f}, young_r={ry:+.3f}, WST_r={wst:+.3f}".format(
                name=row["composite"],
                r=float(row["r_lps_power_residual"]),
                p=float(row["p_lps_power_residual"]),
                pp=float(row["perm_p_power_residual"]),
                rsw=float(row["r_lps_power_residual_plus_wst"]),
                ry=float(row["r_lps_young_only"]) if math.isfinite(float(row["r_lps_young_only"])) else math.nan,
                wst=float(row["r_wst_cohort_age"]),
            )
        )
    summary.extend(
        [
            "",
            "## Guardrail",
            "",
            "- This is still a same-dataset search over already generated features.",
            "- A serious LEMON hit should keep its direction with cohort+age control, remain nontrivial after WST control, survive the power/coherence residual model, and not simply become a WST correlate.",
            "- Young-only survival is the cleanest age-cohort stress test, but it has less power.",
        ]
    )
    (OUT / "SUMMARY.md").write_text("\n".join(summary) + "\n")
    print(f"Wrote {OUT / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
