#!/usr/bin/env python3
"""Domain-level Yreg recovery from frozen shared identity coordinates.

This tests whether broad physiological domains are more recoverable from the
shared identity coordinate than individual/personality domains.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from lemon_identity_yreg_highdim import (
    IDENTITY_NPZ,
    component_reconstruction_r2,
    corr,
    label_columns,
    load_rows,
    matrix,
    reliability_basis,
    write_csv,
    zscore_matrix,
)


ROOT = Path("/Users/ximon/Documents/interp")
OUT = ROOT / "experiments/eeg_coherence/lemon_identity_yreg_domain_recovery"


def domain_for(label: str) -> str:
    if "Anthropometry" in label:
        return "body_size"
    if "Blood_Pressure" in label:
        return "blood_pressure"
    if "Blood_Sample" in label:
        return "blood_chemistry"
    if "LPS__" in label or "WST__" in label or "TMT__" in label or "TAP_" in label or "CVLT" in label or "RWT" in label:
        return "cognition"
    if "Emotion_and_Personality" in label or "META_File" in label:
        return "personality_affect"
    return "other"


def project(avg_raw: np.ndarray, axes: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    return np.nan_to_num((avg_raw - mean) / sd) @ axes


def yreg_predict(Xraw: np.ndarray, Yraw: np.ndarray, rank: int) -> tuple[np.ndarray, np.ndarray]:
    X, _xm, _xs = zscore_matrix(Xraw)
    Y, _ym, _ys = zscore_matrix(Yraw)
    C = (X.T @ Y) / max(len(X) - 1, 1)
    U, S, Vt = np.linalg.svd(C, full_matrices=False)
    r = min(rank, len(S))
    Yhat = X @ U[:, :r] @ np.diag(S[:r]) @ Vt[:r, :]
    return Y, Yhat


def loocv_domain(
    odd: np.ndarray,
    even: np.ndarray,
    avg: np.ndarray,
    Yraw: np.ndarray,
    k: int,
    rank: int,
) -> tuple[np.ndarray, np.ndarray]:
    Yhat = np.full_like(Yraw, np.nan, dtype=float)
    for held in range(len(avg)):
        train = np.ones(len(avg), dtype=bool)
        train[held] = False
        axes, mean, sd, _s = reliability_basis(odd[train], even[train], k)
        Xtr = project(avg[train], axes, mean, sd)
        Xh = project(avg[held : held + 1], axes, mean, sd)
        Xtrz, xm, xs = zscore_matrix(Xtr)
        Ytrz, ym, ys = zscore_matrix(Yraw[train])
        Xhz = np.nan_to_num((Xh - xm) / xs)
        C = (Xtrz.T @ Ytrz) / max(int(train.sum()) - 1, 1)
        U, S, Vt = np.linalg.svd(C, full_matrices=False)
        r = min(rank, len(S))
        Yhat[held] = Xhz @ U[:, :r] @ np.diag(S[:r]) @ Vt[:r, :]
    Yz, _ym, _ys = zscore_matrix(Yraw)
    return Yz, Yhat


def mean_label_corr(Y: np.ndarray, Yhat: np.ndarray) -> float:
    vals = [corr(Y[:, i], Yhat[:, i]) for i in range(Y.shape[1])]
    vals = [v for v in vals if math.isfinite(v)]
    return float(np.nanmean(vals)) if vals else math.nan


def abs_mean_label_corr(Y: np.ndarray, Yhat: np.ndarray) -> float:
    vals = [abs(corr(Y[:, i], Yhat[:, i])) for i in range(Y.shape[1])]
    vals = [v for v in vals if math.isfinite(v)]
    return float(np.nanmean(vals)) if vals else math.nan


def main() -> None:
    rows = load_rows()
    y_cols = label_columns(rows, include_availability=False)
    Yall = matrix(rows, y_cols)
    z = np.load(IDENTITY_NPZ)
    sids = [str(x) for x in z["sid"]]
    order = [sids.index(str(r["sid"])) for r in rows if str(r["sid"]) in sids]
    rows = [r for r in rows if str(r["sid"]) in sids]
    Yall = matrix(rows, y_cols)
    odd = z["odd"][order]
    even = z["even"][order]
    avg = 0.5 * (odd + even)

    domains = sorted(set(domain_for(c) for c in y_cols))
    dims = [3, 6, 10, 20]
    rows_out = []
    label_rows = []
    for domain in domains:
        idx = [i for i, c in enumerate(y_cols) if domain_for(c) == domain]
        if len(idx) < 2:
            continue
        Yraw = Yall[:, idx]
        labels = [y_cols[i] for i in idx]
        for k in dims:
            axes, mean, sd, _s = reliability_basis(odd, even, k)
            Xk = project(avg, axes, mean, sd)
            rank = min(k, 10)
            Yz, Yhat = yreg_predict(Xk, Yraw, rank)
            Yz_cv, Yhat_cv = loocv_domain(odd, even, avg, Yraw, k, rank)
            rows_out.append(
                {
                    "domain": domain,
                    "n_labels": len(idx),
                    "k_identity_dims": k,
                    "rank": rank,
                    "in_sample_r2": component_reconstruction_r2(Yz, Yhat),
                    "loocv_r2": component_reconstruction_r2(Yz_cv, Yhat_cv),
                    "mean_cv_corr": mean_label_corr(Yz_cv, Yhat_cv),
                    "mean_abs_cv_corr": abs_mean_label_corr(Yz_cv, Yhat_cv),
                }
            )
            for j, label in enumerate(labels):
                label_rows.append(
                    {
                        "domain": domain,
                        "k_identity_dims": k,
                        "label": label,
                        "n": int(np.isfinite(Yraw[:, j]).sum()),
                        "cv_corr": corr(Yz_cv[:, j], Yhat_cv[:, j]),
                        "in_sample_corr": corr(Yz[:, j], Yhat[:, j]),
                    }
                )
    label_rows.sort(key=lambda r: (r["domain"], int(r["k_identity_dims"]), -abs(float(r["cv_corr"]) if math.isfinite(float(r["cv_corr"])) else 0.0)))

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "domain_recovery_summary.csv", rows_out)
    write_csv(OUT / "domain_label_reconstruction.csv", label_rows)

    lines = [
        "# LEMON Identity Yreg Domain Recovery",
        "",
        "Question: are shared physiology/metabolic domains more recoverable from broad identity coordinates than personality/affect labels?",
        f"Subjects: {len(rows)}",
        "",
        "## Domain Summary",
        "",
    ]
    for row in rows_out:
        lines.append(
            "- {domain}, K={k_identity_dims}, labels={n_labels}: LOOCV_R2={loocv_r2:+.4f}, "
            "mean_cv_corr={mean_cv_corr:+.3f}, mean_abs_cv_corr={mean_abs_cv_corr:+.3f}, in_R2={in_sample_r2:+.4f}".format(
                **row
            )
        )
    lines.extend(["", "## Best Labels Per Domain At K=6", ""])
    for domain in domains:
        top = [r for r in label_rows if r["domain"] == domain and r["k_identity_dims"] == 6][:12]
        if not top:
            continue
        lines.append(f"### {domain}")
        for row in top:
            lines.append(
                "- {label}: cv_corr={cv_corr:+.3f}, in_sample_corr={in_sample_corr:+.3f}, n={n}".format(
                    label=row["label"],
                    cv_corr=float(row["cv_corr"]),
                    in_sample_corr=float(row["in_sample_corr"]),
                    n=row["n"],
                )
            )
    lines.extend(
        [
            "",
            "## Reading Rule",
            "",
            "Use domain-level LOOCV and mean label correlations to compare shared recoverability. Individual label rows explain which labels drive each domain, but the domain summary is the main test.",
        ]
    )
    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
