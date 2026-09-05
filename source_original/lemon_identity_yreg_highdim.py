#!/usr/bin/env python3
"""Higher-dimensional identity-coordinate Yreg sweep.

This expands the frozen identity coordinate beyond the compact 5-axis summary.
The raw identity tensor has 285 dimensions per split. We derive a reliability
basis from odd/even identity agreement, keep the top K coordinates, and ask how
much of the multivariate phenotype field Y is recoverable as K grows.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from lemon_identity_yreg_object import (
    OUT as LOWDIM_OUT,
    component_reconstruction_r2,
    corr,
    ffloat,
    label_columns,
    load_rows,
    matrix,
    write_csv,
    zscore_matrix,
)


ROOT = Path("/Users/ximon/Documents/interp")
IDENTITY_NPZ = ROOT / "experiments/eeg_coherence/lemon_shared_identity_sex_beam/shared_identity_full.npz"
OUT = ROOT / "experiments/eeg_coherence/lemon_identity_yreg_highdim"


def reliability_basis(odd_raw: np.ndarray, even_raw: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    odd, mean, sd = zscore_matrix(odd_raw)
    even, _, _ = zscore_matrix(even_raw)
    C = (odd.T @ even) / max(len(odd) - 1, 1)
    u, s, vh = np.linalg.svd(C, full_matrices=False)
    axes = []
    for i in range(min(k, len(s))):
        axis = u[:, i] + vh[i, :]
        axis = axis / (np.linalg.norm(axis) + 1e-12)
        axes.append(axis)
    return np.stack(axes, axis=1), mean, sd, s


def project(avg_raw: np.ndarray, axes: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    X = np.nan_to_num((avg_raw - mean) / sd)
    return X @ axes


def yreg_fit_predict(Xraw: np.ndarray, Yraw: np.ndarray, rank: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    X, _xm, _xs = zscore_matrix(Xraw)
    Y, _ym, _ys = zscore_matrix(Yraw)
    C = (X.T @ Y) / max(len(X) - 1, 1)
    U, S, Vt = np.linalg.svd(C, full_matrices=False)
    r = min(rank or len(S), len(S))
    T = X @ U[:, :r]
    Yhat = T @ np.diag(S[:r]) @ Vt[:r, :]
    return Y, Yhat, U, S, Vt


def loocv_predict(
    odd_raw: np.ndarray,
    even_raw: np.ndarray,
    avg_raw: np.ndarray,
    Yraw: np.ndarray,
    k: int,
    rank: int | None = None,
) -> np.ndarray:
    Yhat = np.full_like(Yraw, np.nan, dtype=float)
    for held in range(len(avg_raw)):
        train = np.ones(len(avg_raw), dtype=bool)
        train[held] = False
        axes, mean, sd, _s = reliability_basis(odd_raw[train], even_raw[train], k)
        Xtr = project(avg_raw[train], axes, mean, sd)
        Xh = project(avg_raw[held : held + 1], axes, mean, sd)
        Xtrz, xm, xs = zscore_matrix(Xtr)
        Ytrz, ym, ys = zscore_matrix(Yraw[train])
        Xhz = np.nan_to_num((Xh - xm) / xs)
        C = (Xtrz.T @ Ytrz) / max(int(train.sum()) - 1, 1)
        U, S, Vt = np.linalg.svd(C, full_matrices=False)
        r = min(rank or len(S), len(S))
        Yhat[held] = Xhz @ U[:, :r] @ np.diag(S[:r]) @ Vt[:r, :]
    return Yhat


def main() -> None:
    rows = load_rows()
    y_cols = label_columns(rows, include_availability=False)
    Yraw = matrix(rows, y_cols)
    z = np.load(IDENTITY_NPZ)
    sids = [str(x) for x in z["sid"]]
    row_by_sid = {str(r["sid"]): i for i, r in enumerate(rows)}
    order = [sids.index(str(r["sid"])) for r in rows if str(r["sid"]) in sids]
    rows = [r for r in rows if str(r["sid"]) in sids]
    Yraw = matrix(rows, y_cols)
    odd = z["odd"][order]
    even = z["even"][order]
    first = z["first"][order]
    second = z["second"][order]
    avg = 0.5 * (odd + even)

    dims = [3, 6, 10, 20, 40, 80]
    summary_rows = []
    perf_rows = []
    axis_rows = []
    for k in dims:
        axes, mean, sd, reliability_s = reliability_basis(odd, even, k)
        Xk = project(avg, axes, mean, sd)
        Yz, Yhat, U, S, Vt = yreg_fit_predict(Xk, Yraw)
        # Use low rank capped at min(k, 20) for CV to avoid making K=80 only a
        # high-variance memorizer. This is still intentionally low-assumption.
        cv_rank = min(k, 20)
        Yhat_cv = loocv_predict(odd, even, avg, Yraw, k, rank=cv_rank)
        row = {
            "k_identity_dims": k,
            "cv_rank": cv_rank,
            "in_sample_r2": component_reconstruction_r2(Yz, Yhat),
            "loocv_r2": component_reconstruction_r2(Yz, Yhat_cv),
            "first_reliability_singular": float(reliability_s[0]),
            "kth_reliability_singular": float(reliability_s[min(k - 1, len(reliability_s) - 1)]),
            "first_yreg_singular": float(S[0]) if len(S) else math.nan,
        }
        summary_rows.append(row)
        for axis in range(min(8, len(S))):
            load = Vt[axis]
            order_labels = np.argsort(-np.abs(load))
            for idx in order_labels[:15]:
                axis_rows.append(
                    {
                        "k_identity_dims": k,
                        "axis": axis + 1,
                        "singular_value": float(S[axis]),
                        "label": y_cols[idx],
                        "loading": float(load[idx]),
                        "corr_with_axis_score": corr((zscore_matrix(Xk)[0] @ U[:, axis]), Yz[:, idx]),
                        "n": int(np.isfinite(Yraw[:, idx]).sum()),
                    }
                )
        for idx, col in enumerate(y_cols):
            perf_rows.append(
                {
                    "k_identity_dims": k,
                    "label": col,
                    "n": int(np.isfinite(Yraw[:, idx]).sum()),
                    "cv_corr": corr(Yz[:, idx], Yhat_cv[:, idx]),
                    "in_sample_corr": corr(Yz[:, idx], Yhat[:, idx]),
                }
            )

    perf_rows.sort(key=lambda r: (int(r["k_identity_dims"]), -abs(ffloat(r["cv_corr"]))))
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "highdim_yreg_summary.csv", summary_rows)
    write_csv(OUT / "highdim_yreg_axis_loadings.csv", axis_rows)
    write_csv(OUT / "highdim_yreg_label_reconstruction.csv", perf_rows)

    lines = [
        "# LEMON Higher-Dimensional Identity to Yreg Sweep",
        "",
        "X source: raw all-band/all-electrode split identity tensor, projected into top K odd/even reliability coordinates.",
        "Y: substantive numeric LEMON labels jointly, excluding availability-table flags.",
        f"Subjects: {len(rows)}",
        f"Raw identity dimensions: {odd.shape[1]}",
        f"Y labels: {len(y_cols)}",
        "",
        "## Dimensional Sweep",
        "",
    ]
    for row in summary_rows:
        lines.append(
            "- K={k_identity_dims}: in_R2={in_sample_r2:+.4f}, LOOCV_R2={loocv_r2:+.4f}, "
            "cv_rank={cv_rank}, first_yreg_singular={first_yreg_singular:.4f}, kth_reliability_singular={kth_reliability_singular:.4f}".format(
                **row
            )
        )
    lines.extend(["", "## Axis 1 Top Loadings By K", ""])
    for k in dims:
        lines.append(f"### K={k}")
        top = [r for r in axis_rows if r["k_identity_dims"] == k and r["axis"] == 1][:12]
        for r in top:
            lines.append(
                "- {label}: loading={loading:+.4f}, corr={corr_with_axis_score:+.3f}, n={n}".format(**r)
            )
    lines.extend(["", "## Best Cross-Validated Reconstructions By K", ""])
    for k in dims:
        lines.append(f"### K={k}")
        top = [r for r in perf_rows if r["k_identity_dims"] == k][:15]
        for r in top:
            lines.append(
                "- {label}: cv_corr={cv_corr:+.3f}, in_sample_corr={in_sample_corr:+.3f}, n={n}".format(
                    label=r["label"],
                    cv_corr=ffloat(r["cv_corr"]),
                    in_sample_corr=ffloat(r["in_sample_corr"]),
                    n=r["n"],
                )
            )
    lines.extend(
        [
            "",
            "## Reading Rule",
            "",
            "If higher K raises in-sample R2 but worsens LOOCV R2, it is just unfolding noise. If the same phenotype bundles recur across K and remain cross-validated, they are stronger evidence for a real high-dimensional subject coordinate.",
        ]
    )
    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
