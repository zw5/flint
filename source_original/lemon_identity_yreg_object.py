#!/usr/bin/env python3
"""Build a multivariate Yreg object from frozen identity coordinates.

This is not a per-label mining pass. It treats the label table as a multivariate
response Y and asks what low-rank part of Y is recoverable from the frozen
all-band/all-electrode identity coordinates X.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from lemon_shared_identity_label_scan import IDENTITY, IDENTITY_COLS, load_numeric_labels


ROOT = Path("/Users/ximon/Documents/interp")
OUT = ROOT / "experiments/eeg_coherence/lemon_identity_yreg_object"


def ffloat(value: object) -> float:
    try:
        out = float(str(value).replace(",", "."))
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def zscore_matrix(M: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = np.nanmean(M, axis=0)
    sd = np.nanstd(M, axis=0) + 1e-9
    return np.nan_to_num((M - mean) / sd), mean, sd


def load_rows() -> list[dict[str, object]]:
    identity_rows = {r["sid"]: dict(r) for r in read_csv(IDENTITY)}
    labels = load_numeric_labels()
    rows = []
    for sid, row in identity_rows.items():
        rec = dict(row)
        rec.update(labels.get(sid, {}))
        rows.append(rec)
    rows.sort(key=lambda r: str(r["sid"]))
    return rows


def label_columns(rows: list[dict[str, object]], include_availability: bool = False) -> list[str]:
    identity_keys = set(read_csv(IDENTITY)[0])
    cols = []
    for col in sorted({k for r in rows for k in r} - identity_keys):
        if not include_availability and "Data_Availability_Tables" in col:
            continue
        vals = np.asarray([ffloat(r.get(col)) for r in rows], dtype=float)
        finite = vals[np.isfinite(vals)]
        if len(finite) >= 60 and float(finite.max() - finite.min()) > 1e-12:
            cols.append(col)
    return cols


def matrix(rows: list[dict[str, object]], cols: list[str]) -> np.ndarray:
    return np.asarray([[ffloat(r.get(c)) for c in cols] for r in rows], dtype=float)


def corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if int(mask.sum()) < 5:
        return math.nan
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def component_reconstruction_r2(Y: np.ndarray, Yhat: np.ndarray) -> float:
    ss_res = float(np.nansum((Y - Yhat) ** 2))
    ss_tot = float(np.nansum((Y - np.nanmean(Y, axis=0, keepdims=True)) ** 2))
    return 1.0 - ss_res / max(ss_tot, 1e-12)


def main() -> None:
    rows = load_rows()
    x_cols = IDENTITY_COLS
    y_cols = label_columns(rows, include_availability=False)
    Xraw = matrix(rows, x_cols)
    Yraw = matrix(rows, y_cols)

    # Impute each label by its mean after z-scoring so incomplete labels can
    # still participate in a single multivariate response object.
    X, x_mean, x_sd = zscore_matrix(Xraw)
    Y, y_mean, y_sd = zscore_matrix(Yraw)

    C = (X.T @ Y) / max(len(rows) - 1, 1)
    U, S, Vt = np.linalg.svd(C, full_matrices=False)
    rank = min(len(S), len(x_cols))
    T = X @ U[:, :rank]
    Yreg = T @ np.diag(S[:rank]) @ Vt[:rank, :]

    subject_rows = []
    for i, row in enumerate(rows):
        rec: dict[str, object] = {"sid": row["sid"]}
        for j in range(rank):
            rec[f"yreg_axis{j + 1}"] = float(T[i, j])
        rec["yreg_norm"] = float(np.linalg.norm(Yreg[i]))
        rec["identity_norm"] = float(np.linalg.norm(X[i]))
        subject_rows.append(rec)

    axis_rows = []
    loading_rows = []
    for j in range(rank):
        axis_rows.append(
            {
                "axis": j + 1,
                "singular_value": float(S[j]),
                "identity_loading": ";".join(
                    f"{name}:{U[k, j]:+.4f}" for k, name in enumerate(x_cols)
                ),
                "subject_score_sd": float(np.nanstd(T[:, j])),
            }
        )
        load = Vt[j, :]
        order = np.argsort(-np.abs(load))
        for idx in order[:40]:
            y = Y[:, idx]
            loading_rows.append(
                {
                    "axis": j + 1,
                    "label": y_cols[idx],
                    "loading": float(load[idx]),
                    "corr_with_axis_score": corr(T[:, j], y),
                    "n": int(np.isfinite(Yraw[:, idx]).sum()),
                }
            )

    # Leave-one-subject-out reconstruction: learn the X->Y cross-covariance on
    # all other subjects and reconstruct the held subject's Y vector.
    Yhat_cv = np.full_like(Y, np.nan)
    for held in range(len(rows)):
        train = np.ones(len(rows), dtype=bool)
        train[held] = False
        Xtr, xm, xs = zscore_matrix(Xraw[train])
        Ytr, ym, ys = zscore_matrix(Yraw[train])
        Xh = np.nan_to_num((Xraw[held : held + 1] - xm) / xs)
        Ctr = (Xtr.T @ Ytr) / max(int(train.sum()) - 1, 1)
        Utr, Str, Vttr = np.linalg.svd(Ctr, full_matrices=False)
        rr = min(rank, len(Str))
        Yhat_cv[held] = Xh @ Utr[:, :rr] @ np.diag(Str[:rr]) @ Vttr[:rr, :]

    label_perf = []
    for idx, col in enumerate(y_cols):
        label_perf.append(
            {
                "label": col,
                "n": int(np.isfinite(Yraw[:, idx]).sum()),
                "cv_corr": corr(Y[:, idx], Yhat_cv[:, idx]),
                "in_sample_corr": corr(Y[:, idx], Yreg[:, idx]),
                "axis1_loading": float(Vt[0, idx]),
                "axis2_loading": float(Vt[1, idx]) if rank > 1 else math.nan,
                "axis3_loading": float(Vt[2, idx]) if rank > 2 else math.nan,
            }
        )
    label_perf.sort(key=lambda r: abs(ffloat(r["cv_corr"])), reverse=True)

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "yreg_subject_coordinates.csv", subject_rows)
    write_csv(OUT / "yreg_axes.csv", axis_rows)
    write_csv(OUT / "yreg_axis_label_loadings.csv", loading_rows)
    write_csv(OUT / "yreg_label_reconstruction.csv", label_perf)
    np.savez_compressed(
        OUT / "yreg_object.npz",
        sid=np.asarray([str(r["sid"]) for r in rows]),
        x_cols=np.asarray(x_cols),
        y_cols=np.asarray(y_cols),
        X=X,
        Y=Y,
        U=U,
        S=S,
        Vt=Vt,
        T=T,
        Yreg=Yreg,
        Yhat_cv=Yhat_cv,
        x_mean=x_mean,
        x_sd=x_sd,
        y_mean=y_mean,
        y_sd=y_sd,
    )

    lines = [
        "# LEMON Identity to Yreg Object",
        "",
        "Frozen X: all-band/all-electrode shared identity coordinates.",
        "Y: all substantive numeric LEMON labels jointly, excluding availability-table flags.",
        f"Subjects: {len(rows)}",
        f"Identity dimensions: {len(x_cols)}",
        f"Y labels: {len(y_cols)}",
        f"In-sample full-rank Y reconstruction R2: {component_reconstruction_r2(Y, Yreg):.4f}",
        f"Leave-one-subject-out full-rank Y reconstruction R2: {component_reconstruction_r2(Y, Yhat_cv):.4f}",
        "",
        "## Axes",
        "",
    ]
    for row in axis_rows:
        lines.append(
            f"- axis {row['axis']}: singular={row['singular_value']:.4f}, identity=({row['identity_loading']})"
        )
        top = [r for r in loading_rows if r["axis"] == row["axis"]][:12]
        for load in top:
            lines.append(
                "  - {label}: loading={loading:+.4f}, corr(axis,label)={corr_with_axis_score:+.3f}, n={n}".format(
                    **load
                )
            )
    lines.extend(["", "## Best Cross-Validated Label Reconstructions", ""])
    for row in label_perf[:30]:
        lines.append(
            "- {label}: cv_corr={cv_corr:+.3f}, in_sample_corr={in_sample_corr:+.3f}, n={n}".format(
                label=row["label"],
                cv_corr=ffloat(row["cv_corr"]),
                in_sample_corr=ffloat(row["in_sample_corr"]),
                n=row["n"],
            )
        )
    lines.extend(
        [
            "",
            "## Reading Rule",
            "",
            "Read axes as the phenotype object jointly recoverable from identity coordinates, not as individually selected label hits. Cross-validated reconstruction is the reality check; in-sample loadings describe the object but are not proof of label-level prediction.",
        ]
    )
    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
