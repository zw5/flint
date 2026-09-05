#!/usr/bin/env python3
"""Label-blind shared identity beam test for LEMON sex labels.

The object here is deliberately not a best-column search. For every subject we
use all electrodes and a fixed set of canonical EEG bands. Each band contributes
the leading spatial envelope of a lagged cross-field operator in two time
splits. A reliability axis is learned from odd/even agreement without sex
labels, then that locked axis is tested against the LEMON sex label.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
from scipy.signal import butter, filtfilt, hilbert

import lemon_intelligence_residual_search as stats
from lemon_frequency_sliced_beam_eval import concentration_from_spectrum, static_plv, window_starts


warnings.filterwarnings("ignore")

ROOT = Path("/Users/ximon/Documents/interp")
EXP = ROOT / "experiments/eeg_coherence"
EXPORT = EXP / "rust_export"
META = ROOT / "data/lemon/Behavioural_Data_MPILMBB_LEMON/META_File_IDs_Age_Gender_Education_Drug_Smoke_SKID_LEMON.csv"
OUT = EXP / "lemon_shared_identity_sex_beam"

BANDS = {
    "delta1_4": (1.0, 4.0),
    "theta4_8": (4.0, 8.0),
    "alpha8_13": (8.0, 13.0),
    "beta13_30": (13.0, 30.0),
    "gamma30_45": (30.0, 45.0),
}
PHASE_FRACTION = 0.25


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
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


def cohort(age: float) -> str:
    if age <= 35:
        return "young"
    if age >= 55:
        return "old"
    return "middle"


def load_sex_meta() -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    with META.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        sex_idx = header.index("Gender_ 1=female_2=male")
        for raw in reader:
            if not raw or raw[sex_idx] not in {"1", "2"}:
                continue
            out[raw[0]] = {
                "sex_female": 1.0 if raw[sex_idx] == "1" else 0.0,
                "sex_label": "female" if raw[sex_idx] == "1" else "male",
            }
    return out


def load_manifest(limit: int = 0) -> list[dict[str, object]]:
    sex_meta = load_sex_meta()
    rows: list[dict[str, object]] = []
    for row in read_csv(EXPORT / "manifest.csv"):
        sid = str(row["sid"])
        if sid not in sex_meta:
            continue
        age = ffloat(row.get("age"))
        c = cohort(age)
        if c == "middle":
            continue
        rec: dict[str, object] = dict(row)
        rec["age"] = age
        rec["LPS"] = ffloat(row.get("LPS"))
        rec["WST"] = ffloat(row.get("WST"))
        rec["sfreq"] = ffloat(row.get("sfreq"))
        rec["n_ch"] = int(float(row["n_ch"]))
        rec["n_eo"] = int(float(row["n_eo"]))
        rec["cohort"] = c
        rec["cohort_code"] = 0.0 if c == "young" else 1.0
        rec.update(sex_meta[sid])
        if rec["n_eo"] > 0 and all(math.isfinite(float(rec[k])) for k in ("age", "sfreq")):
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
    if hi <= lo:
        raise ValueError(f"Band {lo}-{hi} is invalid for sfreq={sfreq}")
    b, a = butter(4, [lo / (sfreq / 2.0), hi / (sfreq / 2.0)], btype="band")
    return filtfilt(b, a, X, axis=0)


def lag_for_band(sfreq: float, lo: float, hi: float) -> int:
    center = 0.5 * (lo + hi)
    return max(1, int(round((PHASE_FRACTION / center) * sfreq)))


def cross_operator(Z: np.ndarray, lag: int) -> np.ndarray:
    present = Z[:-lag]
    future = Z[lag:]
    return (future.T @ present.conj()) / max(len(present), 1)


def split_identity(phase: np.ndarray, starts: list[int], win: int, lag: int) -> tuple[np.ndarray, dict[str, float]]:
    ops: list[np.ndarray] = []
    concs: list[float] = []
    effs: list[float] = []
    locks: list[float] = []
    for start in starts:
        Z = phase[start : start + win]
        if len(Z) <= lag + 10:
            continue
        C = cross_operator(Z, lag)
        try:
            u, s, _vh = np.linalg.svd(C, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        concentration, effrank, lock = concentration_from_spectrum(s)
        if not math.isfinite(concentration):
            continue
        # Use the spatial envelope so arbitrary complex phase/sign choices do
        # not masquerade as identity changes.
        env = np.abs(u[:, 0])
        env = env / (np.linalg.norm(env) + 1e-12)
        ops.append(env)
        concs.append(concentration)
        effs.append(effrank)
        locks.append(lock)
    if not ops:
        return np.array([], dtype=float), {"concentration": math.nan, "effrank": math.nan, "lock": math.nan}
    ident = np.nanmean(np.stack(ops, axis=0), axis=0)
    ident = ident / (np.linalg.norm(ident) + 1e-12)
    return ident, {
        "concentration": float(np.nanmean(concs)),
        "effrank": float(np.nanmean(effs)),
        "lock": float(np.nanmean(locks)),
    }


def split_starts(starts: list[int]) -> dict[str, list[int]]:
    half = len(starts) // 2
    return {
        "odd": starts[::2],
        "even": starts[1::2],
        "first": starts[:half],
        "second": starts[half:],
    }


def process_subject(job: tuple[dict[str, object], dict[str, object]]) -> dict[str, object] | None:
    row, cfg = job
    X = load_eo(row)
    if X is None:
        return None
    sfreq = float(row["sfreq"])
    win = int(round(float(cfg["window_sec"]) * sfreq))
    step = int(round(float(cfg["step_sec"]) * sfreq))
    starts = window_starts(len(X), win, step, int(cfg["max_windows"]))
    start_sets = split_starts(starts)
    if min(len(v) for v in start_sets.values()) < 2:
        return None

    rec: dict[str, object] = {
        "sid": row["sid"],
        "age": row["age"],
        "cohort": row["cohort"],
        "cohort_code": row["cohort_code"],
        "sex_female": row["sex_female"],
        "sex_label": row["sex_label"],
        "LPS": row.get("LPS", math.nan),
        "WST": row.get("WST", math.nan),
    }
    vectors = {"odd": [], "even": [], "first": [], "second": []}
    for band_name, (lo, hi) in BANDS.items():
        try:
            Xb = bandpass(X, sfreq, lo, hi)
            analytic = hilbert(Xb, axis=0)
            phase = analytic / (np.abs(analytic) + 1e-9)
        except Exception as exc:
            print(f"  {row['sid']} {band_name}: {exc}", flush=True)
            return None
        lag = lag_for_band(sfreq, lo, hi)
        rec[f"{band_name}_power"] = float(np.nanmean(Xb**2))
        rec[f"{band_name}_static_plv"] = static_plv(phase)
        for split, split_start_list in start_sets.items():
            ident, metrics = split_identity(phase, split_start_list, win, lag)
            if ident.size == 0:
                return None
            vectors[split].append(ident)
            for key, value in metrics.items():
                rec[f"{band_name}_{split}_{key}"] = value

    for split, pieces in vectors.items():
        vec = np.concatenate(pieces)
        vec = vec / (np.linalg.norm(vec) + 1e-12)
        rec[f"{split}_identity"] = vec
    rec["odd_even_identity_cos"] = float(np.dot(rec["odd_identity"], rec["even_identity"]))
    rec["first_second_identity_cos"] = float(np.dot(rec["first_identity"], rec["second_identity"]))
    rec["mean_identity_persistence"] = float(
        np.nanmean([rec["odd_even_identity_cos"], rec["first_second_identity_cos"]])
    )
    rec["mean_band_power"] = float(np.nanmean([rec[f"{name}_power"] for name in BANDS]))
    rec["mean_static_plv"] = float(np.nanmean([rec[f"{name}_static_plv"] for name in BANDS]))
    return rec


def matrix(rows: list[dict[str, object]], key: str) -> np.ndarray:
    return np.stack([np.asarray(r[key], dtype=float) for r in rows], axis=0)


def zmat(X: np.ndarray, mean: np.ndarray | None = None, sd: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if mean is None:
        mean = np.nanmean(X, axis=0)
    if sd is None:
        sd = np.nanstd(X, axis=0) + 1e-9
    return np.nan_to_num((X - mean) / sd), mean, sd


def reliable_axes(
    odd: np.ndarray,
    even: np.ndarray,
    n_components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    Xo, mean, sd = zmat(odd)
    Xe, _, _ = zmat(even, mean, sd)
    C = (Xo.T @ Xe) / max(len(Xo) - 1, 1)
    u, s, vh = np.linalg.svd(C, full_matrices=False)
    axes = []
    for i in range(min(n_components, len(s))):
        axis = u[:, i] + vh[i, :]
        axis = axis / (np.linalg.norm(axis) + 1e-12)
        axes.append(axis)
    return np.stack(axes, axis=0), mean, sd, s, C


def vec(rows: list[dict[str, object]], col: str) -> np.ndarray:
    return np.asarray([ffloat(r.get(col)) for r in rows], dtype=float)


def eval_score(rows: list[dict[str, object]], score_col: str, controls: list[str]) -> dict[str, object]:
    r, p, n = stats.partial_spearman(
        vec(rows, score_col),
        vec(rows, "sex_female"),
        [vec(rows, c) for c in controls],
    )
    return {"score": score_col, "controls": ";".join(controls), "n": n, "r_female": r, "p_female": p}


def auc_score(y: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    pos = score[y == 1.0]
    neg = score[y == 0.0]
    if len(pos) == 0 or len(neg) == 0:
        return math.nan
    wins = 0.0
    total = 0.0
    for p in pos:
        wins += float(np.sum(p > neg)) + 0.5 * float(np.sum(p == neg))
        total += len(neg)
    return wins / max(total, 1.0)


def residualize_matrix(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    if C.size == 0:
        return X - np.nanmean(X, axis=0, keepdims=True)
    B = np.column_stack([np.ones(len(X)), C])
    coef, *_ = np.linalg.lstsq(B, X, rcond=None)
    return X - B @ coef


def ridge_loocv(rows: list[dict[str, object]], feature_cols: list[str], control_cols: list[str], alpha: float = 1.0) -> dict[str, object]:
    X = np.column_stack([vec(rows, c) for c in feature_cols])
    y = vec(rows, "sex_female")
    C = np.column_stack([vec(rows, c) for c in control_cols]) if control_cols else np.empty((len(rows), 0))
    pred = np.full(len(rows), math.nan)
    for held in range(len(rows)):
        train = np.ones(len(rows), dtype=bool)
        train[held] = False
        mask_train = train & np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        if C.size:
            mask_train &= np.all(np.isfinite(C), axis=1)
        if int(mask_train.sum()) <= len(feature_cols) + len(control_cols) + 3:
            continue
        Xtr = X[mask_train]
        ytr = y[mask_train]
        Ctr = C[mask_train] if C.size else np.empty((len(Xtr), 0))
        mu = np.nanmean(Xtr, axis=0)
        sd = np.nanstd(Xtr, axis=0) + 1e-9
        Xtrz = (Xtr - mu) / sd
        Xtstz = (X[held : held + 1] - mu) / sd
        if C.size:
            cmu = np.nanmean(Ctr, axis=0)
            csd = np.nanstd(Ctr, axis=0) + 1e-9
            Ctrz = (Ctr - cmu) / csd
            Ctstz = (C[held : held + 1] - cmu) / csd
        else:
            Ctrz = np.empty((len(Xtrz), 0))
            Ctstz = np.empty((1, 0))
        Xtr_res = residualize_matrix(Xtrz, Ctrz)
        ytr_res = residualize_matrix(ytr[:, None], Ctrz).ravel()
        if C.size:
            Bc = np.column_stack([np.ones(len(Ctrz)), Ctrz])
            coef_c, *_ = np.linalg.lstsq(Bc, Xtrz, rcond=None)
            Xtst_res = Xtstz - np.column_stack([np.ones(1), Ctstz]) @ coef_c
        else:
            Xtst_res = Xtstz - np.nanmean(Xtrz, axis=0, keepdims=True)
        A = Xtr_res.T @ Xtr_res + alpha * np.eye(Xtr_res.shape[1])
        b = Xtr_res.T @ ytr_res
        beta = np.linalg.solve(A, b)
        pred[held] = float((Xtst_res @ beta).ravel()[0])
    mask = np.isfinite(pred) & np.isfinite(y)
    r, p, n = stats.partial_spearman(pred[mask], y[mask], [])
    acc = float(np.mean((pred[mask] >= np.nanmedian(pred[mask])) == (y[mask] == 1.0))) if int(mask.sum()) else math.nan
    return {
        "score": "loocv_ridge_" + "_".join(feature_cols),
        "controls": ";".join(control_cols),
        "n": n,
        "r_female": r,
        "p_female": p,
        "auc": auc_score(y[mask], pred[mask]),
        "median_split_accuracy": acc,
    }


def write_summary(path: Path, rows: list[dict[str, object]], eval_rows: list[dict[str, object]], singular: np.ndarray) -> None:
    counts = {label: sum(str(r["sex_label"]) == label for r in rows) for label in ("female", "male")}
    lines = [
        "# LEMON Shared Identity Sex Beam",
        "",
        "Goal: derive an all-band/all-electrode identity that persists across time without using sex labels, then test whether that identity tracks sex.",
        "",
        f"Subjects: {len(rows)}",
        f"Sex-label counts: {counts}",
        f"Bands: {', '.join(BANDS)}",
        f"Identity dimensions: {len(BANDS)} bands x all electrodes",
        f"Reliability singular values first five: {', '.join(f'{x:.4f}' for x in singular[:5])}",
        "",
        "## Sex Association",
        "",
    ]
    for row in eval_rows:
        lines.append(
            "- {score} | controls={controls}: n={n}, r={r_female:+.3f}, p={p_female:.4g}{extra}".format(
                score=row["score"],
                controls=row["controls"] or "none",
                n=row["n"],
                r_female=ffloat(row["r_female"]),
                p_female=ffloat(row["p_female"]),
                extra=(
                    f", AUC={ffloat(row.get('auc')):.3f}, median-threshold acc={ffloat(row.get('median_split_accuracy')):.3f}"
                    if "auc" in row
                    else ""
                ),
            )
        )
    lines.extend(
        [
            "",
            "## Reading Rule",
            "",
            "The reliable identity axis is learned only from odd/even identity agreement. Sex is used afterward as a test label. "
            "The strict line controls age, cohort, mean band power, and mean static PLV, so it asks whether the shared identity carries sex beyond broad spectral nuisance axes.",
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
    parser.add_argument("--components", type=int, default=5)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    suffix = "full" if not args.limit else f"limit{args.limit}"
    OUT.mkdir(parents=True, exist_ok=True)
    npz_path = OUT / f"shared_identity_{suffix}.npz"
    subject_path = OUT / f"shared_identity_subjects_{suffix}.csv"
    eval_path = OUT / f"shared_identity_eval_{suffix}.csv"
    summary_path = OUT / f"SUMMARY_{suffix}.md"

    manifest = load_manifest(args.limit)
    if npz_path.exists() and subject_path.exists() and not args.force:
        rows = read_csv(subject_path)
        data = np.load(npz_path)
        odd = data["odd"]
        even = data["even"]
        first = data["first"]
        second = data["second"]
    else:
        cfg = {"window_sec": args.window_sec, "step_sec": args.step_sec, "max_windows": args.max_windows}
        jobs = [(row, cfg) for row in manifest]
        rows_obj: list[dict[str, object]] = []
        print(f"Subjects queued: {len(jobs)}; workers={args.workers}", flush=True)
        if args.workers <= 1:
            for i, job in enumerate(jobs, start=1):
                print(f"  processing {i}/{len(jobs)} {job[0]['sid']}", flush=True)
                rec = process_subject(job)
                if rec is not None:
                    rows_obj.append(rec)
        else:
            with Pool(args.workers) as pool:
                for i, rec in enumerate(pool.imap_unordered(process_subject, jobs), start=1):
                    if rec is not None:
                        rows_obj.append(rec)
                    if i == 1 or i % 10 == 0:
                        print(f"  completed {i}/{len(jobs)} ({len(rows_obj)} kept)", flush=True)
        rows_obj.sort(key=lambda r: str(r["sid"]))
        odd = matrix(rows_obj, "odd_identity")
        even = matrix(rows_obj, "even_identity")
        first = matrix(rows_obj, "first_identity")
        second = matrix(rows_obj, "second_identity")
        np.savez_compressed(
            npz_path,
            sid=np.asarray([str(r["sid"]) for r in rows_obj]),
            odd=odd,
            even=even,
            first=first,
            second=second,
        )
        csv_rows = []
        for r in rows_obj:
            csv_rows.append({k: v for k, v in r.items() if not k.endswith("_identity")})
        write_csv(subject_path, csv_rows)
        rows = read_csv(subject_path)

    axes, mean, sd, singular, _C = reliable_axes(odd, even, args.components)
    avg_identity = 0.5 * (odd + even)
    Xavg, _, _ = zmat(avg_identity, mean, sd)
    component_cols = []
    for i, axis in enumerate(axes, start=1):
        col = f"reliable_identity_c{i}"
        component_cols.append(col)
        scores = Xavg @ axis
        for row, score in zip(rows, scores):
            row[col] = float(score)
    for row in rows:
        row["reliable_identity_score"] = ffloat(row.get("reliable_identity_c1"))
    # A second label-blind coordinate: temporal persistence itself.
    for row in rows:
        row["mean_identity_persistence"] = ffloat(row.get("mean_identity_persistence"))

    eval_rows = [
        eval_score(rows, "reliable_identity_score", []),
        eval_score(rows, "reliable_identity_score", ["age", "cohort_code"]),
        eval_score(rows, "reliable_identity_score", ["age", "cohort_code", "mean_band_power", "mean_static_plv"]),
        eval_score(rows, "mean_identity_persistence", []),
        eval_score(rows, "mean_identity_persistence", ["age", "cohort_code"]),
        eval_score(rows, "mean_identity_persistence", ["age", "cohort_code", "mean_band_power", "mean_static_plv"]),
    ]
    for col in component_cols[1:]:
        eval_rows.append(eval_score(rows, col, ["age", "cohort_code", "mean_band_power", "mean_static_plv"]))
    eval_rows.extend(
        [
            ridge_loocv(rows, component_cols, []),
            ridge_loocv(rows, component_cols, ["age", "cohort_code"]),
            ridge_loocv(rows, component_cols, ["age", "cohort_code", "mean_band_power", "mean_static_plv"]),
        ]
    )
    write_csv(eval_path, eval_rows)
    scored_path = OUT / f"shared_identity_scored_subjects_{suffix}.csv"
    write_csv(scored_path, rows)
    write_summary(summary_path, rows, eval_rows, singular)
    print(f"Wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
