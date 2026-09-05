#!/usr/bin/env python3
"""Scan LEMON labels against the frozen shared identity beam coordinates.

This script does not construct or select a new beam object. It reuses the
label-blind all-band/all-electrode identity coordinates from
lemon_shared_identity_sex_beam.py and tests available numeric LEMON labels
against those frozen coordinates.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

import lemon_intelligence_residual_search as stats


ROOT = Path("/Users/ximon/Documents/interp")
LEMON = ROOT / "data/lemon/Behavioural_Data_MPILMBB_LEMON"
IDENTITY = ROOT / "experiments/eeg_coherence/lemon_shared_identity_sex_beam/shared_identity_scored_subjects_full.csv"
OUT = ROOT / "experiments/eeg_coherence/lemon_shared_identity_label_scan"

IDENTITY_COLS = [
    "reliable_identity_c1",
    "reliable_identity_c2",
    "reliable_identity_c3",
    "reliable_identity_c4",
    "reliable_identity_c5",
    "mean_identity_persistence",
]


def ffloat(value: object) -> float:
    try:
        out = float(str(value).replace(",", "."))
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def read_csv_raw(path: Path) -> tuple[list[str], list[list[str]]]:
    with path.open(newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        return header, list(reader)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_identity() -> dict[str, dict[str, object]]:
    rows = read_csv(IDENTITY)
    return {str(r["sid"]): dict(r) for r in rows}


def load_numeric_labels() -> dict[str, dict[str, float]]:
    labels: dict[str, dict[str, float]] = {}
    for path in sorted(LEMON.rglob("*.csv")):
        header, rows = read_csv_raw(path)
        if not header:
            continue
        source = path.relative_to(LEMON).as_posix().replace("/", "__").replace(" ", "_")
        source = source[:-4] if source.endswith(".csv") else source
        for raw in rows:
            if not raw:
                continue
            sid = raw[0].strip()
            if not sid.startswith("sub-"):
                continue
            rec = labels.setdefault(sid, {})
            for idx, name in enumerate(header[1:], start=1):
                if idx >= len(raw):
                    continue
                value = ffloat(raw[idx])
                if not math.isfinite(value):
                    continue
                col = name.strip() or f"col{idx}"
                rec[f"{source}__{col}"] = value
    return labels


def vec(rows: list[dict[str, object]], col: str) -> np.ndarray:
    return np.asarray([ffloat(r.get(col)) for r in rows], dtype=float)


def eval_pair(rows: list[dict[str, object]], label: str, ident_col: str, controls: list[str]) -> dict[str, object]:
    r, p, n = stats.partial_spearman(
        vec(rows, ident_col),
        vec(rows, label),
        [vec(rows, c) for c in controls],
    )
    return {
        "label": label,
        "identity_coord": ident_col,
        "controls": ";".join(controls),
        "n": n,
        "r": r,
        "p": p,
    }


def main() -> None:
    identity = load_identity()
    identity_keys = set(next(iter(identity.values())).keys()) if identity else set()
    labels = load_numeric_labels()
    rows: list[dict[str, object]] = []
    for sid, rec in identity.items():
        merged = dict(rec)
        merged.update(labels.get(sid, {}))
        rows.append(merged)
    rows.sort(key=lambda r: str(r["sid"]))

    label_cols = sorted(
        c
        for c in {k for r in rows for k in r}
        if c not in identity_keys
    )
    usable = []
    for col in label_cols:
        vals = vec(rows, col)
        finite = vals[np.isfinite(vals)]
        if len(finite) >= 60 and float(finite.max() - finite.min()) > 1e-12:
            usable.append(col)

    eval_rows: list[dict[str, object]] = []
    for label in usable:
        strict_controls = ["age", "cohort_code", "mean_band_power", "mean_static_plv"]
        for ident in IDENTITY_COLS:
            eval_rows.append(eval_pair(rows, label, ident, []))
            eval_rows.append(eval_pair(rows, label, ident, strict_controls))
    for mode in ("", "age;cohort_code;mean_band_power;mean_static_plv"):
        subset = [r for r in eval_rows if r["controls"] == mode]
        stats.bh_q(subset, "p", "q")

    best_rows = []
    for label in usable:
        strict = [r for r in eval_rows if r["label"] == label and r["controls"]]
        raw = [r for r in eval_rows if r["label"] == label and not r["controls"]]
        best_strict = min(strict, key=lambda r: ffloat(r["p"])) if strict else None
        best_raw = min(raw, key=lambda r: ffloat(r["p"])) if raw else None
        if best_strict and best_raw:
            best_rows.append(
                {
                    "label": label,
                    "n": best_strict["n"],
                    "best_raw_coord": best_raw["identity_coord"],
                    "best_raw_r": best_raw["r"],
                    "best_raw_p": best_raw["p"],
                    "best_raw_q": best_raw.get("q", ""),
                    "best_strict_coord": best_strict["identity_coord"],
                    "best_strict_r": best_strict["r"],
                    "best_strict_p": best_strict["p"],
                    "best_strict_q": best_strict.get("q", ""),
                }
            )
    best_rows.sort(key=lambda r: ffloat(r["best_strict_p"]))

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "identity_label_eval.csv", eval_rows)
    write_csv(OUT / "identity_label_best.csv", best_rows)

    lines = [
        "# LEMON Shared Identity Label Scan",
        "",
        "Frozen object: all-band/all-electrode shared identity coordinates from `lemon_shared_identity_sex_beam.py`.",
        f"Subjects: {len(rows)}",
        f"Usable numeric labels: {len(usable)}",
        "",
        "## Top Strict Label Associations",
        "",
    ]
    for row in best_rows[:40]:
        lines.append(
            "- {label}: n={n}, {coord}, r={r:+.3f}, p={p:.4g}, q={q:.4g}".format(
                label=row["label"],
                n=row["n"],
                coord=row["best_strict_coord"],
                r=ffloat(row["best_strict_r"]),
                p=ffloat(row["best_strict_p"]),
                q=ffloat(row["best_strict_q"]),
            )
        )
    lines.extend(
        [
            "",
            "## Reading Rule",
            "",
            "These are label scans over a frozen identity object, not a new feature search. Strict rows control age, cohort, mean band power, and mean static PLV. Low q values are stronger candidates; isolated low p values should be treated as exploratory.",
        ]
    )
    (OUT / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT / 'SUMMARY.md'}")


if __name__ == "__main__":
    main()
