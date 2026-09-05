"""Render both cross-validation tables into the report template.

The CSV supplies every number. The Markdown template holds the explanation so
changes to the report's layout do not require editing the evaluation code.
"""

from pathlib import Path
from string import Template

from data_io import read_csv

ROOT = Path(__file__).resolve().parents[1]
SCORE_NAMES = {
    "theta_concentration_broad": "Broad theta concentration",
    "theta_concentration_lateral": "Lateral theta concentration",
    "theta_focus_ratio": "Theta focus ratio",
}
MODEL_NAMES = {
    "score_only": "EEG only",
    "score_plus_covariates": "EEG + covariates",
}


def index_summaries(rows: list[dict[str, str]]) -> dict:
    indexed = {}
    for row in rows:
        key = (row["target"], row["score"], row["model"], row["protocol"])
        assert key not in indexed, f"Duplicate summary: {key}"
        indexed[key] = row
    return indexed


def format_leave_one_out_table(summaries: dict) -> str:
    lines = [
        "| Target | Readout | Model | Pearson r | Spearman r | R² | Baseline R² | ΔR² | RMSE |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    metrics = ("pearson_r", "spearman_r", "r2", "baseline_r2", "delta_r2", "rmse")
    for target in ("LPS", "WST"):
        for score, score_name in SCORE_NAMES.items():
            for model, model_name in MODEL_NAMES.items():
                row = summaries[target, score, model, "leave_one_out"]
                values = " | ".join(
                    f"{float(row[metric + '_median']):+.3f}" for metric in metrics
                )
                lines.append(f"| {target} | {score_name} | {model_name} | {values} |")
    return "\n".join(lines)


def format_repeated_folds_table(summaries: dict) -> str:
    lines = [
        "| Target | Readout | Model | Median Pearson r | Median R² | R² across partitions (5–95%) | Median ΔR² |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for target in ("LPS", "WST"):
        for score, score_name in SCORE_NAMES.items():
            for model, model_name in MODEL_NAMES.items():
                row = summaries[target, score, model, "repeated_5fold"]
                partition_range = (
                    f"[{float(row['r2_q05']):+.3f}, {float(row['r2_q95']):+.3f}]"
                )
                values = [
                    target,
                    score_name,
                    model_name,
                    f"{float(row['pearson_r_median']):+.3f}",
                    f"{float(row['r2_median']):+.3f}",
                    partition_range,
                    f"{float(row['delta_r2_median']):+.3f}",
                ]
                lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    rows = read_csv(ROOT / "results/cross_validation/summary.csv")
    summaries = index_summaries(rows)
    template = Template(
        (Path(__file__).parent / "templates/cross_validation.md").read_text()
    )
    report = template.substitute(
        leave_one_out_table=format_leave_one_out_table(summaries),
        repeated_folds_table=format_repeated_folds_table(summaries),
    )
    (ROOT / "CROSS_VALIDATION.md").write_text(report)
    print(f"Wrote CROSS_VALIDATION.md from all {len(rows)} summary rows.")


if __name__ == "__main__":
    main()
