"""Publication figures from complete retained results and new CV outputs."""

import json
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from data_io import read_csv

ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures"
DOMAIN_NAMES = {
    "blood_chemistry": "Blood chemistry",
    "blood_pressure": "Blood pressure",
    "body_size": "Body measurements",
    "cognition": "Cognition",
    "personality_affect": "Personality / affect",
}
COLORS = ["#087e8b", "#c75b39", "#5167ad", "#8b5aa6", "#ad8325"]
SCORES = {
    "theta_concentration_broad": "Broad theta concentration",
    "theta_concentration_lateral": "Lateral theta concentration",
    "theta_focus_ratio": "Theta focus ratio",
}
COORDINATE_COUNTS = [3, 6, 10, 20]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titleweight": "bold",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.hashsalt": "flint-v030",
    }
)


def save_figure(figure, name: str) -> None:
    """Write matching PNG and SVG files, then release the figure."""
    for extension in ("png", "svg"):
        figure.savefig(
            FIGURE_DIR / f"{name}.{extension}",
            dpi=170,
            bbox_inches="tight",
            metadata={"Date": None} if extension == "svg" else None,
        )
    plt.close(figure)


def add_figure_heading(figure, title: str, subtitle: str) -> None:
    figure.suptitle(title, x=0.07, ha="left", fontsize=17, y=0.995)
    figure.text(0.07, 0.918, subtitle, ha="left", fontsize=10, color="#555")
    figure.tight_layout(rect=(0, 0, 1, 0.90))


def plot_phenotype_distribution(rows) -> None:
    """Plot all 291 fields at six coordinates, retaining signed correlations."""
    figure, axis = plt.subplots(figsize=(10, 5.8))
    random = np.random.default_rng(20260905)
    for i, (domain, title) in enumerate(DOMAIN_NAMES.items()):
        group = [
            row
            for row in rows
            if row["domain"] == domain and int(row["k_identity_dims"]) == 6
        ]
        correlations = np.array([float(row["cv_corr"]) for row in group])
        axis.scatter(
            correlations,
            i + random.uniform(-0.24, 0.24, len(correlations)),
            s=20,
            alpha=0.67,
            c=COLORS[i],
            edgecolors="none",
        )
        axis.plot(
            [np.median(correlations)] * 2, [i - 0.34, i + 0.34], c="#172c38", lw=2
        )
    axis.axvline(0, color="#bbb", lw=1)
    axis.set(
        yticks=range(5),
        yticklabels=[
            f'{v} (n={sum(row["domain"]==domain and int(row["k_identity_dims"])==6 for row in rows)})'
            for domain, v in DOMAIN_NAMES.items()
        ],
        xlabel="Reported leave-one-out Pearson r",
        xlim=(-0.65, 0.65),
    )
    axis.invert_yaxis()
    add_figure_heading(
        figure,
        "The full phenotype map",
        "All 291 recorded fields at K = 6 · ticks mark medians · original mean-imputed target scoring",
    )
    save_figure(figure, "phenotype_distribution")


def plot_domain_generalization(domains) -> None:
    """Compare pooled domain R² at every retained coordinate count."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for (domain, title), color in zip(DOMAIN_NAMES.items(), COLORS):
        rows_by_k = {
            int(row["k_identity_dims"]): row
            for row in domains
            if row["domain"] == domain
        }
        for axis, key in zip(axes, ["loocv_r2", "in_sample_r2"]):
            axis.plot(
                COORDINATE_COUNTS,
                [float(rows_by_k[k][key]) for k in COORDINATE_COUNTS],
                marker="o",
                label=title,
                color=color,
            )
            axis.axhline(0, color="#bbb", lw=0.6)
            axis.set(
                xticks=COORDINATE_COUNTS,
                xlabel="Retained identity coordinates K",
                ylabel="Pooled standardized R²",
            )
    axes[0].set_title("Leave-one-out")
    axes[1].set_title("In sample")
    axes[1].legend(frameon=False, fontsize=9)
    add_figure_heading(
        figure,
        "Generalization across whole domains",
        "Five domains · every tested K · domain R² includes every field, not only favorable rows",
    )
    save_figure(figure, "domain_generalization")


def plot_domain_heatmaps(rows) -> None:
    """Draw one complete label-by-coordinate matrix per domain."""
    for domain, title in DOMAIN_NAMES.items():
        labels = sorted({row["label"] for row in rows if row["domain"] == domain})
        lookup = {
            (row["label"], int(row["k_identity_dims"])): float(row["cv_corr"])
            for row in rows
            if row["domain"] == domain
        }
        matrix = np.array(
            [[lookup[(label, k)] for k in COORDINATE_COUNTS] for label in labels]
        )
        figure, axis = plt.subplots(figsize=(12, max(3.8, 0.235 * len(labels) + 1.6)))
        heatmap = axis.imshow(
            matrix,
            aspect="auto",
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-0.65, vcenter=0, vmax=0.65),
        )
        axis.set(
            xticks=range(4),
            xticklabels=COORDINATE_COUNTS,
            yticks=range(len(labels)),
            yticklabels=[" / ".join(label.split("__")[-2:]) for label in labels],
            xlabel="Retained identity coordinates K",
        )
        axis.tick_params(axis="y", labelsize=8)
        axis.set_title(
            f"{title}: all {len(labels)} fields\nReported leave-one-out Pearson r · source-name order",
            loc="left",
            pad=13,
        )
        if len(labels) <= 40:
            for i in range(len(labels)):
                for j in range(4):
                    axis.text(
                        j,
                        i,
                        f"{matrix[i,j]:+.2f}",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if abs(matrix[i, j]) > 0.42 else "#111",
                    )
        figure.colorbar(
            heatmap,
            ax=axis,
            pad=0.025,
            shrink=min(1, 5 / figure.get_figheight()),
            label="Signed prediction correlation",
        )
        figure.tight_layout()
        save_figure(figure, f"atlas_{domain}")


def plot_dimension_curves(rows) -> None:
    """Plot the six previously published outcomes at all four K settings."""
    chosen = [
        ("Waist", "__Waist_cm"),
        ("Systolic pressure (BP2, left)", "__BP2_left_systole"),
        ("LPS reasoning", "__LPS_1"),
        ("CKD-EPI", "__CKDEPI_in_ml_min_1.73m"),
        ("HbA1c (%)", "__HBA1C_in_%"),
        ("Optimism", "__LOT_Optimism"),
    ]
    figure, axes = plt.subplots(2, 3, figsize=(12, 8))
    for axis, (title, suffix) in zip(axes.flat, chosen):
        matches = [row for row in rows if row["label"].endswith(suffix)]
        assert len(matches) == 4, (title, len(matches))
        rows_by_k = {int(row["k_identity_dims"]): row for row in matches}
        axis.plot(
            COORDINATE_COUNTS,
            [float(rows_by_k[k]["cv_corr"]) for k in COORDINATE_COUNTS],
            "-o",
            c="#087e8b",
            label="Leave-one-out",
        )
        axis.plot(
            COORDINATE_COUNTS,
            [float(rows_by_k[k]["in_sample_corr"]) for k in COORDINATE_COUNTS],
            "--o",
            c="#bd7048",
            label="In sample",
        )
        axis.set(
            title=title,
            xticks=COORDINATE_COUNTS,
            xlabel="K",
            ylabel="Pearson r",
            ylim=(0, 0.85),
        )
        axis.axhline(0, c="#aaa", lw=0.5)
    axes[0, 0].legend(frameon=False, fontsize=9)
    add_figure_heading(
        figure,
        "How the readouts change with dimensionality",
        "Six previously reported examples · all four K settings · 111 cohort rows with original target imputation",
    )
    save_figure(figure, "phenotype_dimension_curves")


def plot_incremental_performance(cv_summaries) -> None:
    """Compare leave-one-out gains with repeated-fold partition ranges."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.8))
    for axis, target in zip(axes, ["LPS", "WST"]):
        for offset, model, color, label in [
            (-0.13, "score_only", "#087e8b", "EEG score vs mean"),
            (
                0.13,
                "score_plus_covariates",
                "#bd7048",
                "EEG + covariates vs covariates",
            ),
        ]:
            for i, name in enumerate(SCORES):
                leave_one_out = next(
                    row
                    for row in cv_summaries
                    if (row["score"], row["target"], row["model"], row["protocol"])
                    == (name, target, model, "leave_one_out")
                )
                repeated_folds = next(
                    row
                    for row in cv_summaries
                    if (row["score"], row["target"], row["model"], row["protocol"])
                    == (name, target, model, "repeated_5fold")
                )
                y = i + offset
                axis.plot(
                    [
                        float(repeated_folds["delta_r2_q05"]),
                        float(repeated_folds["delta_r2_q95"]),
                    ],
                    [y, y],
                    c=color,
                    lw=4,
                    alpha=0.35,
                )
                axis.scatter(
                    float(leave_one_out["delta_r2_median"]),
                    y,
                    c=color,
                    s=55,
                    label=label if i == 0 else None,
                    zorder=4,
                )
        axis.axvline(0, c="#555", lw=1)
        axis.set(
            yticks=range(3),
            yticklabels=list(SCORES.values()),
            xlabel="Improvement in out-of-fold R² over matched baseline",
            title=target,
        )
        axis.invert_yaxis()
        axis.legend(
            frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18)
        )
    add_figure_heading(
        figure,
        "What does the EEG readout add?",
        "New analyses · 111 participants · points: leave-one-out · bars: 5–95% across 20 five-fold partitions (not CIs)",
    )
    save_figure(figure, "cv_incremental_performance")


def plot_participant_predictions(predictions, cv_summaries) -> None:
    """Show broad-score predictions, observed outcomes, and calibration error."""
    figure, axes = plt.subplots(2, 2, figsize=(10, 9))
    for axis, (target, model) in zip(
        axes.flat,
        [
            ("LPS", "score_only"),
            ("LPS", "score_plus_covariates"),
            ("WST", "score_only"),
            ("WST", "score_plus_covariates"),
        ],
    ):
        group = [
            row
            for row in predictions
            if (row["score"], row["target"], row["model"], row["protocol"])
            == ("theta_concentration_broad", target, model, "leave_one_out")
        ]
        x = np.array([float(row["observed"]) for row in group])
        y = np.array([float(row["prediction"]) for row in group])
        row = next(
            row
            for row in cv_summaries
            if (row["score"], row["target"], row["model"], row["protocol"])
            == ("theta_concentration_broad", target, model, "leave_one_out")
        )
        axis.scatter(x, y, c="#087e8b", s=26, alpha=0.7, edgecolors="none")
        bounds = [min(x.min(), y.min()), max(x.max(), y.max())]
        axis.plot(bounds, bounds, "--", c="#bbb", lw=1)
        axis.set(
            xlabel=f"Observed {target}",
            ylabel=f"Out-of-fold predicted {target}",
            title=f"{target}: "
            + ("EEG only" if model == "score_only" else "EEG + covariates"),
        )
        axis.text(
            0.03,
            0.97,
            f"r = {float(row['pearson_r_median']):+.3f}\nR² = {float(row['r2_median']):+.3f}\nBaseline R² = {float(row['baseline_r2_median']):+.3f}",
            transform=axis.transAxes,
            va="top",
            fontsize=9,
        )
    add_figure_heading(
        figure,
        "Participant-level predictions",
        "Broad posterior theta concentration · leave-one-out calibration · identity line shown · scales differ by target",
    )
    save_figure(figure, "cv_predictions")


def plot_partition_stability(runs) -> None:
    """Show each complete five-fold reasoning run and its matched baseline."""
    figure, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    for axis, model in zip(axes, ["score_only", "score_plus_covariates"]):
        for i, (name, label) in enumerate(SCORES.items()):
            group = [
                row
                for row in runs
                if (row["score"], row["target"], row["model"], row["protocol"])
                == (name, "LPS", model, "repeated_5fold")
            ]
            values = [float(row["r2"]) for row in group]
            axis.scatter(
                np.repeat(i, 20) + np.linspace(-0.14, 0.14, 20),
                values,
                s=30,
                c=COLORS[i],
                alpha=0.8,
            )
            baseline = np.median([float(row["baseline_r2"]) for row in group])
            axis.plot(
                [i - 0.22, i + 0.22],
                [baseline, baseline],
                c="#333",
                lw=2,
                label="Median matched baseline" if i == 0 else None,
            )
        axis.set(
            xticks=range(3),
            xticklabels=["Broad", "Lateral", "Focus ratio"],
            ylabel="Out-of-fold R²",
            title="EEG only" if model == "score_only" else "EEG + covariates",
        )
        axis.legend(frameon=False, fontsize=8)
    add_figure_heading(
        figure,
        "Reasoning prediction across 20 data partitions",
        "Five-fold CV · all 111 participants per partition · every run shown · black marks: matched baselines",
    )
    save_figure(figure, "cv_partition_stability")


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rows = read_csv(ROOT / "results/context/phenotypes/domain_label_reconstruction.csv")
    domains = read_csv(ROOT / "results/context/phenotypes/domain_recovery_summary.csv")
    cv_summaries = read_csv(ROOT / "results/cross_validation/summary.csv")
    runs = read_csv(ROOT / "results/cross_validation/run_metrics.csv")
    predictions = read_csv(ROOT / "results/cross_validation/predictions.csv")
    assert len(rows) == 1164
    plot_phenotype_distribution(rows)
    plot_domain_generalization(domains)
    plot_domain_heatmaps(rows)
    plot_dimension_curves(rows)
    plot_incremental_performance(cv_summaries)
    plot_participant_predictions(predictions, cv_summaries)
    plot_partition_stability(runs)
    outputs = [
        "phenotype_distribution",
        "domain_generalization",
        *[f"atlas_{domain}" for domain in DOMAIN_NAMES],
        "phenotype_dimension_curves",
        "cv_incremental_performance",
        "cv_predictions",
        "cv_partition_stability",
    ]
    (FIGURE_DIR / "figure_manifest.json").write_text(
        json.dumps(
            {
                "new_figures": outputs,
                "count": len(outputs),
                "formats": ["png", "svg"],
                "phenotype_source": "retained aggregate outputs",
                "cv_source": "new participant-level cross-validation",
                "all_phenotype_rows_plotted": len(rows),
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(outputs))


if __name__ == "__main__":
    main()
