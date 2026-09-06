"""Recompute the published partial-rank associations from participant tables.

The analysis uses no raw EEG. It reproduces two result families, checks
within-recording reliability, and bootstraps the selected broad theta score.
"""

import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.stats import rankdata, t

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from data_io import numeric_column, read_csv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/retained"
OUT = ROOT / "results/verification"
RELATIVE_TOLERANCE = 1e-9
ABSOLUTE_TOLERANCE = 1e-10
BOOTSTRAP_REPLICATES = 3000
BOOTSTRAP_SEED = 20260905


def residual_ranks(values: np.ndarray, controls: list[np.ndarray]) -> np.ndarray:
    """Remove an intercept and ranked nuisance controls from ranked values."""
    design = np.column_stack(
        [np.ones(len(values))] + [rankdata(control) for control in controls]
    )
    ranked = rankdata(values)
    coefficients = np.linalg.lstsq(design, ranked, rcond=None)[0]
    return ranked - design @ coefficients


def partial_spearman(
    score: np.ndarray, outcome: np.ndarray, controls: list[np.ndarray]
) -> tuple[float, float, int]:
    """Return partial rank r, its historical t-approximation p, and complete n."""
    complete = np.isfinite(score) & np.isfinite(outcome)
    for control in controls:
        complete &= np.isfinite(control)
    score, outcome = score[complete], outcome[complete]
    controls = [control[complete] for control in controls]
    score_residuals = residual_ranks(score, controls)
    outcome_residuals = residual_ranks(outcome, controls)
    correlation = float(np.corrcoef(score_residuals, outcome_residuals)[0, 1])
    degrees_of_freedom = len(score) - len(controls) - 2
    statistic = abs(correlation) * np.sqrt(
        degrees_of_freedom / max(1 - correlation * correlation, 1e-9)
    )
    p_value = float(2 * t.sf(statistic, degrees_of_freedom))
    return correlation, p_value, len(score)


def benjamini_hochberg(p_values) -> np.ndarray:
    """Adjust one table's p-values and restore their original row order."""
    p_values = np.array(p_values)
    order = np.argsort(p_values)
    scaled = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    sorted_q_values = np.minimum.accumulate(scaled[::-1])[::-1]
    q_values = np.empty(len(p_values))
    q_values[order] = np.minimum(sorted_q_values, 1)
    return q_values


def select_feature_rows(features, evaluation, benchmark: bool) -> list[dict]:
    if not benchmark:
        return features
    group = [
        row
        for row in features
        if all(
            row[key] == evaluation[key] for key in ("condition", "band", "channel_set")
        )
    ]
    if evaluation["analysis"].startswith("young"):
        group = [row for row in group if row["cohort"] == "young"]
    return group


def compute_associations(features, evaluation) -> dict:
    score = numeric_column(features, evaluation["feature"])
    reasoning = numeric_column(features, "LPS")
    vocabulary = numeric_column(features, "WST")
    controls = [
        numeric_column(features, key) for key in evaluation["controls"].split(";")
    ]
    reasoning_r, reasoning_p, count = partial_spearman(score, reasoning, controls)
    vocabulary_r, vocabulary_p, _ = partial_spearman(score, vocabulary, controls)
    adjusted_r, adjusted_p, _ = partial_spearman(
        score, reasoning, controls + [vocabulary]
    )
    return dict(
        n=count,
        r_lps=reasoning_r,
        p_lps=reasoning_p,
        r_wst=vocabulary_r,
        p_wst=vocabulary_p,
        r_lps_plus_wst=adjusted_r,
        p_lps_plus_wst=adjusted_p,
    )


def verify_family(features, evaluations, benchmark: bool = False) -> dict:
    computed = [
        compute_associations(
            select_feature_rows(features, evaluation, benchmark), evaluation
        )
        for evaluation in evaluations
    ]
    for p_key, q_key in (("p_lps", "q_lps"), ("p_lps_plus_wst", "q_lps_plus_wst")):
        if q_key in evaluations[0]:
            q_values = benjamini_hochberg([row[p_key] for row in computed])
            for row, q_value in zip(computed, q_values):
                row[q_key] = float(q_value)
    errors = []
    for expected, actual in zip(evaluations, computed):
        for key, value in actual.items():
            reference = float(expected[key])
            if not np.isclose(
                value,
                reference,
                rtol=RELATIVE_TOLERANCE,
                atol=ABSOLUTE_TOLERANCE,
                equal_nan=True,
            ):
                raise AssertionError(
                    (expected["feature"], expected["analysis"], key, value, reference)
                )
            errors.append(abs(value - reference))
    return dict(
        rows=len(evaluations),
        values_compared=len(errors),
        max_absolute_difference=max(errors),
    )


def verify_reliability(participants, reliability_rows) -> int:
    for row in reliability_rows:
        first, second = (
            ("odd", "even") if row["split"] == "odd_even" else ("first", "second")
        )
        prefix = row["transform"] + "_" + row["aperture"] + "__"
        correlation, _, count = partial_spearman(
            numeric_column(participants, prefix + first + "_score"),
            numeric_column(participants, prefix + second + "_score"),
            [],
        )
        assert count == int(row["n"])
        np.testing.assert_allclose(
            [correlation, 2 * correlation / (1 + correlation)],
            [float(row["split_r"]), float(row["spearman_brown"])],
            rtol=RELATIVE_TOLERANCE,
            atol=ABSOLUTE_TOLERANCE,
        )
    return len(reliability_rows)


def bootstrap_association(participants, result) -> dict:
    """Resample participants, reranking and refitting controls in each sample."""
    score = numeric_column(participants, result["feature"])
    outcome = numeric_column(participants, "LPS")
    controls = [
        numeric_column(participants, key) for key in result["controls"].split(";")
    ]
    random = np.random.default_rng(BOOTSTRAP_SEED)
    correlations = []
    for _ in range(BOOTSTRAP_REPLICATES):
        indices = random.integers(0, len(score), len(score))
        correlation, _, _ = partial_spearman(
            score[indices], outcome[indices], [control[indices] for control in controls]
        )
        correlations.append(correlation)
    return {
        "replicates": BOOTSTRAP_REPLICATES,
        "seed": BOOTSTRAP_SEED,
        "percentile_95_interval": np.quantile(correlations, [0.025, 0.975]).tolist(),
        "interpretation": "Conditional on historical score choice; not adjusted for discovery selection",
    }


def plot_associations(participants, evaluations, result) -> None:
    """Plot residual ranks and the retained spatial ablations."""
    controls = [
        numeric_column(participants, key) for key in result["controls"].split(";")
    ]
    score_residuals = residual_ranks(
        numeric_column(participants, result["feature"]), controls
    )
    reasoning_residuals = residual_ranks(numeric_column(participants, "LPS"), controls)
    vocabulary_residuals = residual_ranks(numeric_column(participants, "WST"), controls)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.hashsalt": "flint-v030",
        }
    )
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.1))
    colors = ["#277a83" if r["cohort"] == "young" else "#b86143" for r in participants]
    for axis, target, label, r in [
        (axes[0], reasoning_residuals, "Fluid reasoning (LPS)", result["r_lps"]),
        (axes[1], vocabulary_residuals, "Vocabulary (WST)", result["r_wst"]),
    ]:
        axis.scatter(
            score_residuals, target, c=colors, s=22, alpha=0.75, edgecolors="none"
        )
        slope, intercept = np.polyfit(score_residuals, target, 1)
        domain = np.array([score_residuals.min(), score_residuals.max()])
        axis.plot(domain, slope * domain + intercept, color="#273544", lw=1.5)
        axis.set(
            xlabel="Theta concentration: residual rank",
            ylabel=label + ": residual rank",
            title=f"{label} | partial r = {float(r):.3f}",
        )
    labels = [
        "Raw",
        "Smooth .35",
        "Smooth .50",
        "High-boost .10",
        "High-boost .25",
        "High-boost .50",
        "High-boost .75",
    ]
    transforms = [
        "raw",
        "smooth35",
        "smooth50",
        "deblur10",
        "deblur25",
        "deblur50",
        "deblur75",
    ]
    values = [
        float(
            next(
                r
                for r in evaluations
                if r["feature"] == f"{transform}_posterior_broad__full_score"
                and r["analysis"] == "strict"
            )["r_lps"]
        )
        for transform in transforms
    ]
    axes[2].barh(labels[::-1], values[::-1], color=["#aaa"] * 6 + ["#277a83"])
    axes[2].set(
        xlim=(0, 0.5),
        xlabel="Partial rank correlation with LPS",
        title="Retained spatial ablations",
    )
    figure.suptitle(
        "Eyes-open theta lagged-field concentration and fluid reasoning | 111 LEMON participants",
        fontsize=13,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    for extension in ("png", "svg"):
        figure.savefig(
            ROOT / f"figures/correlate.{extension}",
            dpi=180,
            bbox_inches="tight",
            metadata={"Date": None} if extension == "svg" else None,
        )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    participants = read_csv(DATA / "spatial_subjects_full_full.csv")
    benchmark_features = read_csv(DATA / "bfocus_features_full_grid24.csv")
    spatial_evaluations = read_csv(DATA / "spatial_eval_full_full.csv")
    benchmark_evaluations = read_csv(DATA / "bfocus_eval_full_grid24.csv")
    assert len(participants) == 111 and len({row["sid"] for row in participants}) == 111
    report = {
        "scope": "Statistical recomputation from retained participant features; no raw-EEG rerun",
        "spatial": verify_family(participants, spatial_evaluations),
        "benchmark": verify_family(
            benchmark_features, benchmark_evaluations, benchmark=True
        ),
    }
    report["reliability_rows_verified"] = verify_reliability(
        participants, read_csv(DATA / "spatial_reliability_full_full.csv")
    )
    result = next(
        row
        for row in spatial_evaluations
        if row["feature"] == "raw_posterior_broad__full_score"
        and row["analysis"] == "strict"
    )
    report["headline"] = result
    report["headline_bootstrap"] = bootstrap_association(participants, result)
    plot_associations(participants, spatial_evaluations, result)
    report["passed"] = True
    (OUT / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
