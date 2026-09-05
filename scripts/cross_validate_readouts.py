"""Cross-validate scalar calibrations of three retained EEG scores.

The EEG operators and historical feature choices are fixed. For each score,
calibrate LPS and WST with and without covariates, using identical participant
folds. Scaling, least-squares fits, and baselines use training participants only.
"""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from data_io import numeric_column, read_csv, write_csv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/retained"
OUT = ROOT / "results/cross_validation"
REPEATS = 20
SEED = 20260905
TARGETS = ("LPS", "WST")
MODELS = ("score_only", "score_plus_covariates")
PROTOCOLS = ("leave_one_out", "repeated_5fold")
METRICS = (
    "pearson_r",
    "spearman_r",
    "r2",
    "baseline_r2",
    "delta_r2",
    "rmse",
    "baseline_rmse",
)


@dataclass(frozen=True)
class Readout:
    """A fixed score column and its matching nuisance-control columns."""

    name: str
    rows: list[dict[str, str]]
    feature: str
    controls: list[str]


@dataclass(frozen=True)
class FoldPlan:
    """One whole-cohort evaluation; every participant is held out once."""

    protocol: str
    repeat: int
    test_folds: list[list[int]]


def load_readouts() -> tuple[list[str], list[Readout]]:
    spatial = read_csv(DATA / "spatial_subjects_full_full.csv")
    benchmark = [
        row
        for row in read_csv(DATA / "bfocus_features_full_grid24.csv")
        if (row["condition"], row["band"], row["channel_set"])
        == ("EO", "theta", "posterior")
    ]
    readouts = [
        Readout(
            "theta_concentration_broad",
            spatial,
            "raw_posterior_broad__full_score",
            [
                "cohort_code",
                "age",
                "raw_posterior_broad__theta_power",
                "raw_posterior_broad__theta_static_plv",
            ],
        ),
        Readout(
            "theta_concentration_lateral",
            spatial,
            "raw_lateral_posterior__full_score",
            [
                "cohort_code",
                "age",
                "raw_lateral_posterior__theta_power",
                "raw_lateral_posterior__theta_static_plv",
            ],
        ),
        Readout(
            "theta_focus_ratio",
            benchmark,
            "bfocus_ratio",
            ["cohort_code", "age", "band_power", "static_plv"],
        ),
    ]
    participant_ids = sorted(row["sid"] for row in spatial)
    assert len(participant_ids) == 111 and len(set(participant_ids)) == 111
    return participant_ids, readouts


def make_fold_plans(participant_count: int) -> list[FoldPlan]:
    plans = [FoldPlan("leave_one_out", 0, [[i] for i in range(participant_count)])]
    random = np.random.default_rng(SEED)
    for repeat in range(REPEATS):
        folds = np.array_split(random.permutation(participant_count), 5)
        plans.append(
            FoldPlan("repeated_5fold", repeat, [fold.tolist() for fold in folds])
        )
    return plans


def fold_records(plans: list[FoldPlan], participant_ids: list[str]) -> list[dict]:
    return [
        dict(
            protocol=plan.protocol,
            repeat=plan.repeat,
            fold=fold,
            sid=participant_ids[index],
        )
        for plan in plans
        for fold, test_indices in enumerate(plan.test_folds)
        for index in test_indices
    ]


def fit_predict(
    predictors: np.ndarray,
    outcome: np.ndarray,
    train_indices: np.ndarray,
    test_indices: np.ndarray,
) -> np.ndarray:
    """Fit an intercept and linear coefficients, scaling from training data.

    Zero predictor columns select the training-mean baseline. Scaling prevents
    small EEG power columns from disappearing in the least-squares rank cutoff.
    The 1e-12 scale threshold is preserved from the published implementation.
    """
    if predictors.shape[1] == 0:
        return np.full(len(test_indices), outcome[train_indices].mean())
    mean = predictors[train_indices].mean(0)
    scale = predictors[train_indices].std(0)
    scale = np.where(scale > 1e-12, scale, 1)
    train_design = np.column_stack(
        [np.ones(len(train_indices)), (predictors[train_indices] - mean) / scale]
    )
    test_design = np.column_stack(
        [np.ones(len(test_indices)), (predictors[test_indices] - mean) / scale]
    )
    coefficients = np.linalg.lstsq(train_design, outcome[train_indices], rcond=None)[0]
    return test_design @ coefficients


def prediction_metrics(
    outcome: np.ndarray, prediction: np.ndarray, baseline: np.ndarray
) -> dict:
    """Pool held-out errors across participants using one shared R² denominator."""
    total_sum_squares = float(np.sum((outcome - outcome.mean()) ** 2))
    r2 = 1 - float(np.sum((outcome - prediction) ** 2)) / total_sum_squares
    baseline_r2 = 1 - float(np.sum((outcome - baseline) ** 2)) / total_sum_squares
    return dict(
        pearson_r=float(np.corrcoef(outcome, prediction)[0, 1]),
        spearman_r=float(spearmanr(outcome, prediction).statistic),
        r2=r2,
        baseline_r2=baseline_r2,
        delta_r2=r2 - baseline_r2,
        rmse=float(np.sqrt(np.mean((outcome - prediction) ** 2))),
        baseline_rmse=float(np.sqrt(np.mean((outcome - baseline) ** 2))),
    )


def check_held_label_invariance(predictors: np.ndarray, outcome: np.ndarray) -> None:
    """Changing participant zero's held-out label must not change its prediction."""
    train_indices = np.arange(1, len(outcome))
    test_indices = np.array([0])
    perturbed = outcome.copy()
    perturbed[0] += 10000
    np.testing.assert_array_equal(
        fit_predict(predictors, outcome, train_indices, test_indices),
        fit_predict(predictors, perturbed, train_indices, test_indices),
    )


def evaluate_plan(
    predictors: np.ndarray,
    baseline_predictors: np.ndarray,
    outcome: np.ndarray,
    participant_ids: list[str],
    plan: FoldPlan,
    configuration: dict[str, str],
) -> tuple[list[dict], dict]:
    """Fit each fold and return per-person predictions plus pooled metrics."""
    participant_count = len(participant_ids)
    all_indices = np.arange(participant_count)
    prediction = np.full(participant_count, np.nan)
    baseline = np.full(participant_count, np.nan)
    records = []
    for fold, held_indices in enumerate(plan.test_folds):
        test_indices = np.array(held_indices)
        train_indices = np.setdiff1d(all_indices, test_indices)
        prediction[test_indices] = fit_predict(
            predictors, outcome, train_indices, test_indices
        )
        baseline[test_indices] = fit_predict(
            baseline_predictors, outcome, train_indices, test_indices
        )
        for index in test_indices:
            records.append(
                dict(
                    **configuration,
                    protocol=plan.protocol,
                    repeat=plan.repeat,
                    fold=fold,
                    sid=participant_ids[index],
                    observed=outcome[index],
                    prediction=prediction[index],
                    baseline_prediction=baseline[index],
                )
            )
    assert np.isfinite(prediction).all() and np.isfinite(baseline).all()
    metrics = dict(
        **configuration,
        protocol=plan.protocol,
        repeat=plan.repeat,
        n=participant_count,
        **prediction_metrics(outcome, prediction, baseline),
    )
    return records, metrics


def evaluate_readouts(
    participant_ids: list[str], readouts: list[Readout], plans: list[FoldPlan]
):
    predictions = []
    run_metrics = []
    invariance_checks = 0
    for readout in readouts:
        by_id = {row["sid"]: row for row in readout.rows}
        assert len(by_id) == len(
            readout.rows
        ), f"Duplicate participants in {readout.name}"
        assert set(by_id) == set(participant_ids)
        rows = [by_id[participant_id] for participant_id in participant_ids]
        score = numeric_column(rows, readout.feature)[:, None]
        covariates = np.column_stack(
            [numeric_column(rows, key) for key in readout.controls]
        )
        for target in TARGETS:
            outcome = numeric_column(rows, target)
            for model in MODELS:
                baseline_predictors = (
                    covariates
                    if model == "score_plus_covariates"
                    else np.empty((len(rows), 0))
                )
                predictors = np.column_stack([baseline_predictors, score])
                assert np.isfinite(predictors).all() and np.isfinite(outcome).all()
                check_held_label_invariance(predictors, outcome)
                invariance_checks += 1
                configuration = dict(score=readout.name, target=target, model=model)
                for plan in plans:
                    records, metrics = evaluate_plan(
                        predictors,
                        baseline_predictors,
                        outcome,
                        participant_ids,
                        plan,
                        configuration,
                    )
                    predictions.extend(records)
                    run_metrics.append(metrics)
    return predictions, run_metrics, invariance_checks


def summarize_runs(
    run_metrics: list[dict], readouts: list[Readout], participant_count: int
) -> list[dict]:
    summaries = []
    for readout in readouts:
        for target in TARGETS:
            for model in MODELS:
                for protocol in PROTOCOLS:
                    runs = [
                        row
                        for row in run_metrics
                        if (row["score"], row["target"], row["model"], row["protocol"])
                        == (readout.name, target, model, protocol)
                    ]
                    summary = dict(
                        score=readout.name,
                        target=target,
                        model=model,
                        protocol=protocol,
                        n=participant_count,
                        repeats=len(runs),
                    )
                    for metric in METRICS:
                        values = [run[metric] for run in runs]
                        summary[metric + "_median"] = float(np.median(values))
                        summary[metric + "_q05"] = float(np.quantile(values, 0.05))
                        summary[metric + "_q95"] = float(np.quantile(values, 0.95))
                    summaries.append(summary)
    return summaries


def build_manifest(
    readouts,
    participant_count,
    prediction_count,
    run_count,
    summary_count,
    invariance_checks,
):
    return {
        "scope": "New cross-validation from retained fixed participant EEG scores; no EEG extraction rerun or external cohort",
        "subjects": participant_count,
        "targets": ["LPS", "WST"],
        "scores": [
            dict(name=readout.name, feature=readout.feature, controls=readout.controls)
            for readout in readouts
        ],
        "protocols": {
            "leave_one_out": 111,
            "repeated_5fold": {"folds": 5, "repeats": REPEATS, "seed": SEED},
        },
        "models": [
            "score_only versus training-mean baseline",
            "score_plus_covariates versus covariate-only baseline",
        ],
        "calibration": "ordinary least squares on training subjects; linear scalar readout with intercept; scaling fitted in training fold",
        "selection": "scores fixed from published historical lineage before this run; no choice of best model by outcome; historical score discovery not nested",
        "missing_data": "all inputs and targets finite; no target imputation",
        "target_perturbation_invariance_checks": invariance_checks,
        "prediction_rows": prediction_count,
        "run_metric_rows": run_count,
        "summary_rows": summary_count,
        "r2_definition": "1 - pooled out-of-fold SSE / full evaluation sample centered SST; baseline_r2 uses same denominator",
        "repeat_intervals": "5th to 95th percentiles across random partitions; not confidence intervals or new participants",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    participant_ids, readouts = load_readouts()
    plans = make_fold_plans(len(participant_ids))
    write_csv(OUT / "fold_assignments.csv", fold_records(plans, participant_ids))
    predictions, runs, checks = evaluate_readouts(participant_ids, readouts, plans)
    summaries = summarize_runs(runs, readouts, len(participant_ids))
    write_csv(OUT / "predictions.csv", predictions)
    write_csv(OUT / "run_metrics.csv", runs)
    write_csv(OUT / "summary.csv", summaries)
    manifest = build_manifest(
        readouts,
        len(participant_ids),
        len(predictions),
        len(runs),
        len(summaries),
        checks,
    )
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    for summary in summaries:
        if summary["protocol"] == "leave_one_out":
            print(
                summary["score"],
                summary["target"],
                summary["model"],
                "r",
                round(summary["pearson_r_median"], 3),
                "R2",
                round(summary["r2_median"], 3),
                "deltaR2",
                round(summary["delta_r2_median"], 3),
            )


if __name__ == "__main__":
    main()
