# Cross-validation of Flint EEG readouts

New analyses dated 5 September 2026 · Simon Velez

Three fixed scalar EEG readouts were evaluated for LPS fluid reasoning and WST vocabulary in all 111 participants with complete inputs. Both leave-one-participant-out and 20 repeats of five-fold cross-validation were run. All predictor scaling, calibration coefficients, and baseline fits use training participants only. The EEG operators and retained per-participant features are unchanged.

The readouts were specified from the published lineage before these runs: raw broad-posterior theta concentration, raw lateral-posterior theta concentration, and the earlier posterior-theta focus ratio. Every specified result is reported; no new best-score search or K selection is performed. The historical discovery of these readouts used this cohort and is not made independent by this new evaluation.

## Leave-one-out results

Each participant is predicted once using the other 110 participants. `EEG only` uses a scalar score plus intercept; its baseline is the training-set mean. `EEG + covariates` adds cohort code, age-bin midpoint, band power, and static phase locking; its baseline uses those same four covariates without the score. Calibration is ordinary least squares on the native outcome scale. The covariates are standardized using the training fold, so small power values are not silently dropped through numerical scale differences.

| Target | Readout | Model | Pearson r | Spearman r | R² | Baseline R² | ΔR² | RMSE |
|---|---|---|---:|---:|---:|---:|---:|---:|
| LPS | Broad theta concentration | EEG only | +0.401 | +0.457 | +0.157 | -0.018 | +0.175 | +3.625 |
| LPS | Broad theta concentration | EEG + covariates | +0.666 | +0.661 | +0.441 | +0.445 | -0.004 | +2.951 |
| LPS | Lateral theta concentration | EEG only | +0.408 | +0.455 | +0.163 | -0.018 | +0.181 | +3.610 |
| LPS | Lateral theta concentration | EEG + covariates | +0.664 | +0.654 | +0.439 | +0.443 | -0.005 | +2.957 |
| LPS | Theta focus ratio | EEG only | +0.429 | +0.459 | +0.182 | -0.018 | +0.201 | +3.569 |
| LPS | Theta focus ratio | EEG + covariates | +0.683 | +0.676 | +0.465 | +0.445 | +0.019 | +2.888 |
| WST | Broad theta concentration | EEG only | -0.364 | -0.269 | -0.036 | -0.018 | -0.018 | +8.665 |
| WST | Broad theta concentration | EEG + covariates | -0.033 | -0.053 | -0.081 | -0.080 | -0.001 | +8.852 |
| WST | Lateral theta concentration | EEG only | -0.381 | -0.283 | -0.037 | -0.018 | -0.019 | +8.669 |
| WST | Lateral theta concentration | EEG + covariates | -0.023 | -0.027 | -0.073 | -0.080 | +0.007 | +8.819 |
| WST | Theta focus ratio | EEG only | -0.562 | -0.445 | -0.040 | -0.018 | -0.022 | +8.681 |
| WST | Theta focus ratio | EEG + covariates | -0.130 | -0.089 | -0.081 | -0.080 | -0.002 | +8.853 |

![Contribution over matched baselines](figures/cv_incremental_performance.png)

The EEG-only reasoning correlations are 0.401–0.429, with R² of 0.157–0.182. Full covariate models give higher total correlations, but most of their predictive performance is already available from the baseline covariates. Adding the broad and lateral concentration scores changes leave-one-out R² by approximately -0.004 and -0.005; adding the focus ratio changes it by +0.019. The full-model correlations must not be reported as the contribution of the EEG score.

The WST predictions are near the sample mean and have negative pooled R² in every leave-one-out model. Negative prediction correlations here do not establish an inverse vocabulary biomarker: leave-one-out mean effects and unstable calibration can produce correlations even when predictions convey little useful variation. Calibration and baseline comparisons are therefore reported beside r.

## Repeated five-fold results

Each of 20 random partitions predicts every participant once from the other four folds. The random seed is 20260905, and the exact folds are retained. Intervals below are 5th–95th percentiles across the 20 partitions, not confidence intervals, independent cohorts, or additional participants.

| Target | Readout | Model | Median Pearson r | Median R² | R² across partitions (5–95%) | Median ΔR² |
|---|---|---|---:|---:|---:|---:|
| LPS | Broad theta concentration | EEG only | +0.405 | +0.157 | [+0.090, +0.180] | +0.175 |
| LPS | Broad theta concentration | EEG + covariates | +0.657 | +0.426 | [+0.353, +0.461] | -0.011 |
| LPS | Lateral theta concentration | EEG only | +0.410 | +0.165 | [+0.105, +0.182] | +0.181 |
| LPS | Lateral theta concentration | EEG + covariates | +0.655 | +0.426 | [+0.347, +0.455] | -0.009 |
| LPS | Theta focus ratio | EEG only | +0.435 | +0.187 | [+0.146, +0.198] | +0.204 |
| LPS | Theta focus ratio | EEG + covariates | +0.677 | +0.455 | [+0.383, +0.479] | +0.013 |
| WST | Broad theta concentration | EEG only | -0.205 | -0.040 | [-0.111, -0.014] | -0.026 |
| WST | Broad theta concentration | EEG + covariates | -0.005 | -0.083 | [-0.122, -0.052] | +0.005 |
| WST | Lateral theta concentration | EEG only | -0.208 | -0.041 | [-0.109, -0.014] | -0.026 |
| WST | Lateral theta concentration | EEG + covariates | +0.007 | -0.076 | [-0.109, -0.045] | +0.017 |
| WST | Theta focus ratio | EEG only | -0.247 | -0.046 | [-0.101, -0.016] | -0.031 |
| WST | Theta focus ratio | EEG + covariates | -0.078 | -0.087 | [-0.137, -0.053] | +0.004 |

![Reasoning stability across all partitions](figures/cv_partition_stability.png)

## Participant predictions

![Held-participant predictions](figures/cv_predictions.png)

Every point is an out-of-fold prediction. The broad concentration score is displayed because it was the original headline readout, not because it won this comparison. The axes use the original LPS and WST values; the diagonal is perfect agreement.

## What was computed

- 12 leave-one-out evaluations: three readouts × two targets × two calibration models.
- 240 repeated five-fold evaluations: the same 12 configurations × 20 partitions.
- 27,972 held-participant predictions, including their matched baseline predictions.
- 24 aggregate summary rows covering both protocols.
- 12 target-perturbation checks verify that changing one held participant’s label leaves that participant’s prediction unchanged.

This is new cross-validation of fixed retained EEG features. It is distinct from the original adjusted partial Spearman result: the new calibration uses raw numerical values with linear terms, while the historical statistic correlates nuisance-residualized ranks. It is also distinct from the five-band phenotype extension, whose participant predictions remain unavailable for a fresh rerun.

## Reproduce

```sh
python scripts/cross_validate_readouts.py
python scripts/write_cv_report.py
python scripts/plot_results.py
```

[Run manifest](results/cross_validation/manifest.json) · [Fold assignments](results/cross_validation/fold_assignments.csv) · [All predictions](results/cross_validation/predictions.csv) · [All run metrics](results/cross_validation/run_metrics.csv) · [Summary](results/cross_validation/summary.csv)

R² is 1 minus pooled out-of-fold squared error divided by the evaluation sample’s centered sum of squares. Baselines use the same denominator. ΔR² is model R² minus the corresponding baseline R²; it is not an in-sample change in explained variance.
