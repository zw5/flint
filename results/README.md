# Result tables

All tables are per-participant or per-analysis CSV files. Participant IDs use
the LEMON `sid` convention. Nothing here is raw EEG.

## `retained/` — original participant features and evaluations

Copied unchanged from the research workspace. These are the inputs to
`scripts/verify.py` and `scripts/cross_validate_readouts.py`.

| File | Rows | Contents |
|---|---:|---|
| `spatial_subjects_full_full.csv` | 111 | One row per participant: LPS, WST, age, cohort, and every spatial-transform × aperture score, with odd/even and first/second split scores. |
| `spatial_eval_full_full.csv` | 112 | Partial Spearman associations for each spatial score under each control set, with within-table q-values. |
| `spatial_reliability_full_full.csv` | 28 | Split-half correlations and Spearman–Brown estimates per transform, aperture, and split. |
| `spatial_splits_full_full.csv` | 25 | Historical random-split summaries (retained observations, not freshly reproduced). |
| `bfocus_features_full_grid24.csv` | 2,904 | Participant × condition × band × channel-set features for the predictive-focus benchmark. |
| `bfocus_eval_full_grid24.csv` | 768 | Associations for every benchmark configuration and analysis type, with q-values. |
| `bfocus_splits_full_grid24.csv` | 4 | Historical split summaries for the benchmark. |

Key columns: `r_lps` / `p_lps` / `q_lps` are the partial rank correlation with
fluid reasoning and its p and q values; `r_wst` is the same for vocabulary;
`r_lps_plus_wst` adds vocabulary as a control. `controls` lists the nuisance
columns used, separated by semicolons.

## `cross_validation/` — new predictive calibration (v0.3.0)

Written by `scripts/cross_validate_readouts.py`. Regenerating them should
reproduce these files byte for byte; `tests/test_analysis.py` checks this.

| File | Rows | Contents |
|---|---:|---|
| `fold_assignments.csv` | 2,331 | Which participant is held out in which fold, for leave-one-out and all 20 five-fold repeats. |
| `predictions.csv` | 27,972 | Every held-participant prediction with its matched baseline prediction. |
| `run_metrics.csv` | 252 | Pearson r, Spearman r, R², baseline R², ΔR², and RMSE per whole-cohort run. |
| `summary.csv` | 24 | Median and 5th–95th percentile of each metric per target, score, model, and protocol. |
| `manifest.json` | – | Seed, repeats, and the exact configuration used. |

## `context/phenotypes/` — retained 291-field phenotype sweep

Copied unchanged from the identity-coordinate experiment. These are evidence
extractions, not fresh reruns; the underlying predictions were not recovered.

| File | Contents |
|---|---|
| `domain_label_reconstruction.csv` | 1,164 rows: leave-one-out and in-sample correlation for every field at K = 3, 6, 10, 20. |
| `domain_recovery_summary.csv` | 20 rows: pooled R² per domain and K. |
| `evidence.json` | The four HbA1c rows with their original scoring conventions. |
| `catalog_coverage.json` | Coverage receipt written by `scripts/build_phenotype_catalog.py`. |
| `SUMMARY.md` | Original experiment summary. |

## `verification/`

`verification.json` is written by `scripts/verify.py` and records the number of
values compared, the maximum absolute difference, the headline result, and the
bootstrap interval.
