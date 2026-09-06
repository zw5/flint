# Archived estimator source

These twelve files are the original EEG-extraction and analysis scripts, copied
unchanged from the research workspace. They are preserved as a record of how the
retained participant tables in `results/` were produced. They are **not** a
runnable package: they keep their historical paths, expect an EEG export that is
not distributed here, and depend on libraries outside `requirements.txt`. See
`PROVENANCE.md` for what is and is not reproduced.

## Export stage

| File | Purpose |
|---|---|
| `export_for_rust.py` | Loads LEMON eyes-closed and eyes-open resting EEG, selects the channel set common to all participants, and writes flat binaries plus a participant manifest. |

## Fluid-reasoning lineage

These produced the tables in `results/retained/` and the numbers the verifier
checks.

| File | Purpose |
|---|---|
| `lemon_intelligence_residual_search.py` | Searches extracted feature tables for EEG dynamics that track LPS reasoning while testing cohort, age, vocabulary, and power confounds. |
| `lemon_bfocus_benchmark.py` | Benchmarks the predictive-focus score (leading energy over effective rank) across conditions, bands, channel sets, and null families. Produces `bfocus_*_full_grid24.csv`. |
| `lemon_frequency_sliced_beam_eval.py` | Replaces fixed millisecond lags with fractions of the theta cycle and evaluates spectral concentration. |
| `lemon_beam_hypothesis_dossier.py` | Tests aperture geometry, directional propagation, and nested validity hypotheses on held-out splits. |
| `lemon_beam_reliability_carrier.py` | Estimates each participant's theta carrier, compares fixed with carrier-aligned lags, and measures split reliability. |
| `lemon_spatial_filter_ceiling.py` | Applies fixed outcome-blind spatial filters (smoothing and high-boost) and reports the concentration score under each. Produces `spatial_*_full_full.csv`. |

## Spatial-identity and phenotype lineage

These produced the tables in `results/context/phenotypes/`.

| File | Purpose |
|---|---|
| `lemon_shared_identity_sex_beam.py` | Builds label-blind five-band identity coordinates from leading spatial envelopes and tests them against sex labels. |
| `lemon_shared_identity_label_scan.py` | Scans all numeric LEMON labels against the frozen identity coordinates. |
| `lemon_identity_yreg_object.py` | Treats the label table as a multivariate response and asks what low-rank part is recoverable from identity coordinates. |
| `lemon_identity_yreg_highdim.py` | Sweeps the number of retained identity coordinates (K) using an odd/even reliability basis. |
| `lemon_identity_yreg_domain_recovery.py` | Compares recoverability across physiological, cognitive, and personality domains. Produces the domain tables. |
