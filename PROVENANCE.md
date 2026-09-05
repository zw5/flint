# Release provenance

The release was assembled from the live `/Users/cio/Documents/interp` workspace on 2026-09-05. That workspace root has no Git repository, so no original root commit is claimed.

## Retained inputs

- `experiments/eeg_coherence/lemon_bfocus_benchmark/`: three `*full_grid24.csv` files, containing participant/configuration features, evaluation rows, and split summary rows.
- `experiments/eeg_coherence/lemon_spatial_filter_ceiling/`: four `*full_full.csv` files, containing participant features, evaluations, reliability, and split summaries.
- `scripts/`: the seven exact source files copied into `source_original/`, including the export script, both evaluators, and operator helpers. These retain their historical paths and are not presented as a portable execution package.

No participant table was re-estimated, sharpened, normalized, or selected to improve the publication result. The raw and transformed lanes are all retained. The `full_grid24` benchmark contains both eye conditions, four bands, three channel sets, eight score readouts, and four control/cohort analysis types. Its q-values use all 768 association rows. The spatial table contains seven transforms, two apertures, four readouts, and two control sets, for 112 association rows.

## Fresh work

Command from the workspace root:

```sh
.venv/bin/python experiments/fluid_intelligence_publication_20260905/publication/scripts/verify.py
```

This recomputes rank-residual associations, t-approximation p-values, within-table q-values, and reliability from the copied participant tables. It also generates a conditional bootstrap interval and the figure. The receipt records 7,808 exact association-value matches and 28 checked reliability rows. Tolerances are relative 1e-9 and absolute 1e-10. Statistical values, not document hashes, are the verification evidence.

The historical split summaries are retained observations. Their full random split assignments and original commands are not retained here, so this release does not claim a fresh reproduction of those split summaries. Source defaults alone are not evidence of the exact original invocation.

## Access and selection history

The EEG feature extraction functions do not fit to LPS or WST. Cognitive labels were available to the historical association analyses and informed the wider research program. The release authors inspected all included evaluation tables. The raw score was chosen as the main reported result while retaining the full tables. The descriptive bootstrap interval does not remove historical selection effects.

## Outstanding raw-data provenance

The original `rust_export` directory, its `channels.json` and `manifest.csv`, the source EEG files, the complete participant exclusion ledger, and exact original shell invocations were not found in this checkout. The export source specifies 125 Hz resampling when required and intersection-based channel selection, but source intent alone cannot certify the literal missing export. No raw-EEG reproduction is claimed.

The official LEMON source is https://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON/ . The dataset descriptor is https://doi.org/10.1038/sdata.2018.308 . The official INDI dataset page links the Public Domain Dedication and License (PDDL): https://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON/MPI_LEMON.html . Preserve the LEMON attribution when reusing these derived tables.
