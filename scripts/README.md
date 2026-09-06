# Runnable analysis scripts

Every script reads the retained tables under `results/` and writes reports,
tables, or figures back into the repository. None of them touches raw EEG; the
EEG-to-feature stage is archived separately in `source_original/`.

Run them from the repository root with the virtual environment described in the
top-level README, or use the `Makefile` targets.

| Script | Reads | Writes | Make target |
|---|---|---|---|
| `verify.py` | `results/retained/*.csv` | `results/verification/verification.json`, `figures/correlate.{png,svg}` | `make verify` |
| `cross_validate_readouts.py` | `results/retained/*.csv` | `results/cross_validation/*.csv`, `manifest.json` | `make cv` |
| `write_cv_report.py` | `results/cross_validation/summary.csv`, `templates/cross_validation.md` | `CROSS_VALIDATION.md` | `make report` |
| `build_phenotype_catalog.py` | `results/context/phenotypes/*.csv`, `templates/phenotype_atlas.md` | `PHENOTYPE_ATLAS.md`, `catalog_coverage.json` | `make atlas` |
| `plot_results.py` | `results/**` | `figures/*` except `correlate.*`, `figure_manifest.json` | `make figures` |
| `data_io.py` | – | – | shared CSV helpers imported by the others |

`templates/` holds the Markdown prose for the two generated reports. The
scripts substitute tables into those templates, so wording changes go in the
template and numbers always come from the CSV files.

Order matters only in one place: `write_cv_report.py` and the cross-validation
figures in `plot_results.py` read the outputs of `cross_validate_readouts.py`.
Every script is deterministic; rerunning it should leave `git status` clean.
