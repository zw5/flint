# Theta lagged-field concentration correlates with fluid reasoning in LEMON

**Ximon ([zw5](https://github.com/zw5)) · Research release · 5 September 2026**

In **111 LEMON participants**, a Beam readout of eyes-open theta lagged-field concentration is associated with fluid reasoning: **partial Spearman r = 0.408, p = 1.31 × 10⁻⁵** after adjustment for age cohort, age-bin midpoint, theta power, and static phase-locking value. The vocabulary association is **r = 0.001**. Including vocabulary as an additional control leaves the fluid-reasoning result at **r = 0.408**.

![Fluid reasoning, vocabulary, and retained spatial comparisons](figures/correlate.png)

The score measures spectral concentration of the complex lagged cross-field operator. It retains the full singular-value spectrum when computing entropy. Its lags are quarter-, third-, and half-cycles of each participant's measured theta carrier. The headline uses the original sensor field without a spatial smoothing or high-boost transform.

This release includes the participant feature tables, all **880 association rows**, original estimator source, a research note, and an independent statistical verifier. **7,808 numerical values reproduce exactly from the retained tables**, and all 28 within-recording reliability rows pass numerical checks. This is a statistical reproduction from retained EEG features. The original EEG export is absent from the release workspace; the EEG-to-feature computation was not rerun for this release.

## Read and reproduce

- [Research note](paper.md): operator, cohort, controls, findings, and interpretation.
- [Verification receipt](results/verification/verification.json): exact comparison counts and a new descriptive participant-bootstrap interval.
- [Provenance](PROVENANCE.md): what is retained, what was recomputed, and what remains to recover.
- [Original source](source_original/): unchanged source files for the selected lineage and its direct supporting modules.

```sh
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/verify.py
```

The verifier checks every retained association and its within-table Benjamini–Hochberg q-value, checks the reliability rows, and regenerates the figure. It raises an error on a mismatch. The figure plots residual ranks; the points are participants, and the lines are descriptive fits rather than predicted IQ scores. Teal denotes the younger cohort and rust denotes the older cohort.

## Scope

This is an exploratory association in the original cohort. Within-cohort splits characterize stability; they do not constitute independent-cohort replication or fully nested discovery validation. The historical searches evaluated multiple representations. The reported q-values correct their named tables, not the entire research history. The new bootstrap interval, **[0.186, 0.549]**, is conditional on the selected raw score.

The earlier focus benchmark provides a separate condition/frequency comparison: eyes-open posterior theta focus has r = 0.369, with r = 0.371 in the younger participants alone. Its three null-normalized variants do not reproduce the primary association. These results accompany the positive finding and constrain its mechanistic interpretation.

## Data and credit

LEMON was collected by Babayan and colleagues. Its dataset descriptor identifies LPS-2 subtest 3 as a fluid-intelligence measure and WST as a vocabulary/crystallized-intelligence measure: [Babayan et al., Scientific Data 6, 180308 (2019)](https://doi.org/10.1038/sdata.2018.308). Data access and the PDDL data-sharing link are provided on the [official LEMON page](https://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON/MPI_LEMON.html); the [GWDG distribution](https://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON/) provides EEG and behavioral data. These are derived participant tables, not a redistribution of raw EEG. Retain attribution to the dataset authors when reusing them.

Research authorship is Ximon's. Codex assisted with source inspection, statistical recomputation, figure generation, and release preparation.
