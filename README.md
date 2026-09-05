# An EEG correlate of fluid intelligence

**Ximon ([zw5](https://github.com/zw5))**

A method for measuring the concentration of theta-band EEG dynamics, with a fluid-reasoning correlate in the LEMON dataset.

In **111 participants**, eyes-open posterior theta concentration correlates with fluid reasoning at **partial Spearman r = 0.408** after adjustment for age cohort, age-bin midpoint, theta power, and static phase locking. The association remains **r = 0.408** when vocabulary is included as an additional control.

![EEG concentration, fluid reasoning, and vocabulary in 111 participants](figures/correlate.png)

| Measure | Result |
|---|---:|
| Fluid reasoning · LPS | **r = 0.408**, p = 1.31 × 10⁻⁵ |
| Vocabulary · WST | r = 0.001 |
| Participant-bootstrap 95% interval for the LPS association | [0.186, 0.549] |
| Within-recording consistency · odd/even windows | r = 0.954 |
| Participants | 76 younger adults · 35 older adults |

Correlations with LPS and WST use the same four controls. The bootstrap interval is conditional on the historically selected score. Within-recording consistency describes repeated measurements from one recording.

[**Read the research note →**](paper.md)

## How the method works

The method filters the EEG to **4–8 Hz** and constructs a complex lagged cross-field operator from its analytic phase. The operator describes the relationship between the sensor field at one time and the field a fraction of a theta cycle later.

The score measures how concentrated the operator's singular-value energy is, using normalized spectral entropy. Every singular value contributes. Lags follow each participant's measured theta carrier at quarter-, third-, and half-cycle offsets; the score averages across these lags and recording windows.

The main result uses the untransformed sensor field. The release also includes the earlier fixed-lag Beam focus benchmark, condition and frequency comparisons, null-normalized scores, and spatial smoothing/high-boost ablations.

## Reproduce the statistics

```sh
python3.11 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python scripts/verify.py
```

The verifier recomputes **880 association rows**, including their within-table multiple-comparison corrections, and checks **28 reliability rows**. All **7,808 compared association values match the retained results exactly**. It also regenerates the figure and the participant-bootstrap interval.

| Material | Contents |
|---|---|
| [Research note](paper.md) | Equations, sample, controls, results, and interpretation |
| [Participant features and results](results/retained/) | Complete tables for the two included experiments |
| [Verification](results/verification/verification.json) | Numerical agreement and bootstrap results |
| [Estimator source](source_original/) | Original implementation and supporting source files |
| [Provenance](PROVENANCE.md) | Input history, reproduction scope, and missing raw-data records |

This release reproduces the statistics from retained participant features. The original EEG export must be recovered to repeat EEG-to-feature extraction; the archived estimator files preserve their historical paths and imports.

## What the result establishes

Theta lagged-field concentration is associated with measured fluid reasoning in this LEMON cohort after the specified adjustments. The near-zero vocabulary association provides a useful comparison, although a formal difference between the two associations was not tested.

This is an exploratory result from a research program that examined several representations. Within-cohort splits measure stability; independent-cohort validation remains open. The q-values correct the named result tables, and the bootstrap interval conditions on the chosen score. The research note includes the null-normalized variants that did not retain the primary association.

## Dataset and citation

The human data were collected by **Babayan and colleagues** for the Leipzig Mind-Brain-Body dataset (LEMON). LPS-2 subtest 3 measures fluid reasoning; WST measures vocabulary and crystallized intelligence.

Babayan A, et al. *A mind-brain-body dataset of MRI, EEG, cognition, emotion, and peripheral physiology in young and old adults.* Scientific Data **6**, 180308 (2019). [doi:10.1038/sdata.2018.308](https://doi.org/10.1038/sdata.2018.308).

The [official LEMON page](https://fcon_1000.projects.nitrc.org/indi/retro/MPI_LEMON/MPI_LEMON.html) provides data access and a PDDL data-sharing link. EEG and behavioral files are also available through the [GWDG distribution](https://ftp.gwdg.de/pub/misc/MPI-Leipzig_Mind-Brain-Body-LEMON/). This repository distributes derived participant tables. Please cite the dataset descriptor alongside this work.
