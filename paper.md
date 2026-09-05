# Flint: An EEG correlate of fluid intelligence

*Theta lagged-field concentration and fluid reasoning in LEMON*

Ximon (zw5)

Research note · 5 September 2026

## Abstract

We describe an EEG method based on the spectral concentration of theta-band lagged cross-field geometry and report its association with fluid reasoning in 111 LEMON participants. In the untransformed broad posterior sensor field, concentration correlates with LPS fluid reasoning at partial Spearman r = 0.408 (p = 1.31 × 10⁻⁵), adjusting for cohort, age-bin midpoint, theta power, and static phase locking. The corresponding vocabulary correlation is 0.001, and adding vocabulary as a control preserves the reasoning association. The score uses the complete singular-value energy spectrum and averages across participant-carrier-aligned phase lags and recording windows. Within-recording odd/even and first/second split correlations are 0.954 and 0.948. Spatial smoothing attenuates the association, while high-boost ablations reach approximately 0.44 in individual lag rows. The release provides the participant feature tables, estimator source, and a statistical verifier that reproduces all 7,808 compared values across 880 association rows exactly. The result is an exploratory correlate in the original cohort; independent-cohort validation and a fresh EEG-to-feature reproduction remain open.

## Data and outcome

The analysis uses retained LEMON features from 111 participants with eyes-open EEG, fluid-reasoning and vocabulary scores, and age information: 76 in the younger cohort and 35 in the older cohort. Participants, rather than windows, are the statistical units. LPS denotes the exported first numerical LPS column; the original export code reads column index 1. WST denotes column index 3 of its behavioral table. The dataset descriptor identifies the administered reasoning task as LPS-2 subtest 3 and the vocabulary task as WST [1]. Age is the midpoint of the distributed age bin.

The earlier benchmark additionally includes 131 eyes-closed participants. Comparisons between eyes-open and eyes-closed results therefore have unequal denominators; they are not a paired test of the condition interaction. The historical export required usable eyes-closed data and exported eyes-open data when available. A complete original participant-exclusion ledger and the exact channel manifest are not present in this checkout.

## Operator and readout

The EEG pipeline resamples the preprocessed EEG to 125 Hz when necessary, using a common ordered channel set. For the concentration lineage, the signal is filtered to 4–8 Hz with a fourth-order Butterworth filter applied forward and backward. Its complex analytic phase is obtained by the Hilbert transform and amplitude normalization with the original numerical stabilizer:

```math
Z_j(t)=\frac{\operatorname{Hilbert}(X_{\theta,j})(t)}{|\operatorname{Hilbert}(X_{\theta,j})(t)|+10^{-9}}.
```

The broad posterior aperture is defined in source as channels not passing its frontal-label predicate. This is a sensor-label aperture, not anatomical source localization. The raw lane applies no spatial transform.

Within each window, the lagged cross-field matrix is

```math
C_{\tau}=\frac{1}{T-\tau}\sum_t Z(t+\tau)Z(t)^*,\qquad
p_i=\frac{\sigma_i(C_{\tau})^2}{\sum_j\sigma_j(C_{\tau})^2}.
```

The source does not subtract a temporal mean in this matrix calculation. The concentration readout is

```math
c_{\tau}=1-\frac{-\sum_i p_i\log p_i}{\log d},
```

where d is the number of singular values. Every singular value enters the entropy.

The carrier frequency is the peak of the aperture-averaged Welch spectrum in 4–8 Hz. Lags are the nearest sample to one quarter, one third, and one half of that carrier cycle. Scores average the three lag concentrations within each eight-second window, then average across the selected windows. The retained run used at most 24 windows; the source uses a four-second candidate step and evenly spaced selection across candidates when the cap is exceeded. The source phase field is computed over the full exported condition before window extraction. Consequently, within-recording splits should not be interpreted as independent recording sessions.

The earlier benchmark is a distinct functional of the same lagged spectrum: leading energy fraction divided by entropy effective rank, or z(leading energy fraction) minus z(effective rank). It uses 40, 80, 120, and 160 ms lags. Standardization for the z-difference is performed across participants within each condition/band/aperture configuration before its split analysis.

## Statistical analysis

Rank-transform the score, outcome, and each control. Regress the ranked score and ranked outcome separately on an intercept and the ranked controls. Correlate their residuals. The retained two-sided p-value uses a t approximation with n minus number-of-controls minus two degrees of freedom. Controls are cohort code, age-bin midpoint, theta power, and static phase-locking value measured in the corresponding representation. An additional analysis includes vocabulary; another adds theta carrier frequency. This procedure estimates an adjusted association between participants.

The independent verifier recomputes the entire 112-row spatial table and 768-row condition/band benchmark, including their Benjamini–Hochberg corrections. It also checks 28 reliability rows. A new 3,000-replicate participant bootstrap reranks and refits controls within every bootstrap sample, with seed 20260905. Its percentile 95% interval for the raw headline score is [0.186, 0.549], conditional on historical feature choice.

## Results

| Retained analysis | n | Partial r with LPS | p | q within its table | Partial r with WST |
|---|---:|---:|---:|---:|---:|
| Raw broad-posterior carrier-aligned concentration | 111 | 0.408 | 1.31e-5 | 4.62e-5 | 0.001 |
| Raw lateral-posterior carrier-aligned concentration | 111 | 0.398 | 2.20e-5 | 5.42e-5 | -0.002 |
| Earlier posterior theta focus z-difference, eyes open | 111 | 0.369 | 9.02e-5 | 0.00770 | -0.002 |
| Earlier posterior theta focus z-difference, younger cohort | 76 | 0.371 | 0.00124 | 0.0306 | 0.057 |

Controlling additionally for vocabulary leaves the headline association at 0.408. Adding theta peak frequency gives 0.407. These observations support differentiation between the measured reasoning and vocabulary outcomes within this sample; a near-zero vocabulary correlation alone is not a formal test that the two associations differ.

The raw broad-posterior concentration score has odd/even window split r = 0.954 and first/second split r = 0.948. The corresponding Spearman–Brown values are 0.976 and 0.973. They describe within-recording consistency, with potentially overlapping windows and shared preprocessing, rather than independent-session test–retest reliability. Historical participant splits give median held-half r = 0.376 for the raw broad-posterior score. These repeatedly reuse the original cohort and follow earlier feature exploration.

Spatial ablations preserve the raw lane. Gaussian smoothing reduces the broad-posterior association to 0.223 and 0.196 for its two widths. High-boost variants give aggregate correlations of 0.382–0.424 for the broad posterior aperture; individual lag rows reach 0.438. These are retained exploratory comparisons, not independent confirmations of a superior estimator. The raw result is the headline because it directly characterizes the untransformed field.

The earlier benchmark's three posterior-theta null-normalized focus variants yield r = 0.115 (synchronous time permutation), -0.173 (independent channel circular shifts), and -0.069 (timewise channel permutation). None is significant at the nominal 0.05 level. These variants divide the observed-minus-null difference by a null standard deviation and thus measure a different functional from the observed operator score. Their negative result limits claims of null-normalized temporal specificity; it does not change the observed raw association. A zero-lag full-operator comparison is not independently reproduced in this release.

## Interpretation and scope

The empirical claim is that theta lagged-field spectral concentration covaries with measured fluid reasoning in this cohort after the specified adjustments. It is a correlate of a cognitive outcome. It does not establish causal enhancement, moment-to-moment intelligence measurement, or an independent-cohort prediction result.

The retained implementation developed through multiple representation and estimator experiments. Terms such as “locked” in historical filenames and summaries refer to settings fixed for particular runs; they do not establish preregistration before inspection of LEMON outcomes. Within-table q-values do not account for every earlier research decision. The release preserves the positive observations together with their actual selection history and comparison results.

The next reproduction step is recovery of the original channel/participant export and a fresh EEG-to-feature calculation against these retained rows. A separate cohort evaluated with the frozen raw operator would address external generalization. Both steps can use the same operator and score definition.

## Availability and attribution

The accompanying repository contains all retained participant features for the two tables, their original evaluations and split summaries, unchanged source files, and the statistical verifier. The source files retain their historical absolute paths and supporting project imports; they are provenance records, not a portable raw-data reproduction command. Only the documented verifier is the tested reproduction entry point in this release.

The original human data and their acquisition are credited to Babayan and the LEMON investigators [1].

[1] Babayan A, et al. A mind-brain-body dataset of MRI, EEG, cognition, emotion, and peripheral physiology in young and old adults. Scientific Data 6, 180308 (2019). https://doi.org/10.1038/sdata.2018.308. [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6371893/).
