# Flint: An EEG correlate of fluid intelligence

*Theta lagged-field concentration, fluid reasoning, and broader phenotype readouts in LEMON*

Simon Velez

Research note · 5 September 2026

## Abstract

We describe an EEG method based on the spectral concentration of theta-band lagged cross-field geometry and report its association with fluid reasoning in 111 LEMON participants. In the untransformed broad posterior sensor field, concentration correlates with LPS fluid reasoning at partial Spearman r = 0.408 (p = 1.31 × 10⁻⁵), adjusting for cohort, age-bin midpoint, theta power, and static phase locking. The corresponding vocabulary correlation is 0.001, and adding vocabulary as a control preserves the reasoning association. The score uses the complete singular-value energy spectrum and averages across participant-carrier-aligned phase lags and recording windows. Within-recording odd/even and first/second split correlations are 0.954 and 0.948. Spatial smoothing attenuates the association, while high-boost ablations reach approximately 0.44 in individual lag rows. The release provides the participant feature tables, estimator source, and a statistical verifier that reproduces all 7,808 compared values across 880 association rows exactly. A related five-band spatial-identity analysis contributes 1,164 retained evaluations covering 291 physiological and behavioral fields. These include positive individual prediction correlations for HbA1c, blood pressure, body measurements, reasoning, and selected questionnaire scales, with heterogeneous performance across entire domains. The results are exploratory findings in the original cohort; independent-cohort validation and fresh EEG-to-feature reproduction remain open.

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

## A broader phenotype map

The same research program examined whether a repeatable EEG identity representation could support readouts across several types of measurement. The retained domain-recovery experiment includes 291 numeric phenotype fields: 35 blood-chemistry fields, 11 blood-pressure and pulse fields, four anthropometric fields, 101 cognitive fields, and 140 personality/affect fields. All were evaluated at K = 3, 6, 10, and 20 identity coordinates, producing 1,164 label-level rows and 20 domain summaries. The complete [phenotype atlas](PHENOTYPE_ATLAS.md) renders every result in source-name order.

This analysis has a shared EEG representation and a separate phenotype map for each domain. “Shared” refers to the representation and operator family, not to one fitted scalar that predicts all outcomes. The coordinate count is the same within a given K comparison. The reliability axes and domain map are refitted within each held-participant fold. The following HbA1c section describes the construction in detail.

At K = 6, illustrative retained Pearson correlations are +0.494 for waist circumference, +0.472 for the second left systolic blood-pressure measurement, +0.443 for LPS_1, +0.434 for the CKD-EPI laboratory field, +0.345 for HbA1c percentage, +0.319 for weight, +0.318 for ASAT, +0.308 for hip circumference, +0.303 for UPPS lack of perseverance, and +0.302 for LOT-R optimism. These rows illustrate the breadth of the target field; they are not a list of independently validated biomarkers or the result of corrected per-label hypothesis tests.

| Domain | Fields | Pooled LOOCV R², K = 3 | Pooled LOOCV R², K = 6 |
|---|---:|---:|---:|
| Blood chemistry | 35 | +0.0098 | -0.0105 |
| Blood pressure | 11 | +0.1263 | +0.1008 |
| Body measurements | 4 | +0.1192 | +0.0871 |
| Cognition | 101 | -0.0074 | -0.0535 |
| Personality and affect | 140 | -0.0089 | -0.0500 |

The aggregate results distinguish a domain's shared recoverability from its strongest individual outcomes. At these settings the blood-pressure and body-measurement domains show positive pooled standardized R². The wider cognition and personality/affect fields do not, despite individual positive associations. Reported negative prediction correlations are also important: at K = 6, TSH is -0.454 and the NYC-Q item 28 is -0.368. These indicate an inverted relationship between the fitted prediction and the outcome. They cannot be reinterpreted as successful predictions simply by taking absolute values, or as the direction of a causal EEG–trait relationship.

The recorded labels have different meanings and dependencies. Several blood-pressure rows repeat measurements across sides or occasions; HbA1c has alternate units; cognitive and personality rows can be test components, individual items, or composite scales. The atlas preserves the original names and signs and supplies observed-label counts. Missing target handling follows the common source procedure described below. The sweep does not report multiplicity-adjusted p-values for these 1,164 correlations, and its label inclusion and K exploration were not nested into an external cohort.

### Why the breadth is consequential

The research proposition is that EEG geometry can be used as a common starting point for studying several aspects of a person's functioning. A cognitive task, a blood sample, a cuff measurement, and a questionnaire observe different processes with different instruments. Finding structured relationships to all of them suggests a route toward studying their shared and distinct physiological organization.

The value of the operator approach is that it supplies coordinates for this question. Temporal lags expose relationships across time; spectra describe how those relationships concentrate; spatial envelopes describe where a repeatable mode is expressed; and cross-covariance maps connect those readouts to other measurements. That makes it possible to ask which relationships share a carrier and which require distinct readouts, rather than treating every target as an unrelated prediction problem.

The pattern also poses testable alternatives. Common age structure could account for part of several outcomes at once. Stable metabolic or vascular differences could shape aspects of neural coordination. Other associations could depend on behavior, arousal, or measurement conventions. These explanations have different predictions under covariate adjustment, cross-cohort transfer, and repeated observations of the same person. The published map provides a concrete empirical starting point for distinguishing them.

A reproducible account of those relationships would be a substantial scientific advance: it would connect electrical dynamics at the scalp to a broader description of cognition and physiology. The present evidence motivates that program through actual retained results, with its strongest and weakest domains visible together.

## Physiological extension: HbA1c

### Motivation

A fluid-reasoning correlate asks whether neural organization is informative about an ability measured through behavior. HbA1c adds a different kind of outcome: a laboratory measure of glycated hemoglobin reflecting average glucose exposure over approximately three months [2]. Its association with a short resting EEG observation raises a question about how persistent physiology is expressed in fast neural dynamics.

The proposed explanation is shared physiological organization. Neural activity occurs within metabolic and vascular conditions that can persist beyond the observation window. If those conditions systematically shape coordination, a repeatable EEG coordinate can contain information about a laboratory measurement taken through another modality. This interpretation does not require the recorded brain activity to encode the assay value explicitly. It requires a reproducible statistical relationship between the two observations, and further evidence to distinguish its possible causes.

There is relevant empirical precedent for a connection between glycemic regulation and electrophysiology. Cooray and colleagues studied 28 people with type 2 diabetes and 21 controls; among the patients, a group receiving intensified glycemic control showed changes in cognition, resting EEG alpha activity, and connectivity after two months [3]. That study supports investigating the connection, but it neither validates this Beam estimator nor identifies the cause of the LEMON association.

### Operator lineage and prediction procedure

The retained HbA1c result comes from `lemon_identity_yreg_domain_recovery.py`. Its input is a broader identity representation built from delta, theta, alpha, beta, and gamma bands. For each band and window, the source constructs the complex lagged cross-field operator and takes the absolute spatial envelope of its leading left singular vector. These envelopes are normalized and averaged within recording splits, then assembled across bands. This is a representation derived from the lagged operator; it is distinct from the full-spectrum theta concentration scalar used in the primary analysis.

The reliability basis is obtained from the cross-covariance of separately standardized odd and even identity arrays. Each retained axis combines the corresponding left and right singular vectors. The domain analysis examines K = 3, 6, 10, and 20 coordinates. In each leave-one-out fold, it refits the reliability basis on the remaining participants, projects the held participant using the training transformation, and standardizes the projected coordinates using training statistics.

A second cross-covariance maps the resulting EEG coordinates to the standardized blood-chemistry field. Its singular-value decomposition is retained up to rank min(K, 10). The held participant's EEG coordinates are contracted through that map to produce standardized label predictions. This is the original cross-covariance calculation; it should not be described as ordinary least-squares regression.

### Retained observations

| K | Reported HbA1c leave-one-out Pearson r | In-sample Pearson r |
|---|---:|---:|
| 3 | 0.369105 | 0.418033 |
| 6 | 0.345048 | 0.443139 |
| 10 | 0.319307 | 0.472136 |
| 20 | 0.285914 | 0.533383 |

The source analysis contains 111 participants and reports 109 finite original HbA1c observations. Its standardization helper replaces missing standardized targets with zero. Consequently, the reported correlation includes the two mean-imputed target entries; 109 is the observed-label count, not the literal correlation denominator. The leave-one-out predictions are in training-standardized units, while the final target table is standardized globally. These conventions are preserved in the archived source and should be reproduced explicitly before reporting a complete-case or clinically calibrated estimate.

The analysis does not residualize age, cohort, sex, body composition, or EEG power before estimating this HbA1c relationship. Its Pearson prediction correlation is therefore different from the nuisance-adjusted partial rank correlation reported for fluid reasoning. The percentage and mmol/mol HbA1c columns describe the same assay in different units. They have nearly identical correlations and supply one biological finding.

The dimension sweep is informative: higher K improves the in-sample fit while reducing the held-out correlation across these four settings. The strongest retained association is present at K = 3. Because the choice of K is examined across the original cohort, this is an exploratory comparison rather than a separately validated optimum.

The complete domain tables are included to retain the broader context. Across all 35 blood-chemistry labels, pooled standardized LOOCV R² is 0.0098 at K = 3 and -0.0105 at K = 6. The positive HbA1c row therefore does not imply accurate reconstruction of the full blood panel. The same experiment reports positive domain-level values for blood pressure and body size at these settings, which motivates further work on their shared structure.

### Scientific significance

The consequential possibility is a bridge between cognitive and physiological measurement. The primary finding relates an explicit measure of theta coordination to reasoning. The physiological extension asks whether repeatable spatial organization in the same operator family is informative about a longer-term metabolic measure. Together, they motivate studying how a person's neural dynamics reflect several aspects of their functioning, with different readouts preserving different parts of the operator.

The distinction between readouts is useful. It lets the research ask whether reasoning and metabolic state occupy shared coordinates, separable coordinates, or mixtures of both. A single positive scalar cannot resolve that question. The original operator geometry provides a way to formulate it: compare spectra, spatial envelopes, reliability directions, and cross-modal mappings while retaining the underlying field.

If replicated with frozen definitions, observed-label scoring, relevant covariates, and repeated EEG/laboratory observations, this approach could help characterize how neural coordination relates to persistent physiological conditions and changes within a person. That would extend the value of EEG beyond the immediate electrical trace toward an interpretable description of the organism's state. The current cohort establishes a concrete association worth pursuing; the mechanism, external generalization, and within-person relationship remain separate empirical questions.

The HbA1c tables are retained experimental outputs, not newly recomputed predictions. The original identity tensor, joined phenotype matrix, and per-participant prediction records are absent from this release workspace. Their recovery is needed for a fresh numerical reproduction. HbA1c also has determinants beyond glucose exposure, including red-cell turnover [2], so a causal account must be more specific than a generic claim that EEG reads metabolism.

## Interpretation and scope

The empirical claim is that theta lagged-field spectral concentration covaries with measured fluid reasoning in this cohort after the specified adjustments. It is a correlate of a cognitive outcome. It does not establish causal enhancement, moment-to-moment intelligence measurement, or an independent-cohort prediction result.

The retained implementation developed through multiple representation and estimator experiments. Terms such as “locked” in historical filenames and summaries refer to settings fixed for particular runs; they do not establish preregistration before inspection of LEMON outcomes. Within-table q-values do not account for every earlier research decision. The release preserves the positive observations together with their actual selection history and comparison results.

The next reproduction step is recovery of the original channel/participant export and a fresh EEG-to-feature calculation against these retained rows. A separate cohort evaluated with the frozen raw operator would address external generalization. Both steps can use the same operator and score definition.

## Availability and attribution

The accompanying repository contains all retained participant features for the two primary tables, their original evaluations and split summaries, unchanged source files, and the statistical verifier. The physiological supplement adds the full retained domain-recovery tables, their summary, and the source defining the identity and cross-covariance readouts. The source files retain their historical absolute paths and supporting project imports; they are provenance records, not a portable raw-data reproduction command. Only the documented verifier is the tested reproduction entry point in this release.

The original human data and their acquisition are credited to Babayan and the LEMON investigators [1].

[1] Babayan A, et al. A mind-brain-body dataset of MRI, EEG, cognition, emotion, and peripheral physiology in young and old adults. Scientific Data 6, 180308 (2019). https://doi.org/10.1038/sdata.2018.308. [Open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6371893/).

[2] National Institute of Diabetes and Digestive and Kidney Diseases. [The A1C Test & Diabetes](https://www.niddk.nih.gov/health-information/diagnostic-tests/a1c-test). Accessed 5 September 2026.

[3] Cooray G, Nilsson E, Wahlin A, Laukka EJ, Brismar K, Brismar T. *Effects of intensified metabolic control on CNS function in type 2 diabetes.* Psychoneuroendocrinology (2011). [doi:10.1016/j.psyneuen.2010.06.009](https://doi.org/10.1016/j.psyneuen.2010.06.009).
