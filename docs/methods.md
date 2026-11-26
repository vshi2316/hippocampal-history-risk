# Analysis definitions

## Clinical estimands

The exposure is annual prior hippocampal log-volume loss across MRI1 and MRI2,
with positive values denoting loss. Clinical follow-up starts at MRI2 after the
specified normal-cognition eligibility assessment. Current volume is measured
at MRI2. Model-specific cognition is either the clinical rating or the continuous
cognitive measure specified in the corresponding script. Each model uses its corresponding
complete-case population.

Mayo contributes 925 participants and 177 first-impairment events to the
clinical-rating model, versus 880 participants and 168 events to the
continuous-cognition three-year model. OASIS-3 contributes 306 participants
and 59 events to the primary first-impairment analysis. Its 312-person,
52-event secondary dementia population retains six people excluded from the
primary first-impairment population. Event components retain each cohort's
released diagnostic definitions.

Reassessment requires a later observed normal cognitive assessment. Confirmed
impairment requires consecutive impaired classifications. Reversion and absence
of another classified visit are recorded separately. Confirmation analyses
retain their available follow-up and exact time origins.

Three-year risks are Cox-model marginal standardizations over the analysis
sample's covariates at the fixed cohort-specific prior-loss quartiles. Bootstrap
sampling uses participants. Script 26 retains 300 successful standardization
resamples. Incremental performance uses five repetitions of stratified five-fold
allocation and averages held-out predictions. Its 300 paired participant
resamples hold those predictions fixed. Fixed-capacity
screening uses 1000 paired resamples of those predictions, reranking within
each resample and accounting for censoring.

## Longitudinal depth

Aligned cognitive models reuse the final clinical risk sets and actual repeated
scores. Mayo global cognition and memory retain their z-score scales. OASIS-3
uses administered MMSE points. The pre-event window ends at least one day
before first impairment. Future structural models require at least two later
scans and exclude MRI2 from the outcome slope. Protocol restrictions retain
Mayo manufacturer and OASIS physical-device/series distinctions. The visit-gap
audit uses the latest observed normal classification before first impairment.

## Biological context

ADNI and HABS-HD are independent amyloid-positive populations without dementia.
Their retained samples are 94 participants/277 scans and 57/127. The HABS-HD
stable subset contains 30 participants/60 scans. MRI/PET matching is within
180 days; clinical anchoring uses 90 days and the released date precision.
Whole-hippocampal models retain their original cohort-specific covariance
specification. Tau scales are standardized within cohort, preserving tracer
differences. ADNI influence checks exclude the predefined participant/scans.

The cortical analysis covers 68 bilateral Desikan-Killiany parcels. All
coefficients and BH-adjusted results remain available, including the HABS-HD
regional null findings. Band comparisons use paired within-scan thickness
contrasts and Holm adjustment within the declared contrast family. BigBrain
is a single postmortem anatomical reference. The 2000-spin tests retain
hemisphere-specific mappings and the original seed. Amyloid-negative age
maps use the approximately 30% held-out partition of eligible cognitively
normal amyloid-negative controls. Within each cohort, identifiers were sorted
before sampling 70% into the separate reference-training partition with seed
20261022 plus the cohort index. The held-out participants supply the age-by-time
slope-interaction coefficients, with final counts determined by regional
measurement availability.

## Measurement and matched controls

The FreeSurfer/ASHS comparison uses identical ADNI scans and paired method
contrasts, including participant bootstrap uncertainty. Same-day OASIS repeats
retain 225 people and 295 anchor pairs; 70 participants have alternate
measurements at both anchors. Raw and scan-specific ICV-normalized log-volume
differences are reported separately. The pre-MRI2 participant ICV reference is
the median of daily medians over the eligible MRI1-to-MRI2 interval.
APOE and ventricular comparisons use matched complete
cases. Continuous-interval models use the original declared interval scale.

## Tau-associated cognitive change

The separate cognitive samples contain 182 ADNI and 158 HABS-HD participants.
ADNI uses annual MMSE point change, and HABS-HD uses annual delayed recall
z-score change from SEVLT_DR_ZSCORE. The exposure is temporal tau minus the
mean of posterior cingulate and lateral parietal tau. Contrast, three-region
mean tau and amyloid Centiloids are standardized within each complete cognitive
sample. The linear model adjusts for mean tau and amyloid burden, with age,
sex and baseline cognition. HC3 standard errors accompany the estimates.
Baseline cognition is matched within 365 days of tau PET. The follow-up visit
is selected within 1 to 5 years, nearest three years. Actual elapsed time is
the annualization denominator.

## Participant characteristics

Script 36 calculates descriptive statistics from the final clinical and tau
analysis populations. Clinical summaries retain model-specific cognition,
scan intervals and APOE availability. Tau summaries retain baseline
hippocampal measurements, tracer burden and structural follow-up duration.
