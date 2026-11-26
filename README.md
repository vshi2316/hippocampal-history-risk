# Hippocampal history and clinical risk

Analysis code for **Prior hippocampal loss distinguishes three-year absolute risk of cognitive impairment beyond current volume and cognition in cognitively normal adults**.

Mayo/MCSA and OASIS-3 quantify subsequent clinical risk among adults who remain
cognitively normal at their second MRI. Independent ADNI and HABS-HD populations
provide tau imaging context. The repository contains cohort preparation and
statistical analyses, with participant summaries.

## Analysis coverage

| Cohort | Analysis populations | Scripts |
|---|---|---|
| Mayo/MCSA | 925 participants/177 first-impairment events; 880/168 with continuous cognition; 527 at normal reassessment | 18-24 and 26-31 |
| OASIS-3 | 306 participants/59 first-impairment events; 219 at reassessment; 312/52 for secondary dementia | 13-17 and 22-32 |
| ADNI | 94 participants/277 scans in the primary tau models | 01-12 and 33-35 |
| HABS-HD | 57 participants/127 scans in the primary tau models; 30/60 with stable acquisition | 01-03, 06-12 and 33-35 |

| Paper evidence | Scripts | Output under the results directory |
|---|---|---|
| Clinical events and three-year risks | 13-18, 22, 26 | `cn_progression_62/`, `absolute_risk_65/` |
| Subsequent cognition and hippocampal change | 23-24, 29 | `aligned_longitudinal_68/`, `prediagnostic_64/` |
| Confirmation and normal reassessment | 24, 26, 28 | `confirmed_impairment_67/` |
| Tau, regional anatomy and cognitive change | 01-12, 33-35 | `tau_hippocampal_models/`, `bigbrain_spin_60/`, `tau_recall_models/` |
| Measurement agreement and matched controls | 05, 30-32 | `paired_segmentation/`, `structural_history_70/`, `icv_sensitivity_71/` |
| Held-out performance and fixed-capacity screening | 27, 31 | `history_increment_66/`, `structural_history_70/` |
| Participant characteristics | 36 | `descriptive_statistics/` |

[analysis-map.tsv](docs/analysis-map.tsv) maps the paper results to scripts.
[methods.md](docs/methods.md) records the statistical definitions.
The Q3-minus-Q1 three-year risk contrast averages predictions over the
remaining observed covariates. Model risk tertiles use the complete
multivariable prediction. All reported negative results are retained.

## Data preparation

Obtain authorized extracts through [ADNI](https://adni.loni.usc.edu/data-samples/adni-data/),
[HABS-HD](https://apps.unthsc.edu/itr/), [OASIS](https://sites.wustl.edu/oasisbrains/)
and the [Mayo Clinic Study of Aging](https://www.mayo.edu/research/centers-programs/alzheimers-disease-research-center/research-activities/mayo-clinic-study-aging).
Arrange the dated files and archive members as listed in
[input-files.tsv](docs/input-files.tsv). Each cohort's access conditions apply.
Participant records remain in the authorized local environment.

## Execution

The completed analyses used R 4.4.2 on Windows. Recorded package versions and
their provenance are in [package-versions.tsv](docs/package-versions.tsv).
The installer obtains missing packages from CRAN. Use the recorded versions
for numerical reproduction.

Open R in the repository root. Inputs default to `data/`, and outputs default
to `results/`. Set `HIPPO_HISTORY_DATA` and `HIPPO_HISTORY_RESULTS` before
sourcing the entry point to use other locations.

```r
source("install_packages.R", encoding = "UTF-8")
source("run_all.R", encoding = "UTF-8")
run_study("all")
```

Sourcing `run_all.R` defines the entry point. Groups run in order as `biology`,
`clinical`, `spatial` and `summary`. Spatial analyses require biology outputs.
Participant summaries require the clinical and biology outputs. Use a fresh results
directory for a complete run. The scripts record `sessionInfo_release.txt`.

Spatial analyses download public ENIGMA reference assets from commit
`b08974b55243060cbc1fad12c87048037446e8f7` to calculate parcel correspondence
and spatial rotations. [NOTICE.md](NOTICE.md) identifies the reference assets
and third-party terms.

## Verification

R syntax and script dependencies were checked. Participant summaries were
checked against existing study outputs. The complete cohort analyses were
not rerun during this repository revision.

## Contents and access

The repository contains the analysis code and input specifications, with the
anatomical field mapping required by the models. Local `data/` and `results/`
contents are excluded from Git except for their README files. The code license
and reference attributions are in `LICENSE` and `NOTICE.md`.
