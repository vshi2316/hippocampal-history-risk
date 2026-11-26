# Define the run_study() entry point.
source("config.R",local=TRUE,encoding="UTF-8")
source(file.path(.study$root,"R/preflight.R"),local=TRUE)
run_study <- function(stage="all") {
 stages <- list(
  biology = c("01_tau_raw_tables.R","02_tau_longitudinal_samples.R","03_tau_hippocampal_models.R","04_tau_participant_influence.R","05_paired_segmentation.R","06_tau_clinical_visits.R","07_regional_tau_models.R","08_regional_pipeline_restriction.R","09_memory_mri_coupling.R","10_memory_sensitivity.R","11_tau_recall_inputs.R","12_tau_recall_models.R"),
  clinical = c("13_oasis_raw_clinical_tables.R","13b_oasis_dementia_landmark.R","13c_oasis_event_definition.R","13d_oasis_acquisition_metadata.R","14_oasis_secondary_outcome.R","15_oasis_concordance.R","16_oasis_scan_crosswalk.R","17_oasis_acquisition_context.R","18_mayo_landmark.R","19_mayo_secondary_sensitivity.R","20_mayo_time_varying_effect.R","21_mayo_event_influence.R","22_primary_clinical_events.R","23_repeated_cognitive_visits.R","24_cognitive_reassessment.R","26_three_year_absolute_risk.R","27_out_of_fold_performance.R","28_confirmed_impairment.R","29_aligned_longitudinal_evidence.R","30_repeat_apoe_ventricular_inputs.R","31_matched_controls_and_screening.R","32_icv_and_repeat_decomposition.R"),
  spatial = c("33_histological_reference.R","34_histological_parcel_alignment.R","35_spatial_spin_and_cn_age.R"),
  summary = c("36_descriptive_statistics.R")
)
 selected <- if(identical(stage,"all")) names(stages) else stage
 if(any(!selected%in%names(stages))) stop("Unknown stage. Choose biology, clinical, spatial, summary or all.")
 check_environment()
 if(any(c("biology","clinical")%in%selected)) check_inputs(selected)
 dir.create(.study$output,recursive=TRUE,showWarnings=FALSE)
 for(group in selected) for(name in stages[[group]]) {
  message("Running ",group," / ",name)
  # Evaluate each script in its own environment.
  env <- new.env(parent=environment())
  sys.source(file.path(.study$root,"scripts",name),envir=env,keep.source=TRUE)
 }
 capture.output(sessionInfo(),file=file.path(.study$output,"sessionInfo_release.txt"))
 invisible(.study$output)
}
