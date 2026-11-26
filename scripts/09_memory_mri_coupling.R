# Within-person hippocampal change and delayed-recall trajectories.
local({
  packages <- c("data.table", "lme4", "lmerTest", "sandwich")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages(library(data.table))
  root <- .study$output
  input <- c(
    frozen = file.path(root, "tau_hippocampal_models/run_frozen/analysis_data/Frozen_scan_audit.tsv"),
    longitudinal = file.path(root, "tau_longitudinal_inputs/run_frozen/analysis_data/MRI_longitudinal_analysis.tsv"),
    visits = file.path(root, "clinical_intersection_26/run_frozen/analysis_data/All_valid_clinical_visits_on_frozen_timeline.tsv"),
    landmark = file.path(root, "clinical_intersection_26/run_frozen/analysis_data/Table2_person_level_clinical_intersection.tsv"))
  if (any(!file.exists(input))) stop("Missing input: ", paste(names(input)[!file.exists(input)], collapse = ", "))
  out <- file.path(root, "hippo_memory_36", paste0("run_", "frozen"))
  for (folder in c("tables", "analysis_data", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", quote = TRUE, na = "NA", bom = TRUE)
  log_line <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  save_tsv(data.table(input = names(input), path = unname(input),
                      md5 = unname(tools::md5sum(input))), "Input_manifest.tsv", "audit")

  log_line("Matching frozen MRI and delayed-recall visits")
  frozen <- fread(input[["frozen"]], encoding = "UTF-8")[
    scope %in% c("ADNI_primary", "HABS_primary"),
    .(cohort, id, mri_date, image_id, time_years, tau_z, stable_pipeline)]
  long <- fread(input[["longitudinal"]], encoding = "UTF-8")[,
    .(cohort, id, mri_date, image_id, baseline_date, tau_landmark,
      log_hip_ratio, age_c, sex_female, education_c, apoe4_carrier,
      diagnosis, amyloid_z, scanner, processing_version)]
  for (d in list(frozen, long)) {
    set(d, j = "id", value = as.character(d$id))
    set(d, j = "image_id", value = as.character(d$image_id))
    set(d, j = "mri_date", value = as.Date(d$mri_date))
  }
  long[, `:=`(baseline_date = as.Date(baseline_date), tau_landmark = as.Date(tau_landmark))]
  keys <- c("cohort", "id", "mri_date", "image_id")
  if (anyDuplicated(frozen, by = keys) || anyDuplicated(long, by = keys))
    stop("MRI scan keys must be unique")
  scans <- merge(frozen, long, by = keys, all.x = TRUE, sort = FALSE)
  if (nrow(scans) != 404L || any(!is.finite(scans$log_hip_ratio)) ||
      anyNA(scans$baseline_date) || anyNA(scans$tau_landmark))
    stop("Frozen MRI scans did not align with longitudinal hippocampal data")
  scans[, hip100 := 100 * log_hip_ratio]
  visits <- fread(input[["visits"]], encoding = "UTF-8")[
    domain == "Delayed recall", .(cohort, id, visit_date, score)]
  visits[, `:=`(id = as.character(id), visit_date = as.Date(visit_date))]
  if (anyDuplicated(visits, by = c("cohort", "id", "visit_date")))
    stop("Duplicate delayed-recall visit dates")
  matched <- merge(scans, visits, by = c("cohort", "id"), allow.cartesian = TRUE)
  matched[, gap_days := abs(as.integer(visit_date - mri_date))]
  matched <- matched[gap_days <= 180L]
  setorder(matched, cohort, id, visit_date, gap_days, mri_date, image_id)
  matched <- unique(matched, by = c("cohort", "id", "visit_date"))
  setorder(matched, cohort, id, mri_date, image_id, gap_days, visit_date)
  matched <- unique(matched, by = keys)
  paired_flow <- matched[, .(people = uniqueN(id), observations = .N), by = cohort]
  eligible <- matched[, .(baseline = any(time_years == 0),
                          followup = any(time_years > 0)), by = .(cohort, id)]
  paired <- merge(matched, eligible[baseline == TRUE & followup == TRUE,
                                    .(cohort, id)], by = c("cohort", "id"))
  baseline <- paired[time_years == 0, .(cohort, id, hip0_100 = hip100, memory0 = score)]
  if (anyDuplicated(baseline, by = c("cohort", "id"))) stop("Duplicate paired baseline MRI")
  paired <- merge(paired, baseline, by = c("cohort", "id"))
  reference <- baseline[, .(memory_mean = mean(memory0), memory_sd = sd(memory0),
                            hip_mean = mean(hip0_100), hip_sd = sd(hip0_100),
                            baseline_people = .N), by = cohort]
  if (nrow(reference) != 2L || any(reference$memory_sd <= 0) || any(reference$hip_sd <= 0))
    stop("Baseline reference scale cannot be estimated")
  paired <- merge(paired, reference, by = "cohort")
  paired[, `:=`(memory_z = (score - memory_mean) / memory_sd,
               hip0_z = (hip0_100 - hip_mean) / hip_sd,
               hip_change_pct = hip100 - hip0_100)]
  fields <- c("id", "memory_z", "time_years", "tau_z", "hip0_z", "hip_change_pct",
              "age_c", "amyloid_z", "sex_female", "education_c", "apoe4_carrier",
              "diagnosis", "scanner", "processing_version")
  paired <- paired[complete.cases(paired[, ..fields])]
  keep <- paired[, .(baseline = any(time_years == 0),
                     followup = any(time_years > 0)), by = .(cohort, id)]
  paired <- merge(paired, keep[baseline == TRUE & followup == TRUE,
                               .(cohort, id)], by = c("cohort", "id"))
  setorder(paired, cohort, id, mri_date)
  save_tsv(reference, "Table2_baseline_scales.tsv")
  save_tsv(paired[, .(cohort, id, mri_date, image_id, visit_date, gap_days,
                      time_years, score, memory_z, hip100, hip0_100,
                      hip_change_pct, tau_z, stable_pipeline)],
           "Matched_MRI_memory_visits.tsv", "analysis_data")

  fit_joint <- function(d, add_change) {
    covs <- c("sex_female", "education_c", "apoe4_carrier", "diagnosis",
              "scanner", "processing_version")
    covs <- covs[vapply(covs, function(v) uniqueN(d[[v]]) > 1L, logical(1))]
    for (v in covs) if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d, j = v, value = droplevels(factor(d[[v]])))
    rhs <- c("time_years * tau_z", "time_years * age_c",
             "time_years * amyloid_z", "hip0_z", covs)
    if (add_change) rhs <- c(rhs, "hip_change_pct")
    lmerTest::lmer(as.formula(paste("memory_z ~", paste(c(rhs, "(1 | id)"), collapse = " + "))),
                   data = d, REML = FALSE,
                   control = lme4::lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 100000)))
  }
  get_effect <- function(fit, term, co, model_name, people, observations) {
    cf <- summary(fit)$coefficients
    if (!term %in% rownames(cf)) stop("Missing coefficient: ", co, ", ", term)
    beta <- unname(cf[term, "Estimate"])
    beta_se <- unname(cf[term, "Std. Error"])
    beta_p <- unname(cf[term, "Pr(>|t|)"])
    beta_df <- unname(cf[term, "df"])
    if (!all(is.finite(c(beta, beta_se, beta_p, beta_df))))
      stop("Non-finite coefficient: ", co, ", ", term)
    critical <- qt(0.975, df = beta_df)
    data.table(cohort = co, model = model_name, term = term,
               estimate = beta, se = beta_se, ci_low = beta - critical * beta_se,
               ci_high = beta + critical * beta_se, p = beta_p,
               n_people = people, n_scans = observations,
               singular = lme4::isSingular(fit, tol = 1e-4))
  }
  log_line("Fitting joint hippocampus-memory models")
  effects <- list()
  for (co in c("ADNI", "HABS-HD")) {
    d <- copy(paired[cohort == co])
    fit0 <- fit_joint(copy(d), FALSE)
    fit1 <- fit_joint(copy(d), TRUE)
    effects[[length(effects) + 1L]] <- get_effect(
      fit0, "time_years:tau_z", co, "Tau trajectory", uniqueN(d$id), nrow(d))
    effects[[length(effects) + 1L]] <- get_effect(
      fit1, "time_years:tau_z", co, "Tau plus hippocampal change", uniqueN(d$id), nrow(d))
    effects[[length(effects) + 1L]] <- get_effect(
      fit1, "hip_change_pct", co, "Tau plus hippocampal change", uniqueN(d$id), nrow(d))
    log_line(co, ": ", uniqueN(d$id), " participants, ", nrow(d), " paired scans")
  }
  effects <- rbindlist(effects)
  save_tsv(effects, "Table3_same_person_joint_trajectory_models.tsv")

  log_line("Fitting MRI2-to-future-memory landmark models")
  person <- fread(input[["landmark"]], encoding = "UTF-8")[
    domain == "Delayed recall" & baseline_mri2 == TRUE & future_mri2 == TRUE,
    .(cohort, id, mri1_date, mri1_image_id, mri2_date, mri2_image_id,
      pre_mri2_date, pre_mri2_score, post_mri2_date, post_mri2_score)]
  person[, `:=`(id = as.character(id), mri1_date = as.Date(mri1_date),
                mri2_date = as.Date(mri2_date),
                mri1_image_id = as.character(mri1_image_id),
                mri2_image_id = as.character(mri2_image_id),
                pre_mri2_date = as.Date(pre_mri2_date),
                post_mri2_date = as.Date(post_mri2_date))]
  hip1 <- scans[, .(cohort, id, mri1_date = mri_date,
                    mri1_image_id = image_id, hip1_100 = hip100)]
  hip2 <- scans[, .(cohort, id, mri2_date = mri_date,
                    mri2_image_id = image_id, hip2_100 = hip100)]
  person <- merge(person, hip1, by = c("cohort", "id", "mri1_date", "mri1_image_id"), all.x = TRUE)
  person <- merge(person, hip2, by = c("cohort", "id", "mri2_date", "mri2_image_id"), all.x = TRUE)
  person <- merge(person, scans[time_years == 0,
                                 .(cohort, id, tau_z, age_c, sex_female)],
                  by = c("cohort", "id"), all.x = TRUE)
  person <- merge(person, reference[, .(cohort, memory_mean, memory_sd)],
                  by = "cohort", all.x = TRUE)
  if (anyDuplicated(person, by = c("cohort", "id")) ||
      anyNA(person$hip1_100) || anyNA(person$hip2_100) ||
      any(person$pre_mri2_date > person$mri2_date) ||
      any(person$post_mri2_date < person$mri2_date + 365L))
    stop("MRI2 landmark failed frozen-image or time-order audit")
  person[, `:=`(prior_hip_change_pct = hip2_100 - hip1_100,
                pre_memory_z = (pre_mri2_score - memory_mean) / memory_sd,
                post_memory_z = (post_mri2_score - memory_mean) / memory_sd,
                years_after_mri2 = as.numeric(post_mri2_date - mri2_date) / 365.25)]
  landmark_fields <- c("prior_hip_change_pct", "pre_memory_z", "post_memory_z",
                       "tau_z", "age_c", "sex_female", "years_after_mri2")
  person <- person[complete.cases(person[, ..landmark_fields])]
  save_tsv(person, "MRI2_landmark_memory_analysis.tsv", "analysis_data")
  landmark_effects <- list()
  for (co in c("ADNI", "HABS-HD")) {
    d <- person[cohort == co]
    covs <- c("pre_memory_z", "tau_z", "age_c", "sex_female", "years_after_mri2")
    covs <- covs[vapply(covs, function(v) uniqueN(d[[v]]) > 1L, logical(1))]
    fit <- lm(as.formula(paste("post_memory_z ~ prior_hip_change_pct +",
                               paste(covs, collapse = " + "))), data = d)
    beta <- unname(coef(fit)["prior_hip_change_pct"])
    beta_se <- unname(sqrt(diag(sandwich::vcovHC(fit, type = "HC3")))["prior_hip_change_pct"])
    if (!all(is.finite(c(beta, beta_se))) || beta_se <= 0)
      stop("Landmark effect is not estimable in ", co)
    beta_df <- df.residual(fit)
    critical <- qt(0.975, df = beta_df)
    landmark_effects[[co]] <- data.table(
      cohort = co, term = "prior_hip_change_pct", estimate = beta, se = beta_se,
      ci_low = beta - critical * beta_se, ci_high = beta + critical * beta_se,
      p_hc3 = 2 * pt(-abs(beta / beta_se), df = beta_df), n_people = nrow(d),
      median_years_after_mri2 = median(d$years_after_mri2))
  }
  landmark_effects <- rbindlist(landmark_effects)
  save_tsv(landmark_effects, "Table4_MRI2_to_future_memory_landmark.tsv")

  flow <- rbindlist(list(
    scans[, .(branch = "Joint trajectory", stage = "Frozen Apositive MRI",
              people = uniqueN(id), observations = .N), by = cohort],
    paired_flow[, .(cohort, branch = "Joint trajectory", stage = "MRI-memory paired",
                    people, observations)],
    paired[, .(branch = "Joint trajectory", stage = "Baseline and followup paired",
               people = uniqueN(id), observations = .N), by = cohort],
    person[, .(branch = "MRI2 landmark", stage = "Pre-MRI2 and future memory",
               people = uniqueN(id), observations = .N), by = cohort]
  ), use.names = TRUE)
  save_tsv(flow, "Table1_frozen_memory_sample_flow.tsv")
  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  writeLines(c("Frozen Stage 16 MRI identities and Stage 26 clinical visit dates retained.",
               "A clinical visit is assigned to at most one MRI within 180 days.",
               "Both mixed models use the same paired scans and differ only by hippocampal change.",
               "The MRI2 landmark is a separate, time-ordered model.",
               paste("Output:", out)), file.path(out, "36_COMPLETED.txt"), useBytes = TRUE)
  log_line("Analysis complete: ", out)
  print(flow)
  print(effects)
  print(landmark_effects)
})
