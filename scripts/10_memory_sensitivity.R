# Scanner restrictions and participant influence in the hippocampus-memory model.
local({
  packages <- c("data.table", "lme4", "lmerTest")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages(library(data.table))
  root <- .study$output
  input <- c(
    paired = file.path(root, "hippo_memory_36/run_frozen/analysis_data/Matched_MRI_memory_visits.tsv"),
    scales = file.path(root, "hippo_memory_36/run_frozen/tables/Table2_baseline_scales.tsv"),
    original = file.path(root, "hippo_memory_36/run_frozen/tables/Table3_same_person_joint_trajectory_models.tsv"),
    longitudinal = file.path(root, "tau_longitudinal_inputs/run_frozen/analysis_data/MRI_longitudinal_analysis.tsv"))
  if (any(!file.exists(input))) stop("Missing input: ", paste(names(input)[!file.exists(input)], collapse = ", "))
  out <- file.path(root, "hippo_memory_37", paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", quote = TRUE, na = "NA", bom = TRUE)
  log_line <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  save_tsv(data.table(input = names(input), path = unname(input),
                      md5 = unname(tools::md5sum(input))), "Input_manifest.tsv", "audit")

  matched <- fread(input[["paired"]], encoding = "UTF-8")
  long <- fread(input[["longitudinal"]], encoding = "UTF-8")[,
    .(cohort, id, mri_date, image_id, age_c, sex_female, education_c,
      apoe4_carrier, diagnosis, amyloid_z, scanner, processing_version)]
  scales <- fread(input[["scales"]], encoding = "UTF-8")[,
    .(cohort, hip_mean, hip_sd)]
  for (d in list(matched, long)) {
    set(d, j = "id", value = as.character(d$id))
    set(d, j = "image_id", value = as.character(d$image_id))
    set(d, j = "mri_date", value = as.Date(d$mri_date))
  }
  keys <- c("cohort", "id", "mri_date", "image_id")
  if (anyDuplicated(matched, by = keys) || anyDuplicated(long, by = keys))
    stop("MRI scan keys must be unique")
  d <- merge(matched, long, by = keys, all.x = TRUE, sort = FALSE)
  if (nrow(d) != nrow(matched) || anyNA(d$age_c))
    stop("Stage 36 paired scans did not match Stage 11 covariates")
  d <- merge(d, scales, by = "cohort")
  d[, hip0_z := (hip0_100 - hip_mean) / hip_sd]
  d[, stable_pipeline := as.logical(stable_pipeline)]
  fields <- c("memory_z", "hip0_z", "hip_change_pct", "time_years",
              "tau_z", "age_c", "amyloid_z", "sex_female", "education_c",
              "apoe4_carrier", "diagnosis", "scanner", "processing_version")
  if (any(!complete.cases(d[, ..fields]))) stop("Stage 36 model data are incomplete")
  original <- fread(input[["original"]], encoding = "UTF-8")[
    model == "Tau plus hippocampal change" &
      term %in% c("hip_change_pct", "time_years:tau_z")]

  fit_fixed <- function(x, co, scenario_name) {
    x <- copy(x)
    people <- uniqueN(x$id)
    scans <- nrow(x)
    terms <- c("hip_change_pct", "time_years:tau_z")
    empty <- data.table(cohort = co, scenario = scenario_name, term = terms,
                        estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
                        ci_high = NA_real_, p = NA_real_, n_people = people,
                        n_scans = scans, singular = NA, status = "model_error")
    covs <- c("sex_female", "education_c", "apoe4_carrier", "diagnosis",
              "scanner", "processing_version")
    covs <- covs[vapply(covs, function(v) uniqueN(x[[v]]) > 1L, logical(1))]
    for (v in covs) if (is.character(x[[v]]) || is.factor(x[[v]]))
      set(x, j = v, value = droplevels(factor(x[[v]])))
    rhs <- c("time_years * tau_z", "time_years * age_c",
             "time_years * amyloid_z", "hip0_z", covs,
             "hip_change_pct", "(1 | id)")
    formula <- as.formula(paste("memory_z ~", paste(rhs, collapse = " + ")))
    fit <- tryCatch(suppressWarnings(lmerTest::lmer(
      formula, data = x, REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 100000)))), error = identity)
    if (inherits(fit, "error")) {
      empty[, status := conditionMessage(fit)]
      return(empty)
    }
    cf <- summary(fit)$coefficients
    if (!all(terms %in% rownames(cf))) {
      empty[, status := "target_coefficient_not_estimable"]
      return(empty)
    }
    for (i in seq_along(terms)) {
      term_name <- terms[i]
      beta <- unname(cf[term_name, "Estimate"])
      beta_se <- unname(cf[term_name, "Std. Error"])
      beta_p <- unname(cf[term_name, "Pr(>|t|)"])
      beta_df <- unname(cf[term_name, "df"])
      if (!all(is.finite(c(beta, beta_se, beta_p, beta_df)))) {
        set(empty, i = i, j = "status", value = "nonfinite_coefficient")
        next
      }
      critical <- qt(0.975, beta_df)
      set(empty, i = i,
          j = c("estimate", "se", "ci_low", "ci_high", "p", "singular", "status"),
          value = list(beta, beta_se, beta - critical * beta_se,
                       beta + critical * beta_se, beta_p,
                       lme4::isSingular(fit, tol = 1e-4), "fitted"))
    }
    empty
  }

  log_line("Reproducing Stage 36 full-cohort effects")
  full <- rbindlist(lapply(c("ADNI", "HABS-HD"), function(co)
    fit_fixed(d[cohort == co], co, "Full")))
  check <- merge(full, original[, .(cohort, term, original_estimate = estimate,
                                    original_n_people = n_people,
                                    original_n_scans = n_scans)],
                 by = c("cohort", "term"))
  check[, absolute_difference := abs(estimate - original_estimate)]
  save_tsv(check, "Table1_Stage36_exact_model_reproduction.tsv")
  if (nrow(check) != 4L || any(check$status != "fitted") ||
      any(check$absolute_difference > 1e-6) ||
      any(check$n_people != check$original_n_people) ||
      any(check$n_scans != check$original_n_scans))
    stop("The Stage 36 fixed model was not reproduced; see Table 1")

  habs <- d[cohort == "HABS-HD"]
  stable_ids <- habs[stable_pipeline == TRUE,
                     .(baseline = any(time_years == 0),
                       followup = any(time_years > 0)), by = id][
                         baseline == TRUE & followup == TRUE, id]
  scenarios <- list(
    list(co = "ADNI", name = "Full", data = d[cohort == "ADNI"]),
    list(co = "HABS-HD", name = "Full", data = habs),
    list(co = "HABS-HD", name = "Without HABS_4739", data = habs[id != "HABS_4739"]),
    list(co = "HABS-HD", name = "Stable people, all scans", data = habs[id %in% stable_ids]),
    list(co = "HABS-HD", name = "Stable people, stable scans",
         data = habs[id %in% stable_ids & stable_pipeline == TRUE]))
  scenario_effects <- rbindlist(lapply(scenarios, function(s)
    if (s$name == "Full") full[cohort == s$co] else
      fit_fixed(s$data, s$co, s$name)))
  save_tsv(scenario_effects, "Table2_fixed_model_sample_and_scanner_sensitivity.tsv")
  log_line("Fitting participant leave-one-out models")
  loo_rows <- list()
  index <- 0L
  for (co in c("ADNI", "HABS-HD")) {
    cohort_data <- d[cohort == co]
    for (excluded in unique(cohort_data$id)) {
      index <- index + 1L
      z <- fit_fixed(cohort_data[id != excluded], co, "Leave one out")
      z <- z[term == "hip_change_pct"]
      z[, excluded_id := excluded]
      loo_rows[[index]] <- z
      if (index %% 20L == 0L) log_line("Leave-one-out models: ", index, "/149")
    }
  }
  loo <- rbindlist(loo_rows)
  full_beta <- full[term == "hip_change_pct", .(cohort, full_estimate = estimate)]
  loo <- merge(loo, full_beta, by = "cohort")
  loo[, difference_from_full := estimate - full_estimate]
  save_tsv(loo, "Table3_all_participant_leave_one_out.tsv")
  influence <- loo[status == "fitted", .(
    leave_one_out_fits = .N,
    minimum_estimate = min(estimate), maximum_estimate = max(estimate),
    range_of_difference = max(abs(difference_from_full)),
    positive_estimates = sum(estimate > 0),
    intervals_above_zero = sum(ci_low > 0),
    intervals_below_zero = sum(ci_high < 0)), by = cohort]
  save_tsv(influence, "Table4_leave_one_out_summary.tsv")
  top <- copy(loo[status == "fitted"])
  top[, absolute_difference := abs(difference_from_full)]
  setorder(top, cohort, -absolute_difference)
  save_tsv(top[, head(.SD, 8L), by = cohort], "Table5_largest_participant_influences.tsv")
  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  writeLines(c("Stage 36 coefficients reproduced before any sensitivity fit.",
               "The outcome, covariates, random effect, memory scaling and MRI-visit matching are unchanged.",
               "No participant is deleted from the frozen original result.",
               paste("Output:", out)), file.path(out, "37_COMPLETED.txt"), useBytes = TRUE)
  log_line("Analysis complete: ", out)
  print(scenario_effects[term == "hip_change_pct"])
  print(influence)
})
