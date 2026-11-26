local({
  packages <- c("data.table", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  library(data.table)
  root <- .study$output
  step69 <- file.path(root, "clinical_field_audit_69/run_frozen")
  step70 <- file.path(root, "structural_history_70/run_frozen")
  step65 <- file.path(root, "absolute_risk_65/run_frozen")
  paths <- c(
    completed70 = file.path(step70, "70_COMPLETED.txt"),
    people = file.path(step70, "source_data/MRI2_complete_model_inputs.tsv"),
    scans = file.path(step70, "audit/OASIS_repeat_scan_metadata_matching.tsv"),
    repeats = file.path(step70, "source_data/OASIS_frozen_anchor_repeat_pairs.tsv"),
    repeated_people = file.path(step70, "source_data/OASIS_two_anchor_alternate_history.tsv"),
    anchors = file.path(step69, "audit/OASIS_selected_scan_identity.tsv"),
    anatomy = file.path(step69, "source_data/Table6_exact_scan_anatomical_pairs.tsv"),
    dictionary = file.path(step69, "audit/OASIS_MUSE_dictionary.tsv"),
    muse = file.path(root, "oasis_tau_mri_gate_27/run_frozen/audit/input_csv_cache/MUSE_volume_OASIS_04Nov2026.csv"),
    quartiles = file.path(step65, "tables/Table1_frozen_prior_loss_quartiles.tsv"),
    reference65 = file.path(step65, "tables/Table2_three_year_absolute_risk.tsv"),
    reference70 = file.path(step70, "tables/Table3_same_sample_Cox_coefficients.tsv"),
    reference_repeat70 = file.path(step70, "tables/Table13_repeat_history_Cox_coefficients.tsv")
  )
  if (any(!file.exists(paths))) stop("Missing frozen input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "icv_sensitivity_71", paste0("run_", "frozen"))
  if (dir.exists(out)) stop("Output directory already exists: ", out)
  for (folder in c("tables", "source_data", "audit")) dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, filename, folder = "tables") fwrite(as.data.table(x), file.path(out, folder, filename), sep = "\t", na = "NA")
  number <- function(x) suppressWarnings(as.numeric(x))
  horizon <- 3
  boot_reps <- 500L
  set.seed(71)
  save_tsv(data.table(input = names(paths), file = unname(paths), md5 = unname(tools::md5sum(paths))), "Input_manifest.tsv", "audit")

  cat("Step 71: verify ICV provenance and construct the pre-MRI2 participant reference\n")
  people <- fread(paths[["people"]])[cohort == "OASIS-3"]
  people[, sex := factor(as.character(sex))]
  if (nrow(people) != 306L || sum(people$event) != 59L || sum(people$event == 1L & people$time <= horizon) != 18L || anyDuplicated(people$id))
    stop("Frozen OASIS-3 primary population changed.")
  dictionary <- fread(paths[["dictionary"]])
  dict_icv <- dictionary[VARNAME == "ICV"]
  if (nrow(dict_icv) != 1L || !grepl("cubic mm", dict_icv$UNITS, fixed = TRUE)) stop("Released ICV unit does not match the inspected dictionary.")
  save_tsv(dict_icv, "ICV_dictionary_and_units.tsv", "audit")
  anchors <- fread(paths[["anchors"]])
  if (nrow(anchors) != 612L || anyDuplicated(anchors, by = c("id", "scan"))) stop("Frozen MRI anchors changed.")
  anchor_wide <- dcast(anchors, id ~ scan, value.var = c("Code", "anchor_day", "hippocampus", "ICV"))
  if (!setequal(anchor_wide$id, people$id)) stop("Anchor identities differ from the primary population.")
  scans <- fread(paths[["scans"]])
  muse_source <- fread(paths[["muse"]])
  if (anyDuplicated(scans$Code) || anyDuplicated(muse_source$Code)) stop("Non-unique MUSE Code.")
  source_small <- muse_source[, .(Code, source_present = TRUE, source_ICV = number(ICV), source_hippo = number(Left_Hippocampus) + number(Right_Hippocampus),
    source_day = as.integer(as.Date(EXAMDATE) - as.Date("2000-01-01")))]
  provenance <- merge(scans, source_small, by = "Code", all.x = TRUE)
  if (any(!provenance$source_present %in% TRUE) ||
      any(is.na(provenance$icv) != is.na(provenance$source_ICV)) ||
      any(is.na(provenance$hippo) != is.na(provenance$source_hippo)) ||
      any(is.na(provenance$mri_day) != is.na(provenance$source_day)) ||
      any(abs(provenance$icv - provenance$source_ICV) > 1e-6, na.rm = TRUE) ||
      any(abs(provenance$hippo - provenance$source_hippo) > 1e-6, na.rm = TRUE) || any(provenance$mri_day != provenance$source_day, na.rm = TRUE))
    stop("Step 70 scan values or dates disagree with the pinned MUSE source.")
  save_tsv(provenance, "Released_MUSE_scan_provenance.tsv", "audit")

  reference_icv <- function(scan_data, anchor_data) {
    candidate <- merge(scan_data, anchor_data[, .(id, start_day = anchor_day_MRI1, end_day = anchor_day_MRI2)], by = "id")
    candidate[, eligible := valid_code_day %in% TRUE & exact_metadata %in% TRUE & is.finite(icv) & icv > 0 &
      is.finite(mri_day) & mri_day >= start_day & mri_day <= end_day]
    daily <- candidate[eligible == TRUE, .(daily_ICV = as.numeric(median(icv)), scans_used = .N,
      ICV_min = min(icv), ICV_max = max(icv), Codes = paste(sort(Code), collapse = " | ")), by = .(id, mri_day)]
    reference <- daily[, .(reference_ICV = as.numeric(median(daily_ICV)), dates_used = .N, acquisitions_used = sum(scans_used),
      first_day_used = min(mri_day), last_day_used = max(mri_day)), by = id]
    list(candidates = candidate, daily = daily, reference = reference)
  }
  ref <- reference_icv(scans, anchor_wide)
  save_tsv(ref$candidates, "ICV_reference_scan_candidates.tsv", "audit")
  save_tsv(ref$daily, "Daily_ICV_reference_inputs.tsv", "source_data")
  save_tsv(ref$reference, "Participant_pre_MRI2_ICV_reference.tsv", "source_data")
  if (nrow(ref$reference) != 306L || !setequal(ref$reference$id, people$id) || any(!is.finite(ref$reference$reference_ICV)))
    stop("The frozen reference rule does not cover all 306 participants. Inspect the saved candidate audit; do not silently change the risk set.")
  d <- merge(merge(people, anchor_wide, by = "id"), ref$reference, by = "id")
  if (any(abs(100 * d$hippocampus_MRI2 / d$ICV_MRI2 - d$hippo2_pct_ICV) > 1e-6) ||
      any(abs(100 * log(d$hippocampus_MRI1 / d$hippocampus_MRI2) / d$interval_years - d$hippo_loss_pct_year) > 1e-6))
    stop("Frozen anchors do not reconstruct the existing hippocampal state/history.")
  d[, `:=`(stable_hippo2_pct_ICV = 100 * hippocampus_MRI2 / reference_ICV,
    current_ICV_log_difference_from_reference_pct = 100 * log(ICV_MRI2 / reference_ICV))]
  vent <- fread(paths[["anatomy"]])[cohort == "OASIS-3" & field == "lateral_ventricles" & usable_pair %in% TRUE,
    .(id, ventricular_MRI2_mm3 = value_MRI2)]
  if (anyDuplicated(vent$id)) stop("Duplicate ventricular records.")
  d <- merge(d, vent, by = "id", all.x = TRUE)
  d[, stable_ventricle2_pct_ICV := 100 * ventricular_MRI2_mm3 / reference_ICV]
  if (any(!is.finite(d$stable_hippo2_pct_ICV)) || any(!is.finite(d$stable_ventricle2_pct_ICV))) stop("Incomplete normalization sensitivity inputs.")
  save_tsv(d, "OASIS_frozen_and_reference_normalized_inputs.tsv", "source_data")
  save_tsv(d[, .(participants = .N, events_total = sum(event), events_by_3y = sum(event == 1L & time <= horizon),
    reference_ICV_median_mm3 = median(reference_ICV), dates_used_min = min(dates_used), dates_used_median = median(dates_used),
    acquisitions_used_median = median(acquisitions_used), max_reference_day_after_MRI2 = max(last_day_used - anchor_day_MRI2))],
    "Table1_reference_ICV_coverage.tsv")
  repeat_pairs <- fread(paths[["repeats"]])[status == "MATCHED_REPEAT"]
  repeat_pairs[, `:=`(ICV_log_difference_pct = 100 * log(alternate_ICV / frozen_ICV),
    ICV_arithmetic_change_pct = 100 * (alternate_ICV / frozen_ICV - 1))]
  repeat_pairs <- merge(repeat_pairs, d[, .(id, reference_ICV)], by = "id", all.x = TRUE)
  repeat_pairs[, common_reference_state_log_difference_pct := 100 * log((alternate_hippo / reference_ICV) / (frozen_hippo / reference_ICV))]
  if (any(abs(repeat_pairs$common_reference_state_log_difference_pct - repeat_pairs$raw_hippo_log_difference_pct) > 1e-8)) stop("Common-denominator identity failed.")
  save_tsv(repeat_pairs, "Same_day_ICV_and_volume_provenance_review.tsv", "source_data")
  save_tsv(repeat_pairs[, .(anchor_pairs = .N, participants = uniqueN(id),
    median_absolute_ICV_log_difference_pct = median(abs(ICV_log_difference_pct)),
    ICV_log_difference_SD_pct = sd(ICV_log_difference_pct), ICV_abs_log_difference_over_5pct = sum(abs(ICV_log_difference_pct) > 5),
    ICV_abs_log_difference_over_10pct = sum(abs(ICV_log_difference_pct) > 10),
    ICV_abs_log_difference_over_20pct = sum(abs(ICV_log_difference_pct) > 20),
    raw_hippo_log_difference_SD_pct = sd(raw_hippo_log_difference_pct),
    scan_specific_normalized_SD_pct = sd(normalized_hippo_log_difference_pct),
    common_reference_normalized_SD_pct = sd(common_reference_state_log_difference_pct))], "Table2_same_day_ICV_discrepancy.tsv")

  quartile <- fread(paths[["quartiles"]])[cohort == "OASIS-3"]
  if (nrow(quartile) != 1L || !is.finite(quartile$loss_Q1) || !is.finite(quartile$loss_Q3) || quartile$loss_Q1 >= quartile$loss_Q3)
    stop("Invalid frozen OASIS-3 loss quartiles.")
  q <- c(q1 = quartile$loss_Q1[[1L]], q3 = quartile$loss_Q3[[1L]])
  main_formula <- survival::Surv(time, event) ~ age_MRI2 + sex + education + hippo2_pct_ICV + baseline_score + hippo_loss_pct_year
  apoe_formula <- update(main_formula, . ~ . + apoe_e4)
  vent_formula <- update(apoe_formula, . ~ . + ventricle2_pct_ICV + ventricle_expansion_pct_year)
  fit_cox <- function(form, x) {
    fit <- withCallingHandlers(survival::coxph(form, data = droplevels(as.data.frame(x)), ties = "efron", x = TRUE, model = TRUE),
      warning = function(w) {
        if (grepl("infinite|converg|singular", conditionMessage(w), ignore.case = TRUE)) stop(conditionMessage(w))
      })
    if (any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit)))) stop("Non-finite Cox coefficient or variance: ",
      paste(names(coef(fit))[!is.finite(coef(fit))], collapse = ", "))
    fit
  }
  risk3 <- function(fit, x) {
    x <- as.data.frame(x); reference <- x[1L, , drop = FALSE]
    sf <- survival::survfit(fit, newdata = reference, se.fit = FALSE)
    ix <- which(sf$time <= horizon)
    href <- if (length(ix)) as.numeric(sf$cumhaz)[max(ix)] else 0
    if (href == 0) return(rep(0, nrow(x)))
    lp <- as.numeric(predict(fit, newdata = x, type = "lp", reference = "zero"))
    rlp <- as.numeric(predict(fit, newdata = reference, type = "lp", reference = "zero"))
    -expm1(-exp(log(href) + lp - rlp))
  }
  metrics <- function(fit, x) {
    low <- copy(x); high <- copy(x)
    low[, hippo_loss_pct_year := as.numeric(q[["q1"]])]
    high[, hippo_loss_pct_year := as.numeric(q[["q3"]])]
    r1 <- mean(risk3(fit, low)); r3 <- mean(risk3(fit, high))
    beta <- unname(coef(fit)["hippo_loss_pct_year"])
    c(beta_loss = beta, HR_per_1pct_year = exp(beta), risk_Q1 = r1, risk_Q3 = r3, risk_difference = r3 - r1)
  }
  limits <- function(v) {
    good <- v[is.finite(v)]
    if (!length(good) || length(good) < .9 * length(v)) return(c(NA_real_, NA_real_))
    as.numeric(quantile(good, c(.025, .975), names = FALSE))
  }
  bootstrap_failures <- list()
  analyze_definitions <- function(sets, form, population, contrasts) {
    x <- sets[[1L]]
    for (z in sets) if (!identical(z$id, x$id) || !identical(z$time, x$time) || !identical(z$event, x$event)) stop("Compared definitions have different participants or follow-up.")
    models <- lapply(sets, function(z) fit_cox(form, z))
    point <- lapply(seq_along(sets), function(k) metrics(models[[k]], sets[[k]]))
    names(point) <- names(sets)
    metric_names <- names(point[[1L]])
    boot <- lapply(sets, function(z) matrix(NA_real_, nrow = boot_reps, ncol = length(metric_names), dimnames = list(NULL, metric_names)))
    for (b in seq_len(boot_reps)) {
      draw <- sample.int(nrow(x), nrow(x), replace = TRUE)
      for (k in seq_along(sets)) {
        z <- sets[[k]][draw]
        boot[[k]][b, ] <- tryCatch(metrics(fit_cox(form, z), z), error = function(e) {
          bootstrap_failures[[length(bootstrap_failures) + 1L]] <<- data.table(population, definition = names(sets)[k], replicate = b, error = conditionMessage(e))
          rep(NA_real_, length(metric_names))
        })
      }
    }
    summary_rows <- rbindlist(lapply(seq_along(sets), function(k) rbindlist(lapply(seq_along(metric_names), function(j) {
      interval <- limits(boot[[k]][, j])
      data.table(population, definition = names(sets)[k], participants = nrow(x), events_total = sum(x$event),
        events_by_3y = sum(x$event == 1L & x$time <= horizon), loss_Q1 = q[["q1"]], loss_Q3 = q[["q3"]],
        metric = metric_names[j], estimate = point[[k]][j], CI_low = interval[1L], CI_high = interval[2L],
        bootstrap_requested = boot_reps, bootstrap_valid = sum(is.finite(boot[[k]][, j])),
        CI_status = if (all(is.finite(interval))) "AVAILABLE" else "BOOTSTRAP_STABILITY_REVIEW")
    }))))
    paired_rows <- rbindlist(lapply(seq_len(nrow(contrasts)), function(k) {
      base <- contrasts$base[k]; alternative <- contrasts$alternative[k]
      if (!all(c(base, alternative) %in% names(sets))) stop("Unknown paired contrast.")
      rbindlist(lapply(seq_along(metric_names), function(j) {
        difference <- boot[[alternative]][, j] - boot[[base]][, j]
        interval <- limits(difference)
        data.table(population, base, alternative, participants = nrow(x), events_total = sum(x$event), metric = metric_names[j],
          paired_change = point[[alternative]][j] - point[[base]][j], CI_low = interval[1L], CI_high = interval[2L],
          bootstrap_requested = boot_reps, bootstrap_valid = sum(is.finite(difference)),
          CI_status = if (all(is.finite(interval))) "AVAILABLE" else "BOOTSTRAP_STABILITY_REVIEW")
      }))
    }))
    coefficient_rows <- rbindlist(lapply(seq_along(sets), function(k) {
      s <- summary(models[[k]])$coefficients
      data.table(population, definition = names(sets)[k], participants = nrow(x), events_total = sum(x$event), term = rownames(s),
        beta = s[, "coef"], SE = s[, "se(coef)"], HR = exp(s[, "coef"]),
        CI_low = exp(s[, "coef"] - 1.96 * s[, "se(coef)"]), CI_high = exp(s[, "coef"] + 1.96 * s[, "se(coef)"]), P = s[, "Pr(>|z|)"])
    }))
    list(summary = summary_rows, contrasts = paired_rows, coefficients = coefficient_rows)
  }
  # Verify the reference models before participant bootstrap resampling.
  reference65 <- fread(paths[["reference65"]])[cohort == "OASIS-3" & landmark == "MRI2"]
  reference70 <- fread(paths[["reference70"]])
  reference_repeat70 <- fread(paths[["reference_repeat70"]])
  reproduced <- list()
  original_check <- function(x, form, label, target_HR, target_RD = NA_real_) {
    m <- metrics(fit_cox(form, x), x)
    if (length(target_HR) != 1L || !is.finite(target_HR) || abs(m[["HR_per_1pct_year"]] - target_HR) > 1e-6 ||
        (is.finite(target_RD) && abs(m[["risk_difference"]] - target_RD) > 1e-6)) stop("Frozen model did not reproduce: ", label)
    data.table(model = label, participants = nrow(x), events_total = sum(x$event), reproduced_HR = m[["HR_per_1pct_year"]], reference_HR = target_HR,
      reproduced_RD = m[["risk_difference"]], reference_RD = target_RD, status = "REPRODUCED")
  }
  if (nrow(reference65) != 1L) stop("Invalid Step 65 reference row.")
  reproduced[[1L]] <- original_check(d, main_formula, "Step65_OASIS_MRI2", reference65$HR_per_1pct_year, reference65$risk_difference)
  g <- d[apoe_e4 %in% c(0, 1)]
  if (nrow(g) != 304L || sum(g$event) != 58L) stop("Step 70 APOE-complete population changed.")
  for (comparison_name in c("APOE_adjustment", "Lateral_ventricle_adjustment")) {
    target <- reference70[cohort == "OASIS-3" & comparison == comparison_name & model == "extended" & term == "hippo_loss_pct_year", HR]
    reproduced[[length(reproduced) + 1L]] <- original_check(g, if (comparison_name == "APOE_adjustment") apoe_formula else vent_formula, comparison_name, target)
  }
  save_tsv(rbindlist(reproduced), "Frozen_model_reproduction.tsv", "audit")

  cat("Step 71: compare scan-specific and participant-reference normalization\n")
  stable_version <- function(z, normalize_ventricles = FALSE) {
    z <- copy(z)
    z[, hippo2_pct_ICV := stable_hippo2_pct_ICV]
    if (normalize_ventricles) z[, ventricle2_pct_ICV := stable_ventricle2_pct_ICV]
    z
  }
  pair <- data.table(base = "Original_scan_ICV", alternative = "Participant_reference_ICV")
  main <- analyze_definitions(list(Original_scan_ICV = d, Participant_reference_ICV = stable_version(d)), main_formula, "Primary_306", pair)
  genetic <- analyze_definitions(list(Original_scan_ICV = g, Participant_reference_ICV = stable_version(g)), apoe_formula, "APOE_complete_304", pair)
  ventricular <- analyze_definitions(list(Original_scan_ICV = g, Participant_reference_ICV = stable_version(g, TRUE)), vent_formula, "APOE_and_ventricles_304", pair)
  main_runs <- list(main, genetic, ventricular)
  save_tsv(rbindlist(lapply(main_runs, `[[`, "summary")), "Table3_normalization_sensitivity_estimates.tsv")
  save_tsv(rbindlist(lapply(main_runs, `[[`, "contrasts")), "Table4_paired_normalization_changes.tsv")
  save_tsv(rbindlist(lapply(main_runs, `[[`, "coefficients")), "Table5_normalization_Cox_coefficients.tsv")

  cat("Step 71: separate repeat-history and current-state changes in the same people\n")
  repeated <- fread(paths[["repeated_people"]])[, .(id, alternate_hippo_MRI1, alternate_hippo_MRI2, alternate_ICV_MRI2, alternate_loss_pct_year)]
  if (anyDuplicated(repeated$id)) stop("Duplicate repeat-subset participants.")
  xr <- merge(g, repeated, by = "id")
  if (nrow(xr) != 70L || sum(xr$event) != 14L) stop("The frozen two-anchor repeat population changed.")
  if (any(abs(100 * log(xr$alternate_hippo_MRI1 / xr$alternate_hippo_MRI2) / xr$interval_years - xr$alternate_loss_pct_year) > 1e-6))
    stop("Alternate hippocampal history does not reconstruct.")
  target <- reference_repeat70[model == "base" & term == "hippo_loss_pct_year", HR]
  reproduced[[length(reproduced) + 1L]] <- original_check(xr, apoe_formula, "Step70_original_repeat_subset", target)
  save_tsv(rbindlist(reproduced), "Frozen_model_reproduction.tsv", "audit")
  variants <- list(Original_scan_ICV = xr, Stable_frozen_history_and_state = stable_version(xr))
  variants$History_only <- copy(variants$Stable_frozen_history_and_state)
  variants$History_only[, hippo_loss_pct_year := alternate_loss_pct_year]
  variants$Current_state_only <- copy(variants$Stable_frozen_history_and_state)
  variants$Current_state_only[, hippo2_pct_ICV := 100 * alternate_hippo_MRI2 / reference_ICV]
  variants$History_and_current_state <- copy(variants$Current_state_only)
  variants$History_and_current_state[, hippo_loss_pct_year := alternate_loss_pct_year]
  contrasts <- data.table(base = c("Original_scan_ICV", rep("Stable_frozen_history_and_state", 3L)),
    alternative = c("Stable_frozen_history_and_state", "History_only", "Current_state_only", "History_and_current_state"))
  repeat_run <- analyze_definitions(variants, apoe_formula, "Two_anchor_repeat_70", contrasts)
  save_tsv(repeat_run$summary, "Table6_repeat_history_state_decomposition.tsv")
  save_tsv(repeat_run$contrasts, "Table7_paired_repeat_decomposition_changes.tsv")
  save_tsv(repeat_run$coefficients, "Table8_repeat_decomposition_Cox_coefficients.tsv")
  save_tsv(rbindlist(lapply(names(variants), function(nm) variants[[nm]][, .(definition = nm, id, time, event,
    hippo2_pct_ICV, hippo_loss_pct_year, reference_ICV)])), "Repeat_definition_person_level_inputs.tsv", "source_data")
  if (length(bootstrap_failures)) save_tsv(rbindlist(bootstrap_failures), "Bootstrap_failures.tsv", "audit")
  capture.output(sessionInfo(), file = file.path(out, "audit/R_session_info.txt"))
  writeLines(c("STEP_71_COMPLETED", normalizePath(out, winslash = "/", mustWork = TRUE)), file.path(out, "71_COMPLETED.txt"), useBytes = TRUE)
  cat("Step 71 complete: ", out, "\n", sep = "")
  print(fread(file.path(out, "tables/Table3_normalization_sensitivity_estimates.tsv"))[metric %in% c("HR_per_1pct_year", "risk_difference"),
    .(population, definition, participants, metric, estimate, CI_low, CI_high, CI_status)])
  print(fread(file.path(out, "tables/Table6_repeat_history_state_decomposition.tsv"))[metric == "HR_per_1pct_year",
    .(definition, participants, estimate, CI_low, CI_high, bootstrap_valid, CI_status)])
})
