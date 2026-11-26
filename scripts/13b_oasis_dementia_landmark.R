# Verify OASIS-3 dementia events in the MRI2 risk set.

needed <- c("data.table", "survival")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))

main <- function() {
  code_dir <- .study$output
  input <- file.path(code_dir, "44_OASIS3_full_clinical_link_results")
  out <- file.path(code_dir, "45_OASIS3_frozen_event_link_results")
  required <- c("00_status.txt", "04_MUSE_QC_and_provenance.csv",
                "07_D1_all_visit_stages.csv", "09_B4_CDR_by_visit.csv",
                "12_two_MRI_clinical_landmarks.csv", "19b_event_risk_set.csv")
  paths <- file.path(input, required)
  if (any(!file.exists(paths)))
    stop("Missing Step 44 output: ", paste(required[!file.exists(paths)], collapse = ", "))
  status44 <- readLines(paths[1], warn = FALSE)
  if (!any(status44 == "STATUS=STEP44_COMPLETE")) stop("Step 44 is incomplete")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines("STATUS=STEP45_INCOMPLETE", file.path(out, "00_status.txt"))
  save_csv <- function(x, name) data.table::fwrite(x, file.path(out, name), na = "")

  two <- data.table::fread(paths[5])
  risk <- data.table::fread(paths[6])
  muse <- data.table::fread(paths[2])
  dx <- data.table::fread(paths[3])
  cdr <- data.table::fread(paths[4])
  if (anyDuplicated(two$OASISID) || anyDuplicated(risk$OASISID))
    stop("Step 44 landmark IDs are not unique")
  if (nrow(risk) != 317L || sum(risk$event, na.rm = TRUE) != 53L)
    stop("Step 44 event risk set changed; inspect 19b_event_risk_set.csv")
  extra <- two[, .(OASISID, MRI1_day, MRI1_Code, MRI2_Code,
                   hippo2_pct_ICV, sex, education, baseline_D1_day,
                   baseline_CDR_day, D1_gap_days, CDR_gap_days,
                   first_dementia_day, last_D1_day)]
  risk <- merge(risk, extra, by = "OASISID", all.x = TRUE, sort = FALSE)
  if (nrow(risk) != 317L || anyNA(risk$MRI1_day))
    stop("Step 44 risk and landmark rows do not join one to one")
  risk[, event := as.logical(event)]
  risk[, sex := factor(sex)]
  cn <- risk[baseline_stage == "CN"]
  if (nrow(cn) != 312L || sum(cn$event) != 52L)
    stop("The frozen CN event set changed")
  save_csv(risk, "01_frozen_event_risk_set.csv")

  # Audit the first recorded event against adjacent diagnostic and CDR visits.
  event_rows <- risk[event == TRUE]
  event_audit <- data.table::rbindlist(lapply(seq_len(nrow(event_rows)), function(i) {
    r <- event_rows[i]
    visits <- dx[OASISID == r$OASISID]
    before <- visits[day < r$first_dementia_day &
                       stage %in% c("CN", "MCI", "Dementia")]
    after <- visits[day > r$first_dementia_day &
                      stage %in% c("CN", "MCI", "Dementia")]
    near_cdr <- cdr[OASISID == r$OASISID &
                      abs(day - r$first_dementia_day) <= 180L]
    if (nrow(near_cdr))
      near_cdr <- near_cdr[order(abs(near_cdr$day - r$first_dementia_day))]
    prior_day <- if (nrow(before)) max(before$day) else NA_integer_
    prior_stage <- if (nrow(before)) before$stage[which.max(before$day)] else NA_character_
    next_day <- if (nrow(after)) min(after$day) else NA_integer_
    next_stage <- if (nrow(after)) after$stage[which.min(after$day)] else NA_character_
    data.table::data.table(
      OASISID = r$OASISID, baseline_stage = r$baseline_stage,
      first_dementia_day = r$first_dementia_day,
      years_after_MRI2 = r$followup_days / 365.25,
      prior_D1_day = prior_day, prior_D1_stage = prior_stage,
      diagnosis_interval_days = r$first_dementia_day - prior_day,
      next_D1_day = next_day, next_D1_stage = next_stage,
      next_D1_within_18m = !is.na(next_day) &&
        next_day - r$first_dementia_day <= 548L,
      CDR_day_near_event = if (nrow(near_cdr)) near_cdr$day[1] else NA_integer_,
      CDRSB_near_event = if (nrow(near_cdr)) near_cdr$cdrsb[1] else NA_real_,
      review_flag = (nrow(near_cdr) > 0L && near_cdr$cdrsb[1] == 0) ||
        (!is.na(next_day) && next_day - r$first_dementia_day <= 548L &&
           next_stage %in% c("CN", "MCI")))
  }), fill = TRUE)
  save_csv(event_audit, "02_event_diagnosis_audit.csv")
  save_csv(event_audit[review_flag == TRUE], "02b_diagnosis_review_queue.csv")
  diagnosis_summary <- data.table::data.table(
    metric = c("First recorded dementia events", "Events with CDR within 180 days",
               "Events with CDR-SB zero near diagnosis",
               "Events with another D1 visit within 18 months",
               "Next D1 returns to CN or MCI within 18 months",
               "Diagnosis interval longer than 2 years"),
    count = c(nrow(event_audit), sum(!is.na(event_audit$CDRSB_near_event)),
              sum(event_audit$CDRSB_near_event == 0, na.rm = TRUE),
              sum(event_audit$next_D1_within_18m),
              sum(event_audit$next_D1_within_18m &
                    event_audit$next_D1_stage %in% c("CN", "MCI")),
              sum(event_audit$diagnosis_interval_days > 730L, na.rm = TRUE)))
  save_csv(diagnosis_summary, "03_event_diagnosis_summary.csv")

  # Same-day MUSE runs provide an internal check of measurement stability.
  muse <- muse[eligible == TRUE]
  muse[, hippocampus := as.numeric(Left_Hippocampus) +
                          as.numeric(Right_Hippocampus)]
  runs <- muse[, .(
    n_runs = .N,
    median_hippo = stats::median(hippocampus),
    range_pct = 100 * (max(hippocampus) - min(hippocampus)) /
      stats::median(hippocampus)), by = .(OASISID, mri_day)]
  first_runs <- data.table::copy(runs)
  data.table::setnames(first_runs,
    c("mri_day", "n_runs", "median_hippo", "range_pct"),
    c("MRI1_day", "MRI1_n_runs", "MRI1_median_hippo", "MRI1_range_pct"))
  second_runs <- data.table::copy(runs)
  data.table::setnames(second_runs,
    c("mri_day", "n_runs", "median_hippo", "range_pct"),
    c("MRI2_day", "MRI2_n_runs", "MRI2_median_hippo", "MRI2_range_pct"))
  risk <- merge(risk, first_runs, by = c("OASISID", "MRI1_day"),
                all.x = TRUE, sort = FALSE)
  risk <- merge(risk, second_runs, by = c("OASISID", "MRI2_day"),
                all.x = TRUE, sort = FALSE)
  if (nrow(risk) != 317L || anyNA(risk$MRI1_median_hippo) ||
      anyNA(risk$MRI2_median_hippo))
    stop("MUSE run-level audit did not match all MRI dates")
  risk[, median_run_loss_pct_year :=
    100 * log(MRI1_median_hippo / MRI2_median_hippo) /
      ((MRI2_day - MRI1_day) / 365.25)]
  risk[, max_same_day_range_pct := pmax(MRI1_range_pct, MRI2_range_pct)]
  save_csv(risk[, .(OASISID, baseline_stage, event, MRI1_Code, MRI2_Code,
                    MRI1_n_runs, MRI2_n_runs, MRI1_range_pct, MRI2_range_pct,
                    hippo_loss_pct_year, median_run_loss_pct_year)],
           "04_same_day_run_reliability.csv")
  repeat_summary <- risk[, .(
    subjects = .N, any_repeat = sum(MRI1_n_runs > 1L | MRI2_n_runs > 1L),
    max_range_over_5pct = sum(max_same_day_range_pct > 5),
    median_max_range_pct = stats::median(max_same_day_range_pct),
    median_abs_change_in_loss = stats::median(abs(
      hippo_loss_pct_year - median_run_loss_pct_year))), by = event]
  save_csv(repeat_summary, "05_run_reliability_summary.csv")
  qc_queue <- data.table::rbindlist(list(
    risk[, .(OASISID, event, role = "MRI1", day = MRI1_day,
             MUSE_Code = MRI1_Code, n_runs = MRI1_n_runs,
             same_day_range_pct = MRI1_range_pct)],
    risk[, .(OASISID, event, role = "MRI2", day = MRI2_day,
             MUSE_Code = MRI2_Code, n_runs = MRI2_n_runs,
             same_day_range_pct = MRI2_range_pct)]))
  run_lookup <- unique(muse[, .(OASISID, MUSE_Code = Code,
                                selected_run = run)])
  qc_queue <- merge(qc_queue, run_lookup,
                    by = c("OASISID", "MUSE_Code"), all.x = TRUE,
                    sort = FALSE)
  if (nrow(qc_queue) != 2L * nrow(risk))
    stop("The T1 QC manifest did not retain two scans per participant")
  qc_queue[, priority := ifelse(same_day_range_pct > 5, "Review run variation",
                         ifelse(event, "Event case", "Other risk-set scan"))]
  data.table::setorder(qc_queue, priority, OASISID, role)
  save_csv(qc_queue, "06_original_T1_QC_manifest.csv")

  fit_cox <- function(label, x, exposure, horizon_days = Inf, small = FALSE) {
    terms <- if (small)
      c(exposure, "hippo2_pct_ICV", "age_MRI2") else
      c(exposure, "hippo2_pct_ICV", "age_MRI2", "sex",
        "education", "baseline_CDRSB")
    keep <- c(terms, "followup_days", "event")
    d <- as.data.frame(x[stats::complete.cases(as.data.frame(x[, ..keep]))])
    if ("sex" %in% names(d)) d$sex <- factor(d$sex)
    d$time_years <- pmin(d$followup_days, horizon_days) / 365.25
    d$outcome <- d$event & d$followup_days <= horizon_days
    n_event <- sum(d$outcome)
    result <- data.frame(model = label, exposure = exposure,
      subjects = nrow(d), events = n_event, HR = NA_real_,
      CI_low = NA_real_, CI_high = NA_real_, p = NA_real_,
      status = "Fewer than 15 events")
    if (n_event < 15L || nrow(d) < 50L)
      return(list(row = result, fit = NULL, data = d, terms = terms))
    formula <- stats::as.formula(paste(
      "survival::Surv(time_years, outcome) ~", paste(terms, collapse = " + ")))
    model <- tryCatch(survival::coxph(formula, data = d, ties = "efron",
                                      x = TRUE), error = function(e) e)
    if (inherits(model, "error")) {
      result$status <- conditionMessage(model)
      return(list(row = result, fit = NULL, data = d, terms = terms))
    }
    beta <- stats::coef(model)[exposure]
    variance <- stats::vcov(model)[exposure, exposure]
    if (!is.finite(beta) || !is.finite(variance) || variance <= 0) {
      result$status <- "Unstable exposure estimate"
      return(list(row = result, fit = NULL, data = d, terms = terms))
    }
    se <- sqrt(variance)
    result$HR <- exp(beta)
    result$CI_low <- exp(beta - 1.96 * se)
    result$CI_high <- exp(beta + 1.96 * se)
    result$p <- 2 * stats::pnorm(-abs(beta / se))
    result$status <- "Fitted"
    list(row = result, fit = model, data = d, terms = terms)
  }

  cn <- risk[baseline_stage == "CN"]
  primary <- fit_cox("CN, full follow-up, adjusted MRI2 state", cn,
                     "hippo_loss_pct_year")
  five_year <- fit_cox("CN, five-year administrative censoring", cn,
                       "hippo_loss_pct_year", horizon_days = 1826L, small = TRUE)
  recent_dx <- fit_cox("CN, D1 and CDR within 180 days", cn[
    D1_gap_days <= 180L & CDR_gap_days <= 180L], "hippo_loss_pct_year")
  median_run <- fit_cox("CN, median of same-day MUSE runs", cn,
                        "median_run_loss_pct_year")
  models <- list(primary, five_year, recent_dx, median_run)
  estimates <- do.call(rbind, lapply(models, `[[`, "row"))
  save_csv(estimates, "07_event_association_models.csv")
  save_csv(primary$data[, c("OASISID", "MRI1_Code", "MRI2_Code",
                            "hippo_loss_pct_year", "hippo2_pct_ICV",
                            "followup_days", "event")],
           "07b_primary_model_subjects.csv")
  save_csv(data.table::data.table(
    horizon = c("2 years", "5 years", "All observed follow-up"),
    CN_subjects = nrow(cn),
    CN_events = c(sum(cn$event & cn$followup_days <= 730L),
                  sum(cn$event & cn$followup_days <= 1826L),
                  sum(cn$event)),
    CN_with_followup_or_event = c(
      sum(cn$followup_days >= 730L | cn$event & cn$followup_days <= 730L),
      sum(cn$followup_days >= 1826L | cn$event & cn$followup_days <= 1826L),
      sum(cn$followup_days > 0))), "08_event_horizon_coverage.csv")

  if (!is.null(primary$fit)) {
    coefficients <- summary(primary$fit)$coefficients
    save_csv(data.frame(term = rownames(coefficients), coefficients,
                        row.names = NULL), "09_primary_model_all_terms.csv")
    ph <- tryCatch(survival::cox.zph(primary$fit), error = function(e) e)
    if (inherits(ph, "error")) {
      writeLines(conditionMessage(ph), file.path(out, "10_PH_check_error.txt"))
    } else {
      save_csv(data.frame(term = rownames(ph$table), ph$table,
                          row.names = NULL), "10_proportional_hazards_check.csv")
    }
    reduced_terms <- setdiff(primary$terms, "hippo_loss_pct_year")
    reduced_formula <- stats::as.formula(paste(
      "survival::Surv(time_years, outcome) ~",
      paste(reduced_terms, collapse = " + ")))
    reduced <- survival::coxph(reduced_formula, data = primary$data,
                                ties = "efron")
    lr <- 2 * (primary$fit$loglik[2] - reduced$loglik[2])
    save_csv(data.frame(subjects = nrow(primary$data),
      events = sum(primary$data$outcome), likelihood_ratio_chisq = lr,
      p_incremental = stats::pchisq(lr, df = 1L, lower.tail = FALSE)),
      "11_increment_beyond_MRI2_state.csv")
  }

  writeLines(c("STATUS=STEP45_COMPLETE",
               paste0("FROZEN_SUBJECTS=", nrow(risk)),
               paste0("FROZEN_EVENTS=", sum(risk$event)),
               paste0("CN_SUBJECTS=", nrow(cn)),
               paste0("CN_EVENTS=", sum(cn$event)),
               paste0("PRIMARY_MODEL=", primary$row$status),
               paste0("FIVE_YEAR_MODEL=", five_year$row$status)),
             file.path(out, "00_status.txt"))
  message("Step 45 complete: ", out)
  invisible(out)
}

main()
