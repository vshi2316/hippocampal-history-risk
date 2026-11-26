# Check event definitions in the OASIS-3 cognitively normal risk set.

needed <- c("data.table", "survival")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))

main <- function() {
  code_dir <- .study$output
  input <- file.path(code_dir, "45_OASIS3_frozen_event_link_results")
  out <- file.path(code_dir, "46_OASIS3_event_definition_final_gate_results")
  files <- c("00_status.txt", "01_frozen_event_risk_set.csv",
             "02_event_diagnosis_audit.csv", "07_event_association_models.csv")
  paths <- file.path(input, files)
  if (any(!file.exists(paths)))
    stop("Missing Step 45 output: ", paste(files[!file.exists(paths)], collapse = ", "))
  if (!any(readLines(paths[1], warn = FALSE) == "STATUS=STEP45_COMPLETE"))
    stop("Step 45 is incomplete")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines("STATUS=STEP46_INCOMPLETE", file.path(out, "00_status.txt"))
  save_csv <- function(x, name)
    data.table::fwrite(x, file.path(out, name), na = "")

  risk <- data.table::fread(paths[2])
  audit <- data.table::fread(paths[3])
  old <- data.table::fread(paths[4])
  cn <- risk[baseline_stage == "CN"]
  if (nrow(risk) != 317L || sum(risk$event) != 53L ||
      nrow(cn) != 312L || sum(cn$event) != 52L ||
      anyDuplicated(risk$OASISID) || anyDuplicated(audit$OASISID))
    stop("The frozen Step 45 risk set or diagnosis audit changed")
  audit <- audit[, .(OASISID, prior_D1_day, prior_D1_stage,
                     diagnosis_interval_days, next_D1_stage,
                     next_D1_within_18m, CDRSB_near_event)]
  cn <- merge(cn, audit, by = "OASISID", all.x = TRUE, sort = FALSE)
  if (nrow(cn) != 312L || sum(cn$event & !is.na(cn$prior_D1_day)) != 52L)
    stop("Diagnosis visits did not map to all 52 CN events")
  cn[, event := as.logical(event)]
  cn[, reversion := event & next_D1_within_18m %in% TRUE &
       next_D1_stage %in% c("CN", "MCI")]
  if (sum(cn$reversion) != 7L)
    stop("The seven Step 45 diagnostic reversions changed")
  cn[, `:=`(analysis_event = event,
            analysis_days = as.numeric(followup_days))]
  save_csv(cn[, .(OASISID, event, reversion, baseline_CDRSB,
                  hippo_loss_pct_year, hippo2_pct_ICV, MRI2_day,
                  first_dementia_day, prior_D1_day,
                  diagnosis_interval_days, next_D1_stage,
                  CDRSB_near_event, followup_days)],
           "01_frozen_CN_event_audit.csv")

  fit_cox <- function(label, x, include_cdr = TRUE) {
    terms <- c("hippo_loss_pct_year", "hippo2_pct_ICV", "age_MRI2",
               "sex", "education")
    if (include_cdr) terms <- c(terms, "baseline_CDRSB")
    keep <- c(terms, "analysis_days", "analysis_event")
    complete <- stats::complete.cases(as.data.frame(x[, ..keep]))
    d <- as.data.frame(x[complete])
    d$sex <- factor(d$sex)
    d$analysis_years <- d$analysis_days / 365.25
    row <- data.frame(model = label, subjects = nrow(d),
                      events = sum(d$analysis_event), HR = NA_real_,
                      CI_low = NA_real_, CI_high = NA_real_,
                      p = NA_real_, status = "Fewer than 20 events")
    if (row$events < 20L || nrow(d) < 50L)
      return(list(row = row, data = d))
    formula <- stats::as.formula(paste(
      "survival::Surv(analysis_years, analysis_event) ~",
      paste(terms, collapse = " + ")))
    fit <- tryCatch(survival::coxph(formula, data = d, ties = "efron"),
                    error = function(e) e)
    if (inherits(fit, "error")) {
      row$status <- conditionMessage(fit)
      return(list(row = row, data = d))
    }
    term <- "hippo_loss_pct_year"
    beta <- stats::coef(fit)[term]
    variance <- stats::vcov(fit)[term, term]
    if (!is.finite(beta) || !is.finite(variance) || variance <= 0) {
      row$status <- "Unstable exposure estimate"
      return(list(row = row, data = d))
    }
    se <- sqrt(variance)
    row$HR <- exp(beta)
    row$CI_low <- exp(beta - 1.96 * se)
    row$CI_high <- exp(beta + 1.96 * se)
    row$p <- 2 * stats::pnorm(-abs(beta / se))
    row$status <- "Fitted"
    list(row = row, data = d)
  }

  reference <- fit_cox("Step 45 reference", cn)
  old_ref <- old[model == "CN, full follow-up, adjusted MRI2 state"]
  if (nrow(old_ref) != 1L || reference$row$status != "Fitted" ||
      abs(reference$row$HR - old_ref$HR) > 1e-6)
    stop("The Step 45 Cox reference did not reproduce")

  zero_cdr <- fit_cox("Baseline CDR-SB = 0", cn[baseline_CDRSB == 0],
                      include_cdr = FALSE)
  nonreverting <- data.table::copy(cn)
  nonreverting[reversion == TRUE, analysis_event := FALSE]
  nonreversion_fit <- fit_cox("Reversions censored at first dementia visit",
                             nonreverting)

  midpoint <- data.table::copy(cn)
  midpoint[event == TRUE, interval_start_day := pmax(
    MRI2_day, prior_D1_day)]
  midpoint[event == TRUE, analysis_days := pmax(1,
    (interval_start_day + first_dementia_day) / 2 - MRI2_day)]
  midpoint_fit <- fit_cox("Event time at diagnosis-interval midpoint", midpoint)

  joint <- data.table::copy(midpoint[baseline_CDRSB == 0])
  joint[reversion == TRUE, `:=`(analysis_event = FALSE,
                                analysis_days = as.numeric(followup_days))]
  joint_fit <- fit_cox("Combined diagnosis and baseline check", joint,
                       include_cdr = FALSE)

  results <- do.call(rbind, lapply(list(reference, zero_cdr,
                                         nonreversion_fit, midpoint_fit,
                                         joint_fit), `[[`, "row"))
  save_csv(results, "02_event_definition_sensitivity_models.csv")
  save_csv(data.table::data.table(
    metric = c("Frozen CN subjects", "Frozen first dementia records",
               "Reversion records", "Events with diagnosis gap over 2 years",
               "Baseline CDR-SB = 0 subjects",
               "Baseline CDR-SB = 0 events",
               "Other baseline CDR-SB subjects",
               "Other baseline CDR-SB events"),
    count = c(nrow(cn), sum(cn$event), sum(cn$reversion),
              sum(cn$event & cn$diagnosis_interval_days > 730L,
                  na.rm = TRUE),
              sum(cn$baseline_CDRSB == 0),
              sum(cn$event & cn$baseline_CDRSB == 0),
              sum(cn$baseline_CDRSB != 0),
              sum(cn$event & cn$baseline_CDRSB != 0))),
           "03_event_definition_counts.csv")
  timing <- cn[event == TRUE, .(
    OASISID, years_after_MRI2 = followup_days / 365.25,
    diagnosis_interval_years = diagnosis_interval_days / 365.25,
    midpoint_years_after_MRI2 = midpoint[event == TRUE]$analysis_days / 365.25,
    reversion, baseline_CDRSB, CDRSB_near_event)]
  save_csv(timing, "04_diagnosis_timing_figure_data.csv")

  decision <- data.frame(
    check = results$model,
    subjects = results$subjects, events = results$events,
    HR_above_one = ifelse(results$status == "Fitted", results$HR > 1, NA),
    CI_excludes_one = ifelse(results$status == "Fitted",
                             results$CI_low > 1, NA))
  save_csv(decision, "05_evidence_consistency_gate.csv")

  writeLines(c("STATUS=STEP46_COMPLETE",
               paste0("CN_SUBJECTS=", nrow(cn)),
               paste0("FIRST_DEMENTIA_RECORDS=", sum(cn$event)),
               paste0("REVERSIONS=", sum(cn$reversion)),
               paste0("FITTED_MODELS=", sum(results$status == "Fitted"))),
             file.path(out, "00_status.txt"))
  message("Step 46 complete: ", out)
  invisible(out)
}

main()
