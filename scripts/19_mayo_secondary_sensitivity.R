local({
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Install the survival package before running this script.")

  root <- .study$output
  step43 <- file.path(root, "mcsa_event_43/run_frozen")
  risk_file <- file.path(step43, "tables/Table2_person_level_frozen_risk_sets.tsv")
  result_file <- file.path(step43, "tables/Table4_MCSA_event_associations.tsv")
  oasis_file <- file.path(root, paste0(
    "oasis_event_38/run_frozen/tables/",
    "Table4_primary_increment_and_frozen_reproduction.tsv"))
  zip_file <- file.path(.study$data,"mayo","Demographics_and_Assessments.zip")
  out <- file.path(root, "mcsa_event_44",
                   paste0("run_", "frozen"))
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, filename) {
    write.table(x, file.path(out, "tables", filename), sep = "\t",
                row.names = FALSE, quote = FALSE, na = "")
  }

  message("Reading frozen MCSA CN risk set")
  risk <- read.delim(risk_file, stringsAsFactors = FALSE)
  if ("same_scanner" %in% names(risk))
    names(risk)[names(risk) == "same_scanner"] <- "same_manufacturer"
  cn <- risk[risk$stage == "CN", , drop = FALSE]
  cn$MRI2_date <- as.Date(cn$MRI2_date)
  if (nrow(cn) != 925L || sum(cn$event) != 21L || anyDuplicated(cn$id))
    stop("The frozen Step 43 CN risk set differs from 925 participants / 21 events.")

  raw <- read.csv(unz(zip_file, "MCSA_Data_12Nov2026.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  raw$VisitDate <- as.Date(raw$VisitDate)
  raw$FCOGIMPR <- as.numeric(raw$FCOGIMPR)
  raw$CDRSUM <- as.numeric(raw$CDRSUM)
  visits_by_person <- split(raw, raw$MCSA_ID, drop = TRUE)

  message("Auditing first dementia records and later diagnosis reversions")
  event_rows <- cn[cn$event == 1L, , drop = FALSE]
  event_audit <- do.call(rbind, lapply(seq_len(nrow(event_rows)), function(i) {
    r <- event_rows[i, , drop = FALSE]
    visits <- visits_by_person[[r$id]]
    visits <- visits[!is.na(visits$VisitDate) &
                       visits$VisitDate > r$MRI2_date, , drop = FALSE]
    unknown <- visits$VisitDate[visits$FCOGIMPR %in% 4]
    if (length(unknown))
      visits <- visits[visits$VisitDate < min(unknown), , drop = FALSE]
    visits <- visits[visits$FCOGIMPR %in% c(0, 1, 3), , drop = FALSE]
    visits <- visits[order(visits$VisitDate), , drop = FALSE]
    first <- which(visits$FCOGIMPR == 3)[1L]
    if (is.na(first)) stop("Cannot reconstruct frozen event for ", r$id)
    event_date <- visits$VisitDate[first]
    prior_date <- max(c(r$MRI2_date, visits$VisitDate[visits$VisitDate <
                                                      event_date]))
    data.frame(id = r$id, MRI2_date = r$MRI2_date,
               first_dementia_date = event_date,
               interval_since_last_known_visit_days =
                 as.integer(event_date - prior_date),
               later_reversion = as.integer(any(visits$VisitDate > event_date &
                                                  visits$FCOGIMPR %in% c(0, 1))),
               CDRSB_at_first_dementia = visits$CDRSUM[first],
               hippo_loss_pct_year = r$hippo_loss_pct_year,
               same_manufacturer = r$same_manufacturer)
  }))
  if (sum(event_audit$later_reversion) != 3L ||
      sum(event_audit$interval_since_last_known_visit_days > 730L) != 2L)
    stop("The event-definition audit differs from the frozen 3 reversions / 2 long intervals.")
  save_tsv(event_audit, "Table1_all_21_first_dementia_event_audit.tsv")

  reduced <- survival::Surv(followup_years, event) ~ age_MRI2 + male +
    education + baseline_CDRSB + hippo2_pct_ICV
  full <- survival::Surv(followup_years, event) ~ age_MRI2 + male +
    education + baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year
  fit_models <- function(d, label) {
    base_fit <- survival::coxph(reduced, data = d, ties = "efron")
    full_fit <- survival::coxph(full, data = d, ties = "efron", x = TRUE)
    term <- "hippo_loss_pct_year"
    b <- unname(stats::coef(full_fit)[term])
    se <- sqrt(stats::vcov(full_fit)[term, term])
    lr <- 2 * (as.numeric(logLik(full_fit)) - as.numeric(logLik(base_fit)))
    row <- data.frame(model = label, participants = nrow(d),
                      events = sum(d$event), HR_per_1pct_year = exp(b),
                      CI_low = exp(b - stats::qnorm(0.975) * se),
                      CI_high = exp(b + stats::qnorm(0.975) * se),
                      Wald_p = 2 * stats::pnorm(-abs(b / se)),
                      LR_p = stats::pchisq(lr, df = 1, lower.tail = FALSE),
                      stringsAsFactors = FALSE)
    list(row = row, full = full_fit)
  }

  original <- read.delim(result_file, stringsAsFactors = FALSE)
  original <- original[original$stage == "CN", , drop = FALSE]
  primary <- fit_models(cn, "MCSA CN primary")
  if (abs(primary$row$HR_per_1pct_year - original$HR_per_1pct_year) > 1e-7)
    stop("The fixed primary model does not reproduce the Step 43 HR.")

  reversal_ids <- event_audit$id[event_audit$later_reversion == 1L]
  gap_ids <- event_audit$id[
    event_audit$interval_since_last_known_visit_days > 730L]
  loss_limits <- stats::quantile(cn$hippo_loss_pct_year, c(0.01, 0.99))
  sets <- list(
    "MCSA CN primary" = cn,
    "Exclude later diagnostic reversions" = cn[!cn$id %in% reversal_ids, ],
    "Exclude events with >2y diagnosis interval" = cn[!cn$id %in% gap_ids, ],
    "Exclude either event-definition concern" =
      cn[!cn$id %in% union(reversal_ids, gap_ids), ],
    "Same MRI manufacturer" = cn[cn$same_manufacturer == 1L, ],
    "Trim exposure outside 1st-99th percentiles" =
      cn[cn$hippo_loss_pct_year >= loss_limits[1L] &
           cn$hippo_loss_pct_year <= loss_limits[2L], ]
  )
  message("Fitting fixed-model event and measurement sensitivities")
  sensitivity <- do.call(rbind, lapply(names(sets), function(label)
    fit_models(sets[[label]], label)$row))
  rownames(sensitivity) <- NULL
  save_tsv(sensitivity, "Table2_fixed_model_sensitivities.tsv")

  message("Calculating one-person-out influence for all 925 participants")
  ref_log_hr <- log(primary$row$HR_per_1pct_year)
  loo <- do.call(rbind, lapply(seq_len(nrow(cn)), function(i) {
    observed_warnings <- character()
    fit <- tryCatch(withCallingHandlers(
      survival::coxph(full, data = cn[-i, , drop = FALSE], ties = "efron"),
      warning = function(w) {
        observed_warnings <<- c(observed_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(fit, "error"))
      return(data.frame(id = cn$id[i], removed_event = cn$event[i],
                        HR = NA_real_, Wald_p = NA_real_,
                        delta_log_HR = NA_real_, status = conditionMessage(fit)))
    b <- unname(stats::coef(fit)["hippo_loss_pct_year"])
    se <- sqrt(stats::vcov(fit)["hippo_loss_pct_year",
                                "hippo_loss_pct_year"])
    data.frame(id = cn$id[i], removed_event = cn$event[i],
               HR = exp(b), Wald_p = 2 * stats::pnorm(-abs(b / se)),
               delta_log_HR = b - ref_log_hr,
               status = if (length(observed_warnings))
                 paste(unique(observed_warnings), collapse = " | ") else "OK")
  }))
  loo$status <- gsub("[\r\n\t]", " ", loo$status)
  loo$abs_delta_log_HR <- abs(loo$delta_log_HR)
  save_tsv(loo, "Table3_all_925_one_person_out_estimates.tsv")
  save_tsv(head(loo[order(-loo$abs_delta_log_HR), ], 15L),
           "Table4_largest_person_level_influences.tsv")
  loo_summary <- data.frame(
    attempts = nrow(loo), successful = sum(is.finite(loo$HR)),
    min_HR = min(loo$HR, na.rm = TRUE),
    max_HR = max(loo$HR, na.rm = TRUE),
    max_absolute_log_HR_change = max(loo$abs_delta_log_HR, na.rm = TRUE),
    all_successful_HR_above_1 = all(loo$HR[is.finite(loo$HR)] > 1),
    all_successful_Wald_p_below_0_05 =
      all(loo$Wald_p[is.finite(loo$Wald_p)] < 0.05))
  save_tsv(loo_summary, "Table5_one_person_out_summary.tsv")

  ph <- survival::cox.zph(primary$full, transform = "km")
  ph_table <- data.frame(term = rownames(ph$table), ph$table,
                         row.names = NULL, check.names = FALSE)
  save_tsv(ph_table, "Table6_proportional_hazards_check.tsv")

  oasis <- read.delim(oasis_file, stringsAsFactors = FALSE)
  forest <- rbind(
    data.frame(model = "OASIS-3 CN primary", participants = oasis$participants[1L],
               events = oasis$events[1L],
               HR_per_1pct_year = oasis$HR_per_1pct_year[1L],
               CI_low = oasis$CI_low[1L], CI_high = oasis$CI_high[1L]),
    sensitivity[c("model", "participants", "events", "HR_per_1pct_year",
                  "CI_low", "CI_high")])
  save_tsv(forest, "Table7_OASIS_and_MCSA_forest_source.tsv")
  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  writeLines("Step 44 complete", file.path(out, "44_COMPLETED.txt"))
  message("Analysis complete: ", out)
  print(sensitivity, row.names = FALSE)
  print(loo_summary, row.names = FALSE)
  print(ph_table, row.names = FALSE)
})
