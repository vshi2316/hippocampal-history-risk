local({
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Install the survival package before running this script.")

  root <- .study$output
  risk_file <- file.path(root, paste0(
    "mcsa_event_43/run_frozen/tables/",
    "Table2_person_level_frozen_risk_sets.tsv"))
  step43_file <- file.path(root, paste0(
    "mcsa_event_43/run_frozen/tables/",
    "Table4_MCSA_event_associations.tsv"))
  step44_file <- file.path(root, paste0(
    "mcsa_event_44/run_frozen/tables/",
    "Table6_proportional_hazards_check.tsv"))
  out <- file.path(root, "mcsa_event_45",
                   paste0("run_", "frozen"))
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name) {
    write.table(x, file.path(out, "tables", name), sep = "\t",
                row.names = FALSE, quote = FALSE, na = "")
  }

  message("Reading the frozen MCSA CN risk set")
  risk <- read.delim(risk_file, stringsAsFactors = FALSE)
  cn <- risk[risk$stage == "CN", , drop = FALSE]
  if (nrow(cn) != 925L || sum(cn$event) != 21L || anyDuplicated(cn$id))
    stop("The Step 43 CN risk set differs from 925 participants / 21 events.")

  frozen <- read.delim(step43_file, stringsAsFactors = FALSE)
  frozen <- frozen[frozen$stage == "CN", , drop = FALSE]
  formula_ph <- survival::Surv(followup_years, event) ~ age_MRI2 + male +
    education + baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year
  ph_fit <- survival::coxph(formula_ph, data = cn, ties = "efron", x = TRUE)
  ph_hr <- exp(unname(stats::coef(ph_fit)["hippo_loss_pct_year"]))
  if (abs(ph_hr - frozen$HR_per_1pct_year) > 1e-7)
    stop("The proportional-hazards model does not reproduce Step 43.")

  message("Checking Schoenfeld residuals against follow-up time")
  z_km <- survival::cox.zph(ph_fit, transform = "km")
  z_time <- survival::cox.zph(ph_fit, transform = "identity")
  saved_ph <- read.delim(step44_file, stringsAsFactors = FALSE)
  saved_p <- saved_ph$p[saved_ph$term == "hippo_loss_pct_year"]
  if (abs(z_km$table["hippo_loss_pct_year", "p"] - saved_p) > 1e-7)
    stop("The KM-transform PH test does not reproduce Step 44.")
  ph_table <- rbind(
    data.frame(transform = "KM", term = rownames(z_km$table),
               z_km$table, row.names = NULL, check.names = FALSE),
    data.frame(transform = "Follow-up years", term = rownames(z_time$table),
               z_time$table, row.names = NULL, check.names = FALSE))
  save_tsv(ph_table, "Table1_proportional_hazards_time_diagnostics.tsv")

  message("Estimating one time-varying coefficient for prior loss")
  formula_time <- survival::Surv(followup_years, event) ~ age_MRI2 + male +
    education + baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year +
    tt(hippo_loss_pct_year)
  time_fit <- survival::coxph(
    formula_time, data = cn, ties = "efron", x = TRUE,
    tt = function(x, t, ...) x * log(t / 5))
  terms <- c("hippo_loss_pct_year", "tt(hippo_loss_pct_year)")
  b <- stats::coef(time_fit)[terms]
  v <- stats::vcov(time_fit)[terms, terms, drop = FALSE]
  if (any(!is.finite(b)) || any(!is.finite(v)))
    stop("The time-varying coefficient is not estimable.")
  se <- sqrt(diag(v))
  lr <- 2 * (as.numeric(logLik(time_fit)) - as.numeric(logLik(ph_fit)))
  model_test <- data.frame(participants = nrow(cn), events = sum(cn$event),
                           primary_PH_HR = ph_hr,
                           log_HR_at_5y = unname(b[1L]),
                           log_time_coefficient = unname(b[2L]),
                           log_time_coefficient_SE = unname(se[2L]),
                           log_time_Wald_p = 2 * stats::pnorm(-abs(b[2L] / se[2L])),
                           time_term_LR_chisq = lr,
                           time_term_LR_p = stats::pchisq(lr, 1, lower.tail = FALSE))
  save_tsv(model_test, "Table2_time_varying_coefficient_test.tsv")

  estimate_at <- function(years) {
    scale <- log(years / 5)
    log_hr <- unname(b[1L] + b[2L] * scale)
    var_hr <- v[1L, 1L] + 2 * scale * v[1L, 2L] +
      scale^2 * v[2L, 2L]
    width <- stats::qnorm(0.975) * sqrt(pmax(var_hr, 0))
    data.frame(years_after_MRI2 = years, HR_per_1pct_year = exp(log_hr),
               CI_low = exp(log_hr - width),
               CI_high = exp(log_hr + width),
               participants_at_risk = vapply(years, function(t)
                 sum(cn$followup_years >= t), integer(1)),
               events_observed_by_time = vapply(years, function(t)
                 sum(cn$event == 1L & cn$followup_years <= t), integer(1)))
  }
  times <- estimate_at(c(2, 3, 5, 8, 10))
  save_tsv(times, "Table3_time_specific_loss_associations.tsv")

  bins <- data.frame(start = c(0, 2, 4, 6, 8), end = c(2, 4, 6, 8, 12))
  bins$participants_at_start <- vapply(bins$start, function(t)
    sum(cn$followup_years > t), integer(1))
  bins$events_in_interval <- vapply(seq_len(nrow(bins)), function(i)
    sum(cn$event == 1L & cn$followup_years > bins$start[i] &
          cn$followup_years <= bins$end[i]), integer(1))
  save_tsv(bins, "Table4_event_timing_and_risk_sets.tsv")

  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  writeLines("Step 45 complete", file.path(out, "45_COMPLETED.txt"))
  message("Analysis complete: ", out)
  print(model_test, row.names = FALSE)
  print(times, row.names = FALSE)
  print(ph_table[ph_table$term %in% c("hippo_loss_pct_year", "GLOBAL"), ],
        row.names = FALSE)
})
