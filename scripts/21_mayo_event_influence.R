local({
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Install the survival package before running this script.")

  root <- .study$output
  risk_path <- file.path(root, paste0(
    "mcsa_event_43/run_frozen/tables/",
    "Table2_person_level_frozen_risk_sets.tsv"))
  reference_path <- file.path(root, paste0(
    "mcsa_event_45/run_frozen/tables/",
    "Table2_time_varying_coefficient_test.tsv"))
  influence_path <- file.path(root,
    "mcsa_event_44/run_frozen/tables/Table3_all_925_one_person_out_estimates.tsv")
  influence <- read.delim(influence_path, stringsAsFactors = FALSE)
  eligible <- influence[is.finite(influence$abs_delta_log_HR), , drop = FALSE]
  if (!nrow(eligible)) stop("No finite leave-one-participant-out influence estimates.")
  case_id <- eligible$id[which.max(eligible$abs_delta_log_HR)]
  out <- file.path(root, "mcsa_event_46",
                   paste0("run_", "frozen"))
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name) {
    write.table(x, file.path(out, "tables", name), sep = "\t",
                row.names = FALSE, quote = FALSE, na = "")
  }

  risk <- read.delim(risk_path, stringsAsFactors = FALSE)
  cn <- risk[risk$stage == "CN", , drop = FALSE]
  if (nrow(cn) != 925L || sum(cn$event) != 21L ||
      sum(cn$id == case_id) != 1L || cn$event[cn$id == case_id] != 1L)
    stop("The frozen CN risk set or the targeted event differs from Step 43.")

  formula_ph <- survival::Surv(followup_years, event) ~ age_MRI2 + male +
    education + baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year
  formula_time <- survival::Surv(followup_years, event) ~ age_MRI2 + male +
    education + baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year +
    tt(hippo_loss_pct_year)
  fit_one <- function(d, label) {
    ph <- survival::coxph(formula_ph, data = d, ties = "efron")
    tv <- survival::coxph(formula_time, data = d, ties = "efron",
                          tt = function(x, t, ...) x * log(t / 5))
    terms <- c("hippo_loss_pct_year", "tt(hippo_loss_pct_year)")
    b <- stats::coef(tv)[terms]
    v <- stats::vcov(tv)[terms, terms, drop = FALSE]
    se <- sqrt(diag(v))
    ph_b <- unname(stats::coef(ph)["hippo_loss_pct_year"])
    ph_se <- sqrt(stats::vcov(ph)["hippo_loss_pct_year",
                                   "hippo_loss_pct_year"])
    lr <- 2 * (as.numeric(logLik(tv)) - as.numeric(logLik(ph)))
    model <- data.frame(
      group = label, participants = nrow(d), events = sum(d$event),
      PH_HR = exp(ph_b),
      PH_CI_low = exp(ph_b - stats::qnorm(0.975) * ph_se),
      PH_CI_high = exp(ph_b + stats::qnorm(0.975) * ph_se),
      log_HR_at_5y = unname(b[1L]),
      log_time_coefficient = unname(b[2L]),
      log_time_coefficient_SE = unname(se[2L]),
      log_time_Wald_p = 2 * stats::pnorm(-abs(b[2L] / se[2L])),
      time_term_LR_chisq = lr,
      time_term_LR_p = stats::pchisq(lr, 1, lower.tail = FALSE))
    estimate <- function(years) {
      scale <- log(years / 5)
      log_hr <- unname(b[1L] + b[2L] * scale)
      var_hr <- v[1L, 1L] + 2 * scale * v[1L, 2L] +
        scale^2 * v[2L, 2L]
      width <- stats::qnorm(0.975) * sqrt(pmax(var_hr, 0))
      data.frame(group = label, years_after_MRI2 = years,
                 HR_per_1pct_year = exp(log_hr),
                 CI_low = exp(log_hr - width),
                 CI_high = exp(log_hr + width),
                 participants_at_risk = vapply(years, function(t)
                   sum(d$followup_years >= t), integer(1)),
                 events_observed_by_time = vapply(years, function(t)
                   sum(d$event == 1L & d$followup_years <= t), integer(1)))
    }
    list(model = model, times = estimate(c(2, 3, 5, 8, 10)),
         curve = estimate(seq(2, 10, length.out = 100)))
  }

  message("Fitting the frozen time model with and without the targeted event")
  full <- fit_one(cn, "All frozen CN participants")
  minus_case <- fit_one(cn[cn$id != case_id, , drop = FALSE],
                        "Influence sensitivity")
  ref <- read.delim(reference_path, stringsAsFactors = FALSE)
  if (abs(full$model$time_term_LR_p - ref$time_term_LR_p) > 1e-7 ||
      abs(full$model$log_time_coefficient - ref$log_time_coefficient) > 1e-7)
    stop("The full-sample time model does not reproduce Step 45.")

  models <- rbind(full$model, minus_case$model)
  times <- rbind(full$times, minus_case$times)
  curve <- rbind(full$curve, minus_case$curve)
  save_tsv(models, "Table1_fixed_time_model_influence.tsv")
  save_tsv(times, "Table2_time_specific_HR_with_and_without_case.tsv")
  save_tsv(curve, "Figure46_time_curve_source.tsv")

  writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
  writeLines("Step 46 complete", file.path(out, "46_COMPLETED.txt"))
  message("Analysis complete: ", out)
  print(models, row.names = FALSE)
  print(times, row.names = FALSE)
})
