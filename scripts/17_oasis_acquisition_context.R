local({
  packages <- c("data.table", "survival")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  library(data.table)
  library(survival)

  root <- .study$output
  step40 <- file.path(root, "oasis_event_40/run_frozen")
  inputs <- c(
    risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    crosswalk = file.path(step40, "source_data/Table6_all_312_exact_scan_crosswalk.tsv"),
    prior_models = file.path(step40, "tables/Table7_frozen_model_protocol_sensitivities.tsv"),
    step40_status = file.path(step40, "audit/40_COMPLETED.txt")
  )
  if (any(!file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  if (!"STATUS=STEP40_COMPLETE" %in% readLines(inputs[["step40_status"]], warn = FALSE))
    stop("Step 40 output is incomplete")

  out <- file.path(root, "oasis_event_41",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(x, filename, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, filename), sep = "\t", na = "NA")
  write_tsv(data.table(input = names(inputs), path = unname(inputs),
                       md5 = unname(tools::md5sum(inputs))), "Input_manifest.tsv", "audit")

  message("Reading the frozen MRI2 clinical risk set and exact scanner crosswalk")
  risk <- fread(inputs[["risk"]])[baseline_stage == "CN"]
  risk[, event := as.logical(event)]
  scan <- fread(inputs[["crosswalk"]])
  if (nrow(risk) != 312L || sum(risk$event) != 52L ||
      nrow(scan) != 312L || anyDuplicated(risk$OASISID) ||
      anyDuplicated(scan$OASISID) || !setequal(risk$OASISID, scan$OASISID))
    stop("The frozen 312-person / 52-event sample does not match Step 40")
  if (any(!scan$both_exact | !scan$serial_known | is.na(scan$same_scanner)))
    stop("Scanner identity is unresolved in the Step 40 crosswalk")
  if (any(scan$MRI1_serial == scan$MRI2_serial &
          scan$MRI1_scanner != scan$MRI2_scanner))
    stop("Device serial and model disagree in the crosswalk")

  scan_fields <- c("OASISID", "MRI1_Code", "MRI2_Code", "MRI1_day", "MRI2_day",
                   "MRI1_scanner", "MRI2_scanner", "MRI1_serial", "MRI2_serial",
                   "MRI1_series", "MRI2_series", "MRI1_accession", "MRI2_accession")
  dat <- merge(risk, scan[, ..scan_fields], by = "OASISID", sort = FALSE,
               suffixes = c("", "_scan"))
  if (nrow(dat) != 312L || any(dat$MRI1_Code != dat$MRI1_Code_scan) ||
      any(dat$MRI2_Code != dat$MRI2_Code_scan) ||
      any(dat$MRI1_day != dat$MRI1_day_scan) ||
      any(dat$MRI2_day != dat$MRI2_day_scan))
    stop("The MRI1/MRI2 scans differ between the frozen risk set and Step 40")
  dat[, scanner_changed := as.integer(MRI1_serial != MRI2_serial)]
  dat[, scanner_group := factor(ifelse(scanner_changed == 0L,
                                       "Same physical device", "Changed physical device"),
                                 levels = c("Same physical device", "Changed physical device"))]
  dat[, `:=`(followup_years = followup_days / 365.25,
              MRI_span_years = (MRI2_day - MRI1_day) / 365.25)]
  if (any(dat$MRI_span_years <= 0) || any(dat$followup_years <= 0))
    stop("The MRI2 landmark or follow-up ordering changed")
  if (anyNA(dat[, .(event, age_MRI2, sex, education, baseline_CDRSB,
                    hippo2_pct_ICV, hippo_loss_pct_year, scanner_changed,
                    followup_years, MRI_span_years)]))
    stop("A model input is missing from the frozen risk set")
  dat[, sex := factor(sex)]
  dat[, loss_c := hippo_loss_pct_year - median(hippo_loss_pct_year)]
  loss_reference <- median(dat$hippo_loss_pct_year)
  write_tsv(data.table(participants = nrow(dat), events = sum(dat$event),
                       same_device = sum(dat$scanner_changed == 0L),
                       changed_device = sum(dat$scanner_changed == 1L),
                       loss_center_pct_per_year = loss_reference,
                       time_zero = "MRI2",
                       outcome = "First subsequent all-cause dementia record"),
            "Table1_frozen_landmark_and_exposure_definition.tsv")

  message("Describing device groups and their clinical case mix")
  group_counts <- dat[, .(participants = .N,
                          events = as.integer(sum(event)),
                          observed_event_percent = 100 * mean(event),
                          person_years = sum(followup_years),
                          incidence_per_1000_person_years = 1000 * sum(event) / sum(followup_years),
                          median_followup_years = as.numeric(median(followup_years)),
                          sex_code_2_n = as.integer(sum(sex == levels(sex)[2]))),
                      by = .(scanner_changed, scanner_group)]
  write_tsv(group_counts, "Table2_device_group_event_and_followup_summary.tsv")

  continuous <- c("age_MRI2", "education", "baseline_CDRSB", "hippo2_pct_ICV",
                  "hippo_loss_pct_year", "MRI_span_years", "followup_years")
  case_mix <- rbindlist(lapply(continuous, function(v) {
    x0 <- dat[scanner_changed == 0L][[v]]
    x1 <- dat[scanner_changed == 1L][[v]]
    pooled_sd <- sqrt((sd(x0)^2 + sd(x1)^2) / 2)
    smd <- if (is.finite(pooled_sd) && pooled_sd > 0)
      (mean(x1) - mean(x0)) / pooled_sd else NA_real_
    rbindlist(list(
      data.table(variable = v, group = "Same physical device", n = as.integer(length(x0)),
                 mean = as.numeric(mean(x0)), sd = as.numeric(sd(x0)),
                 median = as.numeric(median(x0)), q1 = as.numeric(quantile(x0, .25)),
                 q3 = as.numeric(quantile(x0, .75)),
                 standardized_mean_difference_changed_minus_same = smd),
      data.table(variable = v, group = "Changed physical device", n = as.integer(length(x1)),
                 mean = as.numeric(mean(x1)), sd = as.numeric(sd(x1)),
                 median = as.numeric(median(x1)), q1 = as.numeric(quantile(x1, .25)),
                 q3 = as.numeric(quantile(x1, .75)),
                 standardized_mean_difference_changed_minus_same = smd)
    ))
  }))
  write_tsv(case_mix, "Table3_device_group_case_mix.tsv")
  scanner_pairs <- dat[, .(participants = .N, events = as.integer(sum(event))),
                       by = .(MRI1_scanner, MRI1_serial, MRI2_scanner,
                              MRI2_serial, scanner_changed)][order(-participants)]
  write_tsv(scanner_pairs, "Table4_exact_device_transition_counts.tsv")

  message("Testing the loss-by-device interaction in the same 312 people")
  model_data <- as.data.frame(dat)
  f_frozen <- Surv(followup_years, event) ~ age_MRI2 + sex + education +
    baseline_CDRSB + hippo2_pct_ICV + loss_c
  f_common <- update(f_frozen, . ~ . + scanner_changed)
  f_interaction <- update(f_common, . ~ . + loss_c:scanner_changed)
  fit_frozen <- coxph(f_frozen, data = model_data, ties = "efron")
  fit_common <- coxph(f_common, data = model_data, ties = "efron")
  fit_interaction <- coxph(f_interaction, data = model_data, ties = "efron")

  prior <- fread(inputs[["prior_models"]])
  ref0 <- prior[model == "Frozen primary CN model"]
  ref1 <- prior[model == "Scanner change adjusted"]
  if (nrow(ref0) != 1L || nrow(ref1) != 1L ||
      abs(exp(coef(fit_frozen)["loss_c"]) - ref0$HR) > 1e-6 ||
      abs(exp(coef(fit_common)["loss_c"]) - ref1$HR) > 1e-6)
    stop("The frozen Step 40 Cox estimates did not reproduce")

  b <- coef(fit_interaction)
  v <- vcov(fit_interaction)
  term <- "loss_c:scanner_changed"
  beta0 <- unname(b["loss_c"])
  beta_int <- unname(b[term])
  beta1 <- beta0 + beta_int
  se0 <- sqrt(v["loss_c", "loss_c"])
  se1 <- sqrt(v["loss_c", "loss_c"] + v[term, term] +
                2 * v["loss_c", term])
  se_int <- sqrt(v[term, term])
  lr <- 2 * (fit_interaction$loglik[2] - fit_common$loglik[2])
  interaction <- data.table(
    contrast = c("Same-device loss association", "Changed-device loss association",
                 "Changed-versus-same association ratio"),
    log_HR = c(beta0, beta1, beta_int), se = c(se0, se1, se_int),
    HR = exp(c(beta0, beta1, beta_int)),
    CI_low = exp(c(beta0, beta1, beta_int) - 1.96 * c(se0, se1, se_int)),
    CI_high = exp(c(beta0, beta1, beta_int) + 1.96 * c(se0, se1, se_int)),
    Wald_p = 2 * pnorm(-abs(c(beta0, beta1, beta_int) /
                              c(se0, se1, se_int))),
    participants = as.integer(c(group_counts$participants[group_counts$scanner_changed == 0L],
                                group_counts$participants[group_counts$scanner_changed == 1L], nrow(dat))),
    events = as.integer(c(group_counts$events[group_counts$scanner_changed == 0L],
                          group_counts$events[group_counts$scanner_changed == 1L], sum(dat$event))))
  write_tsv(interaction, "Table5_direct_interaction_and_group_effects.tsv")
  model_tests <- data.table(
    test = "Loss-by-device interaction, 1 df likelihood-ratio test",
    participants = nrow(dat), events = sum(dat$event),
    LR_chisq = lr, LR_p = pchisq(lr, df = 1, lower.tail = FALSE),
    Wald_interaction_p = interaction$Wald_p[3],
    logLik_common_slope = fit_common$loglik[2],
    logLik_interaction = fit_interaction$loglik[2])
  write_tsv(model_tests, "Table6_nested_Cox_interaction_test.tsv")
  write_tsv(data.table(term = rownames(cox.zph(fit_interaction)$table),
                       cox.zph(fit_interaction)$table),
            "Table7_proportional_hazards_check.tsv")
  write_tsv(dat[, .(OASISID, scanner_group, scanner_changed, event,
                    followup_years, age_MRI2, sex, education, baseline_CDRSB,
                    hippo2_pct_ICV, hippo_loss_pct_year, MRI_span_years,
                    MRI1_Code, MRI2_Code, MRI1_scanner, MRI2_scanner,
                    MRI1_serial, MRI2_serial, MRI1_accession, MRI2_accession)],
            "Table8_participant_level_source_data.tsv", "source_data")

  writeLines(capture.output(sessionInfo()), file.path(out, "audit/R_sessionInfo.txt"))
  writeLines("STATUS=STEP41_COMPLETE", file.path(out, "audit/41_COMPLETED.txt"))
  message("Step 41 complete: ", out)
  print(group_counts)
  print(interaction)
  print(model_tests)
})
