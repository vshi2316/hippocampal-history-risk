local({
  packages <- c("data.table", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))

  library(data.table)
  library(survival)
  input <- file.path(.study$output,"45_OASIS3_frozen_event_link_results")
  sensitivity_file <- file.path(.study$output,"46_OASIS3_event_definition_final_gate_results/02_event_definition_sensitivity_models.csv")
  files <- c(
    status = file.path(input, "00_status.txt"),
    risk = file.path(input, "01_frozen_event_risk_set.csv"),
    subjects = file.path(input, "07b_primary_model_subjects.csv"),
    original_model = file.path(input, "07_event_association_models.csv"),
    original_increment = file.path(input, "11_increment_beyond_MRI2_state.csv"),
    sensitivity = sensitivity_file
  )
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) stop("Missing frozen input: ", paste(missing_files, collapse = ", "))
  if (!"STATUS=STEP45_COMPLETE" %in% readLines(files[["status"]], warn = FALSE))
    stop("Step 45 output is incomplete")

  out <- file.path(file.path(.study$output,"oasis_event_38"),
                   paste0("run_", "frozen"))
  for (subdir in c("tables", "source_data", "audit"))
    dir.create(file.path(out, subdir), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, subdir = "tables")
    fwrite(as.data.table(x), file.path(out, subdir, name), sep = "\t", na = "NA")

  save_tsv(data.table(input = names(files), file = unname(files),
                      md5 = unname(tools::md5sum(files))), "Input_manifest.tsv", "audit")

  risk <- fread(files[["risk"]])
  original <- fread(files[["original_model"]])
  original_increment <- fread(files[["original_increment"]])
  original_subjects <- fread(files[["subjects"]])
  sensitivity <- fread(files[["sensitivity"]])
  required <- c("OASISID", "baseline_stage", "MRI1_day", "MRI2_day",
                "hippo_loss_pct_year", "hippo2_pct_ICV", "age_MRI2", "sex",
                "education", "baseline_CDRSB", "followup_days", "event")
  if (!all(required %in% names(risk)))
    stop("Frozen risk set is missing: ", paste(setdiff(required, names(risk)), collapse = ", "))
  if (anyDuplicated(risk$OASISID)) stop("Duplicate OASISID in frozen risk set")
  risk[, event := as.logical(event)]
  if (nrow(risk) != 317L || sum(risk$event) != 53L)
    stop("Frozen full risk set differs from 317 participants / 53 events")

  cn <- copy(risk[baseline_stage == "CN"])
  if (nrow(cn) != 312L || sum(cn$event) != 52L ||
      sum(cn$baseline_CDRSB == 0.5) != 4L)
    stop("Frozen CN risk set differs from 312 participants / 52 events / four CDR-SB 0.5")
  if (any(cn$MRI1_day >= cn$MRI2_day) ||
      any(cn$followup_days <= 0) ||
      any(cn$event & cn$first_dementia_day <= cn$MRI2_day))
    stop("MRI and diagnosis time ordering changed")
  save_tsv(data.table(stage = c("All frozen participants", "CN at MRI2"),
                      participants = c(nrow(risk), nrow(cn)),
                      first_dementia_records = c(sum(risk$event), sum(cn$event))),
           "Table1_frozen_risk_set_flow.tsv")
  save_tsv(cn[, .N, by = .(baseline_CDRSB)][order(baseline_CDRSB)],
           "Table2_CN_baseline_CDRSB.tsv")

  model_fields <- c("OASISID", "followup_days", "event", "age_MRI2",
                    "sex", "education", "baseline_CDRSB",
                    "hippo2_pct_ICV", "hippo_loss_pct_year")
  analysis <- as.data.frame(cn[complete.cases(cn[, ..model_fields])])
  analysis$sex <- factor(analysis$sex)
  analysis$followup_years <- analysis$followup_days / 365.25
  if (nrow(analysis) != 312L || sum(analysis$event) != 52L ||
      !setequal(analysis$OASISID, original_subjects$OASISID))
    stop("Complete-case model set differs from the frozen Step 45 subject set")

  m0 <- coxph(Surv(followup_years, event) ~ age_MRI2 + sex + education +
                baseline_CDRSB, data = analysis, ties = "efron")
  m1 <- coxph(Surv(followup_years, event) ~ age_MRI2 + sex + education +
                baseline_CDRSB + hippo2_pct_ICV, data = analysis, ties = "efron")
  m2 <- coxph(Surv(followup_years, event) ~ age_MRI2 + sex + education +
                baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year,
              data = analysis, ties = "efron")

  term <- "hippo_loss_pct_year"
  beta <- unname(coef(m2)[term])
  se <- sqrt(vcov(m2)[term, term])
  hr <- exp(beta)
  ci <- exp(beta + c(-1, 1) * 1.96 * se)
  wald_p <- 2 * pnorm(-abs(beta / se))
  lr <- 2 * (m2$loglik[2] - m1$loglik[2])
  lr_p <- pchisq(lr, df = 1, lower.tail = FALSE)
  ref <- original[model == "CN, full follow-up, adjusted MRI2 state"]
  if (nrow(ref) != 1L || nrow(original_increment) != 1L ||
      abs(hr - ref$HR) > 1e-6 ||
      max(abs(ci - c(ref$CI_low, ref$CI_high))) > 1e-6 ||
      abs(wald_p - ref$p) > 1e-6 ||
      abs(lr - original_increment$likelihood_ratio_chisq) > 1e-6 ||
      abs(lr_p - original_increment$p_incremental) > 1e-6)
    stop("Cox model does not reproduce the frozen Step 45 estimates")

  model_list <- list("Clinical variables" = m0,
                     "Clinical + MRI2 volume" = m1,
                     "Clinical + MRI2 volume + prior loss" = m2)
  model_table <- rbindlist(lapply(names(model_list), function(label) {
    fit <- model_list[[label]]
    conc <- summary(fit)$concordance
    data.table(model = label, participants = fit$n, events = fit$nevent,
               apparent_c_index = unname(conc[1]), apparent_c_index_se = unname(conc[2]),
               log_likelihood = unname(fit$loglik[2]), AIC = AIC(fit))
  }))
  save_tsv(model_table, "Table3_same_sample_nested_Cox_models.tsv")
  save_tsv(data.table(contrast = "Prior loss beyond MRI2 volume and clinical variables",
                      participants = nrow(analysis), events = sum(analysis$event),
                      HR_per_1pct_year = hr, CI_low = ci[1], CI_high = ci[2],
                      Wald_p = wald_p, LR_chisq = lr, LR_p = lr_p,
                      apparent_delta_c_index =
                        model_table$apparent_c_index[3] - model_table$apparent_c_index[2]),
           "Table4_primary_increment_and_frozen_reproduction.tsv")
  save_tsv(data.table(term = rownames(cox.zph(m2)$table), cox.zph(m2)$table),
           "Table5_proportional_hazards_check.tsv")

  sensitivity <- sensitivity[status == "Fitted",
                             .(model, subjects, events, HR, CI_low, CI_high, p)]
  save_tsv(sensitivity, "Table6_event_definition_sensitivity.tsv")

  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  writeLines("STATUS=STEP38_COMPLETE", file.path(out, "38_COMPLETED.txt"))
  message("Step 38 complete: ", out)
  print(model_table)
})
