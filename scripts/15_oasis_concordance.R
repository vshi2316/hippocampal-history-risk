local({
  packages <- c("data.table", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))

  library(data.table)
  library(survival)
  step45 <- file.path(.study$output,"45_OASIS3_frozen_event_link_results")
  step38 <- file.path(.study$output,"oasis_event_38/run_frozen")
  inputs <- c(
    risk = file.path(step45, "01_frozen_event_risk_set.csv"),
    subjects = file.path(step45, "07b_primary_model_subjects.csv"),
    model38 = file.path(step38, "tables/Table3_same_sample_nested_Cox_models.tsv"),
    primary38 = file.path(step38, "tables/Table4_primary_increment_and_frozen_reproduction.tsv"),
    complete38 = file.path(step38, "38_COMPLETED.txt")
  )
  if (any(!file.exists(inputs))) stop("Missing frozen input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  if (!"STATUS=STEP38_COMPLETE" %in% readLines(inputs[["complete38"]], warn = FALSE))
    stop("Step 38 is incomplete")

  out <- file.path(file.path(.study$output,"oasis_bootstrap_39"),
                   paste0("run_", "frozen"))
  for (part in c("tables", "source_data", "audit"))
    dir.create(file.path(out, part), recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(x, filename, part = "tables")
    fwrite(as.data.table(x), file.path(out, part, filename), sep = "\t", na = "NA")
  write_tsv(data.table(input = names(inputs), file = unname(inputs),
                       md5 = unname(tools::md5sum(inputs))),
            "Input_manifest.tsv", "audit")

  risk <- fread(inputs[["risk"]])
  frozen_ids <- fread(inputs[["subjects"]])$OASISID
  reference <- fread(inputs[["model38"]])
  reference_primary <- fread(inputs[["primary38"]])
  cn <- risk[baseline_stage == "CN"]
  cn[, event := as.logical(event)]
  fields <- c("OASISID", "event", "followup_days", "age_MRI2", "sex",
              "education", "baseline_CDRSB", "hippo2_pct_ICV", "hippo_loss_pct_year")
  analysis <- as.data.frame(cn[complete.cases(cn[, ..fields])])
  analysis$sex <- factor(analysis$sex)
  analysis$followup_years <- analysis$followup_days / 365.25
  if (nrow(analysis) != 312L || sum(analysis$event) != 52L ||
      anyDuplicated(analysis$OASISID) || !setequal(analysis$OASISID, frozen_ids))
    stop("The 312-person / 52-event model set differs from the frozen Step 45 cohort")

  base_formula <- Surv(followup_years, event) ~ age_MRI2 + sex + education +
    baseline_CDRSB + hippo2_pct_ICV
  full_formula <- update(base_formula, . ~ . + hippo_loss_pct_year)
  fit_base <- coxph(base_formula, data = analysis, ties = "efron")
  fit_full <- coxph(full_formula, data = analysis, ties = "efron")

  paired_c <- function(base, full, target) {
    lp_base <- as.numeric(predict(base, newdata = target, type = "lp"))
    lp_full <- as.numeric(predict(full, newdata = target, type = "lp"))
    if (any(!is.finite(lp_base)) || any(!is.finite(lp_full)))
      stop("Non-finite risk score")
    scores <- survival::concordance(
      Surv(followup_years, event) ~ lp_base + lp_full,
      data = target, reverse = TRUE)$concordance
    unname(scores)
  }

  apparent <- paired_c(fit_base, fit_full, analysis)
  frozen_c <- reference$apparent_c_index[match(
    c("Clinical + MRI2 volume", "Clinical + MRI2 volume + prior loss"),
    reference$model)]
  if (anyNA(frozen_c) || max(abs(apparent - frozen_c)) > 1e-6 ||
      abs(diff(apparent) - reference_primary$apparent_delta_c_index[1]) > 1e-6)
    stop("The apparent C-index does not reproduce the frozen Step 38 results")

  set.seed(20261113)
  B <- 1000L
  boot <- data.table(iteration = seq_len(B), status = rep("NOT_RUN", B),
                     unique_people = rep(NA_integer_, B),
                     events_in_sample = rep(NA_integer_, B),
                     base_train = rep(NA_real_, B), full_train = rep(NA_real_, B),
                     base_test = rep(NA_real_, B), full_test = rep(NA_real_, B))
  for (b in seq_len(B)) {
    indices <- sample.int(nrow(analysis), nrow(analysis), replace = TRUE)
    train <- analysis[indices, , drop = FALSE]
    boot[b, `:=`(unique_people = uniqueN(train$OASISID),
                events_in_sample = as.integer(sum(train$event)))]
    result <- tryCatch({
      b_base <- suppressWarnings(coxph(base_formula, data = train, ties = "efron"))
      b_full <- suppressWarnings(coxph(full_formula, data = train, ties = "efron"))
      if (any(!is.finite(coef(b_base))) || any(!is.finite(coef(b_full))))
        stop("Non-finite Cox coefficient")
      list(train = paired_c(b_base, b_full, train),
           test = paired_c(b_base, b_full, analysis))
    }, error = function(e) e)
    if (inherits(result, "error")) {
      boot[b, status := conditionMessage(result)]
    } else {
      boot[b, `:=`(status = "OK", base_train = result$train[1],
                  full_train = result$train[2], base_test = result$test[1],
                  full_test = result$test[2])]
    }
    if (b %% 100L == 0L) message("Bootstrap ", b, "/", B)
  }
  write_tsv(boot, "Table0_raw_bootstrap_replicates.tsv")
  boot[, `:=`(delta_train = full_train - base_train,
              delta_test = full_test - base_test,
              optimism_base = base_train - base_test,
              optimism_full = full_train - full_test)]
  boot[, optimism_delta := delta_train - delta_test]
  write_tsv(boot, "Table1_all_bootstrap_replicates.tsv")
  valid <- boot[status == "OK" & is.finite(optimism_delta)]
  failures <- boot[status != "OK", .N, by = status][order(-N)]
  write_tsv(failures, "Table2_bootstrap_fit_status.tsv")
  if (nrow(valid) < 0.9 * B)
    stop("Fewer than 90% of bootstrap replicates yielded finite paired C-index values")

  mean_optimism <- c(mean(valid$optimism_base), mean(valid$optimism_full))
  corrected <- apparent - mean_optimism
  apparent_delta <- diff(apparent)
  corrected_delta <- diff(corrected)
  percentile <- unname(quantile(valid$delta_train, c(0.025, 0.975)))
  corrected_basic <- 2 * apparent_delta - rev(percentile) -
    mean(valid$optimism_delta)
  discrimination <- data.table(
    metric = c("Clinical + MRI2 volume", "Clinical + MRI2 volume + prior loss",
               "Increment from prior loss"),
    apparent = c(apparent, apparent_delta),
    mean_bootstrap_optimism = c(mean_optimism, mean(valid$optimism_delta)),
    optimism_corrected = c(corrected, corrected_delta),
    bootstrap_95CI_low = c(NA_real_, NA_real_, corrected_basic[1]),
    bootstrap_95CI_high = c(NA_real_, NA_real_, corrected_basic[2]),
    successful_replicates = nrow(valid), attempted_replicates = B)
  write_tsv(discrimination, "Table3_optimism_corrected_C_index.tsv")
  write_tsv(data.table(
    estimand = "Increment from prior hippocampal loss",
    apparent_delta_C = apparent_delta,
    paired_bootstrap_apparent_95CI_low = percentile[1],
    paired_bootstrap_apparent_95CI_high = percentile[2],
    mean_optimism = mean(valid$optimism_delta),
    corrected_delta_C = corrected_delta,
    corrected_basic_95CI_low = corrected_basic[1],
    corrected_basic_95CI_high = corrected_basic[2],
    interval_method = "Basic bootstrap interval for apparent delta, shifted by mean optimism"),
    "Table4_increment_uncertainty.tsv")
  write_tsv(valid[, .(iteration, delta_train, delta_test, optimism_delta)],
            "Figure39_bootstrap_source.tsv", "source_data")
  write_tsv(discrimination, "Figure39_increment_estimates.tsv", "source_data")

  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  writeLines("STATUS=STEP39_COMPLETE", file.path(out, "39_COMPLETED.txt"))
  message("Step 39 complete: ", out)
  print(discrimination)
  print(failures)
})
