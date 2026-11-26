local({
  packages <- c("data.table", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
  library(data.table)
  root <- .study$output
  paths <- c(
    clinical = file.path(root, "cn_progression_62/run_frozen/tables/Table1_person_level_clinical_course.tsv"),
    updated = file.path(root, "prediagnostic_64/run_frozen/tables/Table4_second_cognitive_visit_landmark.tsv"),
    updated_model = file.path(root, "prediagnostic_64/run_frozen/tables/Table5_updated_landmark_first_impairment.tsv"),
    mayo_risk = file.path(root, "mcsa_event_43/run_frozen/tables/Table2_person_level_frozen_risk_sets.tsv"),
    mayo_zip = file.path(.study$data,"mayo","Demographics_and_Assessments.zip"),
    oasis_risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    oasis_cognition = file.path(.study$output,"44_OASIS3_full_clinical_link_results/09_B4_CDR_by_visit.csv")
  )
  if (any(!file.exists(paths)))
    stop("Missing frozen input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "absolute_risk_65",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  save_tsv(data.table(input = names(paths), file = unname(paths),
                      md5 = unname(tools::md5sum(paths))),
           "Input_manifest.tsv", "audit")

  clinical <- fread(paths[["clinical"]])
  updated <- fread(paths[["updated"]])
  reference <- fread(paths[["updated_model"]])
  if (nrow(clinical[cohort == "Mayo"]) != 925L ||
      nrow(clinical[cohort == "OASIS-3"]) != 312L ||
      nrow(updated[cohort == "Mayo"]) != 527L ||
      nrow(updated[cohort == "OASIS-3"]) != 219L)
    stop("Frozen analysis populations changed.")

  mayo_risk <- fread(paths[["mayo_risk"]])[stage == "CN",
                       .(id, MRI2_date = as.Date(MRI2_date))]
  mayo_raw <- as.data.table(read.csv(unz(paths[["mayo_zip"]],
                                      "MCSA_Data_12Nov2026.csv"),
                                    stringsAsFactors = FALSE, check.names = FALSE))
  mayo_raw[, `:=`(VisitDate = as.Date(VisitDate),
                  pzglobal = suppressWarnings(as.numeric(pzglobal)))]
  mayo_baseline <- merge(mayo_raw[!is.na(VisitDate) & is.finite(pzglobal),
                                    .(id = MCSA_ID, VisitDate,
                                      baseline_score = pzglobal)],
                         mayo_risk, by = "id")
  mayo_baseline <- mayo_baseline[,
    .(baseline_score = median(baseline_score), MRI2_date = MRI2_date[1L]),
    by = .(id, VisitDate)]
  mayo_baseline[, gap_days := as.numeric(MRI2_date - VisitDate)]
  mayo_baseline <- mayo_baseline[gap_days >= 0 & gap_days <= 90]
  setorder(mayo_baseline, id, gap_days, VisitDate)
  mayo_baseline <- unique(mayo_baseline, by = "id")[,
                             .(cohort = "Mayo", id, baseline_score)]

  oasis_risk <- fread(paths[["oasis_risk"]])[baseline_stage == "CN",
                         .(id = OASISID, MRI2_day)]
  oasis_visits <- fread(paths[["oasis_cognition"]])
  oasis_baseline <- merge(oasis_visits[is.finite(mmse) & mmse >= 0 &
                                          mmse <= 30,
                                          .(id = OASISID, day,
                                            baseline_score = mmse)],
                          oasis_risk, by = "id")
  oasis_baseline[, gap_days := MRI2_day - day]
  oasis_baseline <- oasis_baseline[gap_days >= 0 & gap_days <= 365]
  setorder(oasis_baseline, id, gap_days, day)
  oasis_baseline <- unique(oasis_baseline, by = "id")[,
                              .(cohort = "OASIS-3", id, baseline_score)]
  cognition <- rbindlist(list(mayo_baseline, oasis_baseline))
  mri2 <- merge(clinical[prior_MCI_record == 0L], cognition,
                by = c("cohort", "id"))
  mri2[, `:=`(landmark = "MRI2",
              time = impairment_time_years,
              event = as.integer(impairment_event))]
  updated[, `:=`(landmark = "Second cognitive visit",
                 time = risk_years,
                 event = as.integer(impairment_event))]
  if (anyDuplicated(mri2, by = c("cohort", "id")) ||
      anyDuplicated(updated, by = c("cohort", "id")) ||
      any(mri2$time <= 0) || any(updated$time <= 0))
    stop("Risk-set identity or event time changed.")
  if (nrow(mri2[cohort == "Mayo"]) != 880L ||
      nrow(mri2[cohort == "OASIS-3"]) != 306L)
    stop("MRI2 baseline cognition coverage changed.")
  if (any(!complete.cases(mri2[, .(age_MRI2, sex, education,
                                    hippo2_pct_ICV, baseline_score,
                                    hippo_loss_pct_year, time, event)])) ||
      any(!complete.cases(updated[, .(age_landmark, sex, education,
                                       hippo2_pct_ICV, baseline_score,
                                       score, time_years,
                                       hippo_loss_pct_year, time, event)])))
    stop("A required model field is missing.")

  cutpoints <- clinical[prior_MCI_record == 0L,
                        .(loss_Q1 = unname(quantile(hippo_loss_pct_year, 0.25)),
                          loss_Q3 = unname(quantile(hippo_loss_pct_year, 0.75))),
                        by = cohort]
  save_tsv(cutpoints, "Table1_frozen_prior_loss_quartiles.tsv")

  formulas <- list(
    MRI2 = survival::Surv(time, event) ~ age_MRI2 + sex + education +
      hippo2_pct_ICV + baseline_score + hippo_loss_pct_year,
    `Second cognitive visit` = survival::Surv(time, event) ~
      age_landmark + sex + education + hippo2_pct_ICV +
      baseline_score + score + time_years + hippo_loss_pct_year
  )
  horizon <- 3
  boot_reps <- 300L
  set.seed(65)

  predict_risk <- function(fit, newdata) {
    hazard <- survival::basehaz(fit, centered = FALSE)
    i <- which(hazard$time <= horizon)
    h0 <- if (length(i)) hazard$hazard[max(i)] else 0
    if (h0 == 0) return(rep(0, nrow(newdata)))
    lp <- as.numeric(predict(fit, newdata = newdata,
                             type = "lp", reference = "zero"))
    -expm1(-exp(log(h0) + lp))
  }
  standardized_risks <- function(fit, d, q1, q3) {
    low <- copy(as.data.table(d))
    high <- copy(low)
    low[, hippo_loss_pct_year := q1]
    high[, hippo_loss_pct_year := q3]
    r1 <- mean(predict_risk(fit, low))
    r3 <- mean(predict_risk(fit, high))
    c(risk_Q1 = r1, risk_Q3 = r3, risk_difference = r3 - r1)
  }
  censor_survival <- function(fit, times, left_limit = FALSE) {
    index <- findInterval(times - if (left_limit) 1e-9 else 0, fit$time)
    c(1, fit$surv)[index + 1L]
  }
  brier_3y <- function(d, risk) {
    censor_fit <- survival::survfit(
      survival::Surv(time, 1L - event) ~ 1, data = d)
    occurred <- d$event == 1L & d$time <= horizon
    survived <- d$time > horizon
    weights <- numeric(nrow(d))
    weights[occurred] <- 1 / pmax(
      censor_survival(censor_fit, d$time[occurred], TRUE), 1e-8)
    weights[survived] <- 1 / pmax(
      censor_survival(censor_fit, rep(horizon, sum(survived))), 1e-8)
    mean(weights * (as.numeric(occurred) - risk)^2)
  }
  observed_3y <- function(d) {
    sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = d)
    at <- summary(sf, times = horizon, extend = TRUE)
    c(observed = 1 - at$surv[1L],
      CI_low = 1 - at$upper[1L], CI_high = 1 - at$lower[1L])
  }

  analyze <- function(d, cohort_name, landmark) {
    x <- as.data.frame(copy(d))
    x$sex <- factor(x$sex)
    form <- formulas[[landmark]]
    fit <- survival::coxph(form, data = x, ties = "efron", x = TRUE)
    q <- cutpoints[which(cutpoints[["cohort"]] == cohort_name)]
    if (nrow(q) != 1L)
      stop("Expected one frozen quartile row for ", cohort_name)
    point <- standardized_risks(fit, x, q$loss_Q1, q$loss_Q3)
    beta <- unname(coef(fit)["hippo_loss_pct_year"])
    se <- sqrt(vcov(fit)["hippo_loss_pct_year", "hippo_loss_pct_year"])
    if (landmark == "Second cognitive visit") {
      old <- reference[which(reference[["cohort"]] == cohort_name)]
      if (nrow(old) != 1L ||
          abs(exp(beta) - old$HR_per_1pct_year) > 1e-6)
        stop("Step 64 updated-landmark HR did not reproduce: ", cohort_name)
    }
    predicted <- predict_risk(fit, x)
    apparent_brier <- brier_3y(x, predicted)
    overall <- observed_3y(x)

    boot <- vector("list", boot_reps)
    boot_errors <- character(boot_reps)
    for (j in seq_len(boot_reps)) {
      draw <- sample.int(nrow(x), replace = TRUE)
      train <- x[draw, , drop = FALSE]
      result <- tryCatch({
        refit <- survival::coxph(form, data = train, ties = "efron",
                                 x = TRUE, model = TRUE)
        risk_train <- predict_risk(refit, train)
        risk_original <- predict_risk(refit, x)
        c(standardized_risks(refit, train, q$loss_Q1, q$loss_Q3),
          brier_optimism = brier_3y(x, risk_original) -
            brier_3y(train, risk_train))
      }, error = function(e) {
        boot_errors[j] <<- conditionMessage(e)
        rep(NA_real_, 4L)
      })
      boot[[j]] <- result
    }
    boot <- as.data.table(do.call(rbind, boot))
    setnames(boot, c("risk_Q1", "risk_Q3", "risk_difference",
                     "brier_optimism"))
    boot <- boot[complete.cases(boot) &
                   risk_Q1 >= 0 & risk_Q1 <= 1 &
                   risk_Q3 >= 0 & risk_Q3 <= 1]
    if (nrow(boot) < 250L) {
      reasons <- sort(table(boot_errors[nzchar(boot_errors)]),
                      decreasing = TRUE)
      stop("Fewer than 250 valid bootstrap fits: ", cohort_name, " / ",
           landmark, "; valid=", nrow(boot), "; leading error=",
           if (length(reasons)) names(reasons)[1L] else "non-finite result")
    }
    interval <- function(v) unname(quantile(v, c(0.025, 0.975)))
    q1_ci <- interval(boot$risk_Q1)
    q3_ci <- interval(boot$risk_Q3)
    diff_ci <- interval(boot$risk_difference)
    summary_row <- data.table(
      cohort = cohort_name, landmark = landmark, participants = nrow(x),
      events_total = sum(x$event),
      events_by_3y = sum(x$event == 1L & x$time <= horizon),
      known_3y_status = sum(x$time >= horizon |
                              (x$event == 1L & x$time <= horizon)),
      loss_Q1 = q$loss_Q1, loss_Q3 = q$loss_Q3,
      HR_per_1pct_year = exp(beta),
      HR_CI_low = exp(beta - 1.96 * se),
      HR_CI_high = exp(beta + 1.96 * se),
      risk_Q1 = point["risk_Q1"], risk_Q1_CI_low = q1_ci[1L],
      risk_Q1_CI_high = q1_ci[2L],
      risk_Q3 = point["risk_Q3"], risk_Q3_CI_low = q3_ci[1L],
      risk_Q3_CI_high = q3_ci[2L],
      risk_difference = point["risk_difference"],
      difference_CI_low = diff_ci[1L],
      difference_CI_high = diff_ci[2L],
      observed_3y = overall["observed"],
      mean_predicted_3y = mean(predicted),
      apparent_Brier_3y = apparent_brier,
      optimism_corrected_Brier_3y =
        apparent_brier + mean(boot$brier_optimism),
      bootstrap_success = nrow(boot))

    order <- rank(predicted, ties.method = "first")
    group <- pmin(3L, ceiling(3 * order / length(order)))
    groups <- rbindlist(lapply(1:3, function(g) {
      indices <- which(group == g)
      z <- x[indices, , drop = FALSE]
      obs <- observed_3y(z)
      data.table(cohort = cohort_name, landmark = landmark,
                 risk_tertile = g, participants = nrow(z),
                 events_by_3y = sum(z$event == 1L & z$time <= horizon),
                 known_3y_status = sum(z$time >= horizon |
                   (z$event == 1L & z$time <= horizon)),
                 mean_predicted_3y = mean(predicted[indices]),
                 observed_3y = obs["observed"],
                 observed_CI_low = obs["CI_low"],
                 observed_CI_high = obs["CI_high"])
    }))
    list(summary = summary_row, calibration = groups,
         person = data.table(cohort = cohort_name, landmark = landmark,
                             id = x$id, time = x$time,
                             event = x$event, predicted_3y = predicted,
                             risk_tertile = group), bootstrap = boot)
  }

  work <- list(
    analyze(mri2[cohort == "Mayo"], "Mayo", "MRI2"),
    analyze(mri2[cohort == "OASIS-3"], "OASIS-3", "MRI2"),
    analyze(updated[cohort == "Mayo"], "Mayo", "Second cognitive visit"),
    analyze(updated[cohort == "OASIS-3"], "OASIS-3", "Second cognitive visit")
  )
  results <- rbindlist(lapply(work, `[[`, "summary"))
  calibration <- rbindlist(lapply(work, `[[`, "calibration"))
  person <- rbindlist(lapply(work, `[[`, "person"))
  save_tsv(results, "Table2_three_year_absolute_risk.tsv")
  save_tsv(calibration, "Table3_three_year_calibration.tsv")
  save_tsv(person, "Table4_person_level_apparent_risk.tsv")
  writeLines(capture.output(sessionInfo()),
             file.path(out, "audit", "R_sessionInfo.txt"))
  message("Step 65 complete: ", out)
  print(results, row.names = FALSE)
  print(calibration, row.names = FALSE)
})
