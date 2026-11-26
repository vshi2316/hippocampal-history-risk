local({
  packages <- c("data.table", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "))
  library(data.table)

  root <- .study$output
  paths <- c(
    clinical = file.path(root, "cn_progression_62/run_frozen/tables/Table1_person_level_clinical_course.tsv"),
    updated = file.path(root, "prediagnostic_64/run_frozen/tables/Table4_second_cognitive_visit_landmark.tsv"),
    mayo_risk = file.path(root, "mcsa_event_43/run_frozen/tables/Table2_person_level_frozen_risk_sets.tsv"),
    mayo_zip = file.path(.study$data,"mayo","Demographics_and_Assessments.zip"),
    oasis_risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    oasis_visits = file.path(.study$output,"44_OASIS3_full_clinical_link_results/07_D1_all_visit_stages.csv"),
    oasis_cognition = file.path(.study$output,"44_OASIS3_full_clinical_link_results/09_B4_CDR_by_visit.csv")
  )
  if (any(!file.exists(paths)))
    stop("Missing frozen input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "confirmed_impairment_67",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  save_tsv(data.table(input = names(paths), file = unname(paths),
                      md5 = unname(tools::md5sum(paths))),
           "Input_manifest.tsv", "audit")

  clinical <- fread(paths[["clinical"]])
  updated <- fread(paths[["updated"]])
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
                  FCOGIMPR = suppressWarnings(as.numeric(FCOGIMPR)),
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
  oasis_cog <- fread(paths[["oasis_cognition"]])
  oasis_baseline <- merge(oasis_cog[is.finite(mmse) & mmse >= 0 & mmse <= 30,
                                    .(id = OASISID, day,
                                      baseline_score = mmse)],
                          oasis_risk, by = "id")
  oasis_baseline[, gap_days := MRI2_day - day]
  oasis_baseline <- oasis_baseline[gap_days >= 0 & gap_days <= 365]
  setorder(oasis_baseline, id, gap_days, day)
  oasis_baseline <- unique(oasis_baseline, by = "id")[,
                              .(cohort = "OASIS-3", id, baseline_score)]
  mri2 <- merge(clinical[prior_MCI_record == 0L],
                rbindlist(list(mayo_baseline, oasis_baseline)),
                by = c("cohort", "id"))
  mri2[, `:=`(landmark = "MRI2", origin_day = 0,
              original_time = impairment_time_years,
              original_event = as.integer(impairment_event))]
  updated[, `:=`(landmark = "Second cognitive visit",
                 origin_day = time_years * 365.25,
                 original_time = risk_years,
                 original_event = as.integer(impairment_event))]
  if (nrow(mri2[cohort == "Mayo"]) != 880L ||
      nrow(mri2[cohort == "OASIS-3"]) != 306L ||
      anyDuplicated(mri2, by = c("cohort", "id")) ||
      anyDuplicated(updated, by = c("cohort", "id")))
    stop("Step 65 risk-set identity changed.")

  mayo_visits <- merge(mayo_raw[!is.na(VisitDate),
                                 .(id = MCSA_ID, VisitDate, FCOGIMPR)],
                       mayo_risk, by = "id", allow.cartesian = TRUE)
  mayo_visits[, day := as.numeric(VisitDate - MRI2_date)]
  mayo_visits <- mayo_visits[day > 0]
  first_unknown <- mayo_visits[FCOGIMPR == 4,
                               .(first_unknown_day = min(day)), by = id]
  mayo_visits <- merge(mayo_visits, first_unknown, by = "id", all.x = TRUE)
  mayo_visits <- mayo_visits[(is.na(first_unknown_day) |
                               day < first_unknown_day) &
                               FCOGIMPR %in% c(0, 1, 3)]
  mayo_visits[, stage := c("CN", "MCI", "Dementia")[
    match(FCOGIMPR, c(0, 1, 3))]]
  mayo_visits <- mayo_visits[, .(cohort = "Mayo", id, day, stage)]
  oasis_visits <- merge(fread(paths[["oasis_visits"]])[
    , .(id = OASISID, day, stage)], oasis_risk, by = "id",
    allow.cartesian = TRUE)
  oasis_visits[, day := day - MRI2_day]
  oasis_visits <- oasis_visits[day > 0,
                               .(cohort = "OASIS-3", id, day, stage)]
  visits <- rbindlist(list(mayo_visits, oasis_visits))
  if (anyDuplicated(visits, by = c("cohort", "id", "day")))
    stop("Multiple clinical stages on one visit date.")
  setkey(visits, cohort, id)

  classify_course <- function(r) {
    v <- copy(visits[.(r$cohort, r$id)])
    v <- v[day > r$origin_day]
    setorder(v, day)
    valid <- v$stage %in% c("CN", "MCI", "Dementia")
    impaired <- v$stage %in% c("MCI", "Dementia")
    if (!any(valid)) stop("No eligible follow-up visit: ", r$cohort, " / ", r$id)
    first <- which(impaired)[1L]
    first_day <- if (is.na(first)) NA_real_ else v$day[first]
    if (r$original_event == 1L &&
        (is.na(first_day) ||
         abs((first_day - r$origin_day) / 365.25 - r$original_time) > 1 / 365.25))
      stop("First event disagrees with Step 65: ", r$cohort, " / ", r$id)
    if (r$original_event == 0L &&
        any(impaired & v$day <= r$origin_day + r$original_time * 365.25 + 1))
      stop("Frozen non-event has an earlier impairment: ", r$cohort, " / ", r$id)
    second <- which(impaired[-1L] & impaired[-length(impaired)]) + 1L
    confirmed <- length(second) > 0L
    second <- if (confirmed) second[1L] else NA_integer_
    confirmed_day <- if (confirmed) v$day[second] else NA_real_
    last_valid_day <- max(v$day[valid])
    end_day <- if (confirmed) confirmed_day else last_valid_day
    if (end_day <= r$origin_day)
      stop("Non-positive follow-up: ", r$cohort, " / ", r$id)
    later_valid <- if (is.na(first)) 0L else
      sum(valid & v$day > first_day)
    reverted <- if (is.na(first)) FALSE else
      any(v$stage == "CN" & v$day > first_day &
            (!confirmed | v$day < confirmed_day))
    data.table(cohort = r$cohort, landmark = r$landmark, id = r$id,
               origin_day = r$origin_day,
               original_event = r$original_event,
               original_time = r$original_time,
               confirmed_event = as.integer(confirmed),
               confirmed_time = (end_day - r$origin_day) / 365.25,
               confirmed_stage_pair = if (confirmed)
                 paste(v$stage[second - 1L], v$stage[second], sep = " -> ")
                 else NA_character_,
               first_impairment_day = first_day,
               first_impairment_stage = if (is.na(first)) NA_character_
                 else v$stage[first],
               confirmation_lag_years = if (confirmed)
                 (confirmed_day - first_day) / 365.25 else NA_real_,
               reverted_before_confirmation = as.integer(reverted),
               valid_visits_after_first = later_valid,
               unconfirmed_at_last_visit = as.integer(!is.na(first) &&
                                                         !confirmed &&
                                                         later_valid == 0L),
               valid_visits_after_landmark = sum(valid))
  }

  risk_sets <- rbindlist(list(mri2, updated), fill = TRUE)
  course <- rbindlist(lapply(seq_len(nrow(risk_sets)), function(i)
    classify_course(risk_sets[i])))
  setkey(course, cohort, landmark, id)
  risk_sets <- merge(risk_sets, course,
                     by = c("cohort", "landmark", "id"),
                     suffixes = c("", "_course"))
  save_tsv(course, "Table1_person_level_confirmed_course.tsv")
  flow <- course[, .(
    participants = .N,
    first_impairment_events = sum(original_event),
    first_MCI = sum(first_impairment_stage == "MCI", na.rm = TRUE),
    first_direct_dementia = sum(first_impairment_stage == "Dementia", na.rm = TRUE),
    confirmed_events = sum(confirmed_event),
    confirmed_by_3y = sum(confirmed_event == 1L & confirmed_time <= 3),
    known_3y_status = sum(confirmed_time >= 3 |
                            (confirmed_event == 1L & confirmed_time <= 3)),
    confirmed_by_5y = sum(confirmed_event == 1L & confirmed_time <= 5),
    known_5y_status = sum(confirmed_time >= 5 |
                            (confirmed_event == 1L & confirmed_time <= 5)),
    first_with_later_classified_visit = sum(original_event == 1L &
                                             valid_visits_after_first > 0L),
    reverted_before_confirmation = sum(reverted_before_confirmation),
    unconfirmed_at_last_visit = sum(unconfirmed_at_last_visit),
    median_valid_followup_visits = as.numeric(median(valid_visits_after_landmark)),
    median_confirmation_lag_years = median(confirmation_lag_years,
                                           na.rm = TRUE)
  ), by = .(cohort, landmark)]
  save_tsv(flow, "Table2_event_confirmation_flow.tsv")
  save_tsv(course[confirmed_event == 1L,
                  .(participants = .N),
                  by = .(cohort, landmark, confirmed_stage_pair)],
           "Table2b_confirmed_stage_pairs.tsv")

  formulas <- list(
    MRI2 = survival::Surv(time, event) ~ age_MRI2 + sex + education +
      hippo2_pct_ICV + baseline_score + hippo_loss_pct_year,
    `Second cognitive visit` = survival::Surv(time, event) ~
      age_landmark + sex + education + hippo2_pct_ICV +
      baseline_score + score + time_years + hippo_loss_pct_year
  )
  cutpoints <- clinical[prior_MCI_record == 0L,
    .(loss_Q1 = unname(quantile(hippo_loss_pct_year, 0.25)),
      loss_Q3 = unname(quantile(hippo_loss_pct_year, 0.75))), by = cohort]
  horizon <- 3
  boot_reps <- 300L
  set.seed(67)
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
    c(risk_Q1 = mean(predict_risk(fit, low)),
      risk_Q3 = mean(predict_risk(fit, high)))
  }
  fit_one <- function(d, cohort_name, landmark) {
    x <- as.data.frame(copy(d))
    x$sex <- factor(x$sex)
    x$time <- x$original_time
    x$event <- x$original_event
    needed <- unique(c(all.vars(formulas[[landmark]]), "id"))
    if (any(!complete.cases(x[, needed, drop = FALSE])))
      stop("Incomplete model fields: ", cohort_name, " / ", landmark)
    first_fit <- survival::coxph(formulas[[landmark]], data = x,
                                 ties = "efron", x = TRUE, model = TRUE)
    first_beta <- unname(coef(first_fit)["hippo_loss_pct_year"])
    first_se <- sqrt(vcov(first_fit)["hippo_loss_pct_year",
                                      "hippo_loss_pct_year"])
    x$time <- x$confirmed_time
    x$event <- x$confirmed_event
    confirmed_fit <- survival::coxph(formulas[[landmark]], data = x,
                                     ties = "efron", x = TRUE, model = TRUE)
    beta <- unname(coef(confirmed_fit)["hippo_loss_pct_year"])
    se <- sqrt(vcov(confirmed_fit)["hippo_loss_pct_year",
                                   "hippo_loss_pct_year"])
    q <- cutpoints[cohort == cohort_name]
    point <- standardized_risks(confirmed_fit, x, q$loss_Q1, q$loss_Q3)
    boot <- replicate(boot_reps, {
      train <- x[sample.int(nrow(x), replace = TRUE), , drop = FALSE]
      tryCatch({
        refit <- survival::coxph(formulas[[landmark]], data = train,
                                 ties = "efron", x = TRUE, model = TRUE)
        risks <- standardized_risks(refit, train, q$loss_Q1, q$loss_Q3)
        c(risks, risk_difference = unname(risks[2L] - risks[1L]))
      }, error = function(e) rep(NA_real_, 3L))
    })
    boot <- t(boot)
    colnames(boot) <- c("risk_Q1", "risk_Q3", "risk_difference")
    boot <- boot[complete.cases(boot), , drop = FALSE]
    if (nrow(boot) < 250L)
      stop("Fewer than 250 valid bootstrap fits: ", cohort_name,
           " / ", landmark, "; valid=", nrow(boot))
    ci <- apply(boot, 2L, quantile, c(0.025, 0.975))
    sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = x)
    observed <- 1 - summary(sf, times = horizon, extend = TRUE)$surv[1L]
    data.table(cohort = cohort_name, landmark = landmark,
               participants = nrow(x),
               first_events = sum(x$original_event),
               confirmed_events = sum(x$event),
               confirmed_events_by_3y = sum(x$event == 1L & x$time <= horizon),
               first_HR_per_1pct_year = exp(first_beta),
               first_HR_CI_low = exp(first_beta - 1.96 * first_se),
               first_HR_CI_high = exp(first_beta + 1.96 * first_se),
               confirmed_HR_per_1pct_year = exp(beta),
               HR_CI_low = exp(beta - 1.96 * se),
               HR_CI_high = exp(beta + 1.96 * se),
               HR_p = 2 * pnorm(-abs(beta / se)),
               risk_Q1 = point[1L], risk_Q1_CI_low = ci[1L, "risk_Q1"],
               risk_Q1_CI_high = ci[2L, "risk_Q1"],
               risk_Q3 = point[2L], risk_Q3_CI_low = ci[1L, "risk_Q3"],
               risk_Q3_CI_high = ci[2L, "risk_Q3"],
               risk_difference = point[2L] - point[1L],
               difference_CI_low = ci[1L, "risk_difference"],
               difference_CI_high = ci[2L, "risk_difference"],
               observed_3y = observed,
               bootstrap_success = nrow(boot))
  }
  results <- rbindlist(list(
    fit_one(risk_sets[cohort == "Mayo" & landmark == "MRI2"], "Mayo", "MRI2"),
    fit_one(risk_sets[cohort == "OASIS-3" & landmark == "MRI2"], "OASIS-3", "MRI2"),
    fit_one(risk_sets[cohort == "Mayo" & landmark == "Second cognitive visit"],
            "Mayo", "Second cognitive visit"),
    fit_one(risk_sets[cohort == "OASIS-3" & landmark == "Second cognitive visit"],
            "OASIS-3", "Second cognitive visit")
  ))
  save_tsv(results, "Table3_confirmed_event_association_and_3y_risk.tsv")
  writeLines(capture.output(sessionInfo()),
             file.path(out, "audit/R_sessionInfo.txt"))
  cat("Step 67 complete: ", out, "\n", sep = "")
  print(flow)
  print(results)
})
