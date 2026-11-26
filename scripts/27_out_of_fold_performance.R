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
    oasis_cognition = file.path(.study$output,"44_OASIS3_full_clinical_link_results/09_B4_CDR_by_visit.csv")
  )
  if (any(!file.exists(paths)))
    stop("Missing frozen input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "history_increment_66",
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

  mri2 <- merge(clinical[prior_MCI_record == 0L],
                rbindlist(list(mayo_baseline, oasis_baseline)),
                by = c("cohort", "id"))
  mri2[, `:=`(time = impairment_time_years,
              event = as.integer(impairment_event))]
  updated[, `:=`(time = risk_years,
                 event = as.integer(impairment_event))]
  if (nrow(mri2[cohort == "Mayo"]) != 880L ||
      nrow(mri2[cohort == "OASIS-3"]) != 306L ||
      anyDuplicated(mri2, by = c("cohort", "id")) ||
      anyDuplicated(updated, by = c("cohort", "id")))
    stop("Step 65 risk-set identity changed.")

  formulas <- list(
    MRI2 = list(
      state = survival::Surv(time, event) ~ age_MRI2 + sex + education +
        hippo2_pct_ICV + baseline_score,
      history = survival::Surv(time, event) ~ age_MRI2 + sex + education +
        hippo2_pct_ICV + baseline_score + hippo_loss_pct_year),
    `Second cognitive visit` = list(
      state = survival::Surv(time, event) ~ age_landmark + sex + education +
        hippo2_pct_ICV + baseline_score + score + time_years,
      history = survival::Surv(time, event) ~ age_landmark + sex + education +
        hippo2_pct_ICV + baseline_score + score + time_years +
        hippo_loss_pct_year)
  )
  horizon <- 3
  thresholds <- c(0.05, 0.10, 0.15)
  folds <- 5L
  repeats <- 5L
  boot_reps <- 300L
  set.seed(66)

  predict_risk <- function(fit, newdata) {
    hazard <- survival::basehaz(fit, centered = FALSE)
    i <- which(hazard$time <= horizon)
    h0 <- if (length(i)) hazard$hazard[max(i)] else 0
    if (h0 == 0) return(rep(0, nrow(newdata)))
    lp <- as.numeric(predict(fit, newdata = newdata,
                             type = "lp", reference = "zero"))
    -expm1(-exp(log(h0) + lp))
  }
  censor_survival <- function(fit, times, left_limit = FALSE) {
    index <- findInterval(times - if (left_limit) 1e-9 else 0, fit$time)
    c(1, fit$surv)[index + 1L]
  }
  score_pair <- function(d, state, history) {
    censor_fit <- survival::survfit(
      survival::Surv(time, 1L - event) ~ 1, data = d)
    cases <- which(d$event == 1L & d$time <= horizon)
    controls <- which(d$time > horizon)
    if (!length(cases) || !length(controls)) return(rep(NA_real_, 9L))
    wc <- 1 / pmax(censor_survival(censor_fit, d$time[cases], TRUE), 1e-8)
    wn <- 1 / pmax(censor_survival(censor_fit,
                                    rep(horizon, length(controls))), 1e-8)
    brier <- function(risk) {
      (sum(wc * (1 - risk[cases])^2) +
         sum(wn * risk[controls]^2)) / nrow(d)
    }
    auc <- function(risk) {
      pair_score <- vapply(risk[cases], function(z)
        sum(wn * (z > risk[controls]) +
              0.5 * wn * (z == risk[controls])), numeric(1))
      sum(wc * pair_score) / (sum(wc) * sum(wn))
    }
    net_benefit <- function(risk, threshold) {
      (sum(wc * (risk[cases] >= threshold)) -
         threshold / (1 - threshold) *
         sum(wn * (risk[controls] >= threshold))) / nrow(d)
    }
    state_brier <- brier(state)
    history_brier <- brier(history)
    state_auc <- auc(state)
    history_auc <- auc(history)
    c(state_Brier = state_brier, history_Brier = history_brier,
      state_AUC = state_auc, history_AUC = history_auc,
      Brier_gain = state_brier - history_brier,
      AUC_gain = history_auc - state_auc,
      setNames(vapply(thresholds, function(t)
        net_benefit(history, t) - net_benefit(state, t), numeric(1)),
        paste0("net_benefit_gain_", 100 * thresholds, "pct")))
  }

  analyze <- function(d, cohort_name, landmark) {
    x <- as.data.frame(copy(d))
    x$sex <- factor(x$sex)
    form <- formulas[[landmark]]
    needed <- unique(c(all.vars(form$history), "id"))
    if (any(!complete.cases(x[, needed, drop = FALSE])) ||
        any(x$time <= 0)) stop("Missing field or invalid follow-up: ",
                               cohort_name, " / ", landmark)
    stratum <- ifelse(x$event == 1L & x$time <= horizon, "early_event",
               ifelse(x$time > horizon, "known_control", "early_censor"))
    state <- history <- numeric(nrow(x))
    fold_rows <- vector("list", repeats)
    for (repeat_id in seq_len(repeats)) {
      fold <- integer(nrow(x))
      for (s in unique(stratum)) {
        ids <- which(stratum == s)
        fold[ids] <- sample(rep(seq_len(folds), length.out = length(ids)))
      }
      pred_state <- pred_history <- rep(NA_real_, nrow(x))
      for (k in seq_len(folds)) {
        train <- x[fold != k, , drop = FALSE]
        test <- x[fold == k, , drop = FALSE]
        fit_state <- survival::coxph(form$state, data = train,
                                     ties = "efron", x = TRUE, model = TRUE)
        fit_history <- survival::coxph(form$history, data = train,
                                       ties = "efron", x = TRUE, model = TRUE)
        pred_state[fold == k] <- predict_risk(fit_state, test)
        pred_history[fold == k] <- predict_risk(fit_history, test)
      }
      if (any(!is.finite(pred_state)) || any(!is.finite(pred_history)))
        stop("Non-finite out-of-fold prediction: ", cohort_name, " / ", landmark)
      state <- state + pred_state / repeats
      history <- history + pred_history / repeats
      fold_rows[[repeat_id]] <- data.table(cohort = cohort_name,
                                     landmark = landmark, repeat_id = repeat_id,
                                     id = x$id, fold = fold)
    }
    point <- score_pair(x, state, history)
    boot <- replicate(boot_reps, {
      ids <- sample.int(nrow(x), replace = TRUE)
      score_pair(x[ids, , drop = FALSE], state[ids], history[ids])
    })
    if (any(!is.finite(boot)))
      stop("Bootstrap score unavailable: ", cohort_name, " / ", landmark)
    ci <- t(apply(boot, 1L, quantile, c(0.025, 0.975)))
    overall <- survival::survfit(survival::Surv(time, event) ~ 1, data = x)
    observed <- 1 - summary(overall, times = horizon, extend = TRUE)$surv[1L]
    result <- data.table(cohort = cohort_name, landmark = landmark,
                         participants = nrow(x), events_total = sum(x$event),
                         events_by_3y = sum(x$event == 1L & x$time <= horizon),
                         known_3y_status = sum(x$time > horizon |
                           (x$event == 1L & x$time <= horizon)),
                         observed_3y = observed,
                         mean_state_risk = mean(state),
                         mean_history_risk = mean(history))
    for (metric in names(point)) {
      result[, (metric) := point[[metric]]]
      result[, (paste0(metric, "_CI_low")) := ci[metric, 1L]]
      result[, (paste0(metric, "_CI_high")) := ci[metric, 2L]]
    }
    person <- data.table(cohort = cohort_name, landmark = landmark,
                         id = x$id, time = x$time, event = x$event,
                         state_risk_oof = state, history_risk_oof = history,
                         risk_change = history - state)
    list(summary = result, person = person,
         folds = rbindlist(fold_rows))
  }

  runs <- list(
    analyze(mri2[cohort == "Mayo"], "Mayo", "MRI2"),
    analyze(mri2[cohort == "OASIS-3"], "OASIS-3", "MRI2"),
    analyze(updated[cohort == "Mayo"], "Mayo", "Second cognitive visit"),
    analyze(updated[cohort == "OASIS-3"], "OASIS-3", "Second cognitive visit")
  )
  results <- rbindlist(lapply(runs, `[[`, "summary"))
  people <- rbindlist(lapply(runs, `[[`, "person"))
  fold_assignment <- rbindlist(lapply(runs, `[[`, "folds"))
  save_tsv(results, "Table1_paired_out_of_fold_increment.tsv")
  save_tsv(people, "Table2_person_level_out_of_fold_risk.tsv")
  save_tsv(fold_assignment, "Fold_assignment.tsv", "audit")

  cat("Step 66 complete: ", out, "\n", sep = "")
  print(results)
})
