local({
  needed_packages <- c("data.table", "survival")
  missing_packages <- needed_packages[!vapply(needed_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
  library(data.table)

  root <- .study$output
  audit69 <- file.path(root, "clinical_field_audit_69/run_frozen")
  inputs <- c(
    completed69 = file.path(audit69, "69_COMPLETED.txt"),
    risk = file.path(audit69, "source_data/Frozen_MRI2_risk_with_baseline_score.tsv"),
    apoe = file.path(audit69, "source_data/Table4_person_level_APOE_availability.tsv"),
    anatomy = file.path(audit69, "source_data/Table6_exact_scan_anatomical_pairs.tsv"),
    intervals = file.path(audit69, "source_data/MRI_interval_in_final_risk_set.tsv"),
    scan_identity = file.path(audit69, "audit/OASIS_selected_scan_identity.tsv"),
    scanner = file.path(audit69, "source_data/OASIS_frozen_scan_metadata.tsv"),
    muse = file.path(root, "oasis_tau_mri_gate_27/run_frozen/audit/input_csv_cache/MUSE_volume_OASIS_04Nov2026.csv"),
    mr_metadata = file.path(root, "oasis_tau_mri_gate_27/run_frozen/audit/input_csv_cache/OASIS3_MR_scan_summary.csv"),
    quartiles = file.path(root, "absolute_risk_65/run_frozen/tables/Table1_frozen_prior_loss_quartiles.tsv"),
    oof = file.path(root, "history_increment_66/run_frozen/tables/Table2_person_level_out_of_fold_risk.tsv"),
    increment66 = file.path(root, "history_increment_66/run_frozen/tables/Table1_paired_out_of_fold_increment.tsv")
  )
  required_inputs <- inputs
  if (any(!file.exists(required_inputs))) stop("Missing frozen input: ", paste(required_inputs[!file.exists(required_inputs)], collapse = "; "))
  out <- file.path(root, "structural_history_70", paste0("run_", "frozen"))
  if (dir.exists(out)) stop("Output directory already exists: ", out)
  for (folder in c("tables", "source_data", "audit")) dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables") fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  number <- function(x) suppressWarnings(as.numeric(x))
  horizon <- 3
  model_boot_reps <- 500L
  screening_boot_reps <- 1000L
  screening_fraction <- 0.20
  set.seed(70)

  save_tsv(data.table(input = names(required_inputs), file = unname(required_inputs), md5 = unname(tools::md5sum(required_inputs))), "Input_manifest.tsv", "audit")

  risk <- fread(inputs[["risk"]])
  if (anyDuplicated(risk, by = c("cohort", "id"))) stop("Step 69 risk table has duplicated people.")
  full_counts <- risk[, .(participants = .N, events = sum(impairment_event)), by = cohort]
  if (nrow(risk[cohort == "Mayo"]) != 925L || sum(risk[cohort == "Mayo", impairment_event]) != 177L ||
      nrow(risk[cohort == "OASIS-3"]) != 306L || sum(risk[cohort == "OASIS-3", impairment_event]) != 59L)
    stop("Frozen first-impairment populations changed.")
  d <- copy(risk[cognition_complete %in% TRUE])
  d[, `:=`(time = impairment_time_years, event = as.integer(impairment_event))]
  d[, sex := factor(sex)]
  if (nrow(d[cohort == "Mayo"]) != 880L || sum(d[cohort == "Mayo", event]) != 168L ||
      nrow(d[cohort == "OASIS-3"]) != 306L || sum(d[cohort == "OASIS-3", event]) != 59L)
    stop("Cognition-complete MRI2 population changed.")
  apoe <- fread(inputs[["apoe"]])[, .(cohort, id, apoe_e4)]
  intervals <- fread(inputs[["intervals"]])[, .(cohort, id, interval_years)]
  if (anyDuplicated(apoe, by = c("cohort", "id")) || anyDuplicated(intervals, by = c("cohort", "id"))) stop("Duplicate APOE or interval rows.")
  d <- merge(merge(d, apoe, by = c("cohort", "id"), all.x = TRUE), intervals, by = c("cohort", "id"), all.x = TRUE)
  if (any(!is.finite(d$interval_years) | d$interval_years <= 0)) stop("Invalid MRI interval.")
  anatomy <- fread(inputs[["anatomy"]])[field == "lateral_ventricles" & usable_pair %in% TRUE,
    .(cohort, id, ventricle2_pct_ICV = 100 * value_MRI2 / ICV_MRI2, ventricle_expansion_pct_year = annual_log_change_pct)]
  if (anyDuplicated(anatomy, by = c("cohort", "id"))) stop("Duplicate ventricular pair.")
  d <- merge(d, anatomy, by = c("cohort", "id"), all.x = TRUE)
  d[, interval_log2 := log2(interval_years / 2)]
  d[, model_complete := complete.cases(.SD), .SDcols = c("time", "event", "age_MRI2", "sex", "education", "hippo2_pct_ICV", "baseline_score", "hippo_loss_pct_year")]
  if (any(!d$model_complete) || any(d$time <= 0)) stop("Missing primary model variables or invalid follow-up.")
  qtable <- fread(inputs[["quartiles"]])
  if (anyDuplicated(qtable$cohort) || !all(c("loss_Q1", "loss_Q3") %in% names(qtable))) stop("Invalid frozen quartiles.")
  q_for <- function(cohort_name) {
    z <- qtable[which(qtable[["cohort"]] == cohort_name)]
    if (nrow(z) != 1L || !is.finite(z$loss_Q1) || !is.finite(z$loss_Q3) || z$loss_Q1 >= z$loss_Q3) stop("Invalid quartile row: ", cohort_name)
    c(q1 = z$loss_Q1[[1L]], q3 = z$loss_Q3[[1L]])
  }
  save_tsv(d, "MRI2_complete_model_inputs.tsv", "source_data")
  save_tsv(rbindlist(list(full_counts[, population := "Full_first_impairment"],
    d[, .(participants = .N, events = sum(event)), by = cohort][, population := "Cognition_complete_MRI2"],
    d[apoe_e4 %in% c(0, 1), .(participants = .N, events = sum(event)), by = cohort][, population := "APOE_complete_MRI2"])), "Table1_frozen_populations.tsv")

  state_formula <- survival::Surv(time, event) ~ age_MRI2 + sex + education + hippo2_pct_ICV + baseline_score
  history_formula <- update(state_formula, . ~ . + hippo_loss_pct_year)
  apoe_formula <- update(history_formula, . ~ . + apoe_e4)
  vent_formula <- update(apoe_formula, . ~ . + ventricle2_pct_ICV + ventricle_expansion_pct_year)
  interval_formula <- update(apoe_formula, . ~ . + interval_log2)
  interaction_formula <- update(interval_formula, . ~ . + hippo_loss_pct_year:interval_log2)
  fit_cox <- function(form, x) {
    # Cohorts use different sex codes; remove empty levels again in each bootstrap draw.
    fit_data <- droplevels(as.data.frame(x))
    fit <- withCallingHandlers(survival::coxph(form, data = fit_data, ties = "efron", x = TRUE, model = TRUE),
      warning = function(w) {
        if (grepl("infinite|converg|singular", conditionMessage(w), ignore.case = TRUE)) stop(conditionMessage(w))
      })
    if (any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit)))) {
      bad_terms <- unique(c(names(coef(fit))[!is.finite(coef(fit))],
        colnames(vcov(fit))[colSums(!is.finite(vcov(fit))) > 0L]))
      stop("Non-finite Cox coefficient or variance: ", paste(bad_terms, collapse = ", "))
    }
    fit
  }
  risk3 <- function(fit, x) {
    reference <- as.data.frame(x[1L, ])
    sf <- survival::survfit(fit, newdata = reference, se.fit = FALSE)
    pos <- which(sf$time <= horizon)
    href <- if (length(pos)) as.numeric(sf$cumhaz)[max(pos)] else 0
    if (href == 0) return(rep(0, nrow(x)))
    lp <- as.numeric(predict(fit, newdata = as.data.frame(x), type = "lp", reference = "zero"))
    reference_lp <- as.numeric(predict(fit, newdata = reference, type = "lp", reference = "zero"))
    -expm1(-exp(log(href) + lp - reference_lp))
  }
  standardized <- function(fit, x, q) {
    low <- copy(as.data.table(x)); high <- copy(low)
    low[, hippo_loss_pct_year := as.numeric(q[["q1"]])]
    high[, hippo_loss_pct_year := as.numeric(q[["q3"]])]
    r1 <- mean(risk3(fit, low)); r3 <- mean(risk3(fit, high))
    c(beta_loss = unname(coef(fit)["hippo_loss_pct_year"]), risk_Q1 = r1, risk_Q3 = r3, risk_difference = r3 - r1)
  }
  coefficients <- function(fit, cohort_name, comparison, model_name, x) {
    tab <- as.data.frame(summary(fit)$coefficients)
    data.table(cohort = cohort_name, comparison, model = model_name, participants = nrow(x), events_total = sum(x$event),
      term = rownames(tab), beta = tab[["coef"]], SE = tab[["se(coef)"]],
      HR = exp(tab[["coef"]]), CI_low = exp(tab[["coef"]] - 1.96 * tab[["se(coef)"]]),
      CI_high = exp(tab[["coef"]] + 1.96 * tab[["se(coef)"]]), P = tab[["Pr(>|z|)"]])
  }
  boot_errors <- list()
  ci <- function(z) {
    valid <- z[is.finite(z)]
    if (length(valid) < 0.9 * length(z)) return(c(NA_real_, NA_real_))
    as.numeric(quantile(valid, c(.025, .975), names = FALSE))
  }
  pair_models <- function(x, cohort_name, comparison, base_form, extended_form, alternate = NULL) {
    q <- q_for(cohort_name)
    other <- if (is.null(alternate)) copy(x) else copy(alternate)
    if (!identical(x$id, other$id) || !identical(x$time, other$time) || !identical(x$event, other$event)) stop("Paired model participants or follow-up differ.")
    base_fit <- fit_cox(base_form, x); extended_fit <- fit_cox(extended_form, other)
    p0 <- standardized(base_fit, x, q); p1 <- standardized(extended_fit, other, q)
    errors <- rep("", model_boot_reps)
    draws <- t(vapply(seq_len(model_boot_reps), function(b) {
      idx <- sample.int(nrow(x), nrow(x), replace = TRUE)
      tryCatch({
        z0 <- x[idx]; z1 <- other[idx]
        v0 <- standardized(fit_cox(base_form, z0), z0, q)
        v1 <- standardized(fit_cox(extended_form, z1), z1, q)
        c(v0, v1)
      }, error = function(e) { errors[b] <<- conditionMessage(e); rep(NA_real_, 8L) })
    }, numeric(8L)))
    base_draws <- draws[, 1:4, drop = FALSE]; ext_draws <- draws[, 5:8, drop = FALSE]
    colnames(base_draws) <- colnames(ext_draws) <- names(p0)
    rows <- rbindlist(lapply(seq_along(p0), function(j) {
      c0 <- ci(base_draws[, j]); c1 <- ci(ext_draws[, j]); cd <- ci(ext_draws[, j] - base_draws[, j])
      data.table(cohort = cohort_name, comparison, metric = names(p0)[j], participants = nrow(x), events_total = sum(x$event),
        events_by_3y = sum(x$event == 1L & x$time <= horizon), loss_Q1 = q[["q1"]], loss_Q3 = q[["q3"]],
        base = p0[j], base_CI_low = c0[1L], base_CI_high = c0[2L], extended = p1[j], extended_CI_low = c1[1L], extended_CI_high = c1[2L],
        paired_change = p1[j] - p0[j], change_CI_low = cd[1L], change_CI_high = cd[2L],
        bootstrap_requested = model_boot_reps, bootstrap_valid = sum(is.finite(base_draws[, j]) & is.finite(ext_draws[, j])),
        CI_status = if (all(is.finite(c(c0, c1, cd)))) "AVAILABLE" else "BOOTSTRAP_STABILITY_REVIEW")
    }))
    if (any(nzchar(errors))) boot_errors[[length(boot_errors) + 1L]] <<- data.table(cohort = cohort_name, comparison, replicate = which(nzchar(errors)), error = errors[nzchar(errors)])
    list(summary = rows, coefficients = rbindlist(list(coefficients(base_fit, cohort_name, comparison, "base", x),
      coefficients(extended_fit, cohort_name, comparison, "extended", other))), base_fit = base_fit, extended_fit = extended_fit)
  }

  cat("Step 70: identical-sample clinical comparisons\n")
  paired <- list(); grids <- list()
  for (cohort_name in c("Mayo", "OASIS-3")) {
    x <- d[d[["cohort"]] == cohort_name & apoe_e4 %in% c(0, 1)]
    paired[[paste0(cohort_name, "_APOE")]] <- pair_models(x, cohort_name, "APOE_adjustment", history_formula, apoe_formula)
    paired[[paste0(cohort_name, "_interval")]] <- pair_models(x, cohort_name, "Continuous_MRI_interval_interaction", interval_formula, interaction_formula)
    fit <- paired[[paste0(cohort_name, "_interval")]]$extended_fit
    term <- intersect(c("hippo_loss_pct_year:interval_log2", "interval_log2:hippo_loss_pct_year"), names(coef(fit)))
    if (length(term) != 1L) stop("Cannot identify the MRI interval interaction.")
    vals <- as.numeric(quantile(x$interval_years, c(.25, .5, .75), names = FALSE))
    grids[[cohort_name]] <- rbindlist(lapply(seq_along(vals), function(k) {
      v <- setNames(rep(0, length(coef(fit))), names(coef(fit)))
      v["hippo_loss_pct_year"] <- 1; v[term] <- log2(vals[k] / 2)
      beta <- sum(v * coef(fit)); se <- sqrt(as.numeric(t(v) %*% vcov(fit) %*% v))
      data.table(cohort = cohort_name, interval_percentile = c(25L, 50L, 75L)[k], interval_years = vals[k],
        participants = nrow(x), events_total = sum(x$event), HR_per_1pct_year = exp(beta), CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se),
        interaction_beta = unname(coef(fit)[term]), interaction_P = summary(fit)$coefficients[term, "Pr(>|z|)"])
    }))
    if (cohort_name == "OASIS-3") {
      xv <- x[is.finite(ventricle2_pct_ICV) & is.finite(ventricle_expansion_pct_year)]
      paired[["OASIS_ventricles"]] <- pair_models(xv, cohort_name, "Lateral_ventricle_adjustment", apoe_formula, vent_formula)
    }
  }
  save_tsv(rbindlist(lapply(paired, `[[`, "summary")), "Table2_same_sample_standardized_risk_comparisons.tsv")
  save_tsv(rbindlist(lapply(paired, `[[`, "coefficients")), "Table3_same_sample_Cox_coefficients.tsv")
  grid <- rbindlist(grids)
  ip <- unique(grid[, .(cohort, interaction_P)])
  ip[, interaction_P_Holm_two_cohorts := p.adjust(interaction_P, method = "holm")]
  save_tsv(merge(grid, ip, by = c("cohort", "interaction_P")), "Table4_continuous_interval_conditional_HR.tsv")

  # Descriptive groups retain all three-year censoring and continuous-model comparisons.
  km3 <- function(x) {
    sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = as.data.frame(x), conf.type = "log-log")
    if (max(x$time) < horizon && tail(sf$surv, 1L) > 0) return(c(risk = NA_real_, low = NA_real_, high = NA_real_))
    idx <- which(sf$time <= horizon)
    if (!length(idx)) return(c(risk = 0, low = 0, high = 0))
    j <- max(idx)
    if (max(x$time) < horizon && sf$surv[j] > 0) return(c(risk = NA_real_, low = NA_real_, high = NA_real_))
    c(risk = 1 - sf$surv[j], low = 1 - sf$upper[j], high = 1 - sf$lower[j])
  }

  cat("Step 70: paired 20% capacity comparison using existing out-of-fold predictions\n")
  oof <- fread(inputs[["oof"]])
  if (anyDuplicated(oof, by = c("cohort", "landmark", "id")) || any(!is.finite(oof$state_risk_oof)) || any(!is.finite(oof$history_risk_oof))) stop("Invalid Step 66 predictions.")
  expected <- data.table(cohort = c("Mayo", "OASIS-3", "Mayo", "OASIS-3"), landmark = c("MRI2", "MRI2", "Second cognitive visit", "Second cognitive visit"), n = c(880L, 306L, 527L, 219L), e = c(168L, 59L, 76L, 44L), e3 = c(85L, 18L, 32L, 18L))
  actual <- oof[, .(n = .N, e = sum(event), e3 = sum(event == 1L & time <= horizon)), by = .(cohort, landmark)]
  check <- merge(expected, actual, by = c("cohort", "landmark"), suffixes = c("_expected", "_actual"), all = TRUE)
  if (nrow(check) != 4L || anyNA(check) || any(check$n_expected != check$n_actual | check$e_expected != check$e_actual | check$e3_expected != check$e3_actual)) stop("Step 66 evaluation populations changed.")
  joined <- merge(oof[landmark == "MRI2"], d[, .(cohort, id, expected_time = time, expected_event = event)], by = c("cohort", "id"), all = TRUE)
  if (anyNA(joined[, .(time, event, expected_time, expected_event)]) || any(abs(joined$time - joined$expected_time) > 1e-8 | joined$event != joined$expected_event)) stop("MRI2 out-of-fold outcomes disagree with the frozen risk set.")
  step_survival <- function(sf, times, left_limit = FALSE) {
    idx <- findInterval(times, sf$time)
    if (left_limit) {
      hit <- idx > 0L & sf$time[pmax(idx, 1L)] == times
      idx[hit] <- idx[hit] - 1L
    }
    c(1, sf$surv)[idx + 1L]
  }
  select_capacity <- function(score, id) {
    selected <- rep(FALSE, length(score))
    k <- as.integer(round(screening_fraction * length(score)))
    selected[order(-score, as.character(id), seq_along(score))[seq_len(k)]] <- TRUE
    selected
  }
  screen_metrics <- function(x) {
    state <- select_capacity(x$state_risk_oof, x$id)
    history <- select_capacity(x$history_risk_oof, x$id)
    cases <- x$event == 1L & x$time <= horizon
    if (!any(cases)) stop("No three-year events in bootstrap sample.")
    cf <- survival::survfit(survival::Surv(time, 1L - event) ~ 1, data = as.data.frame(x))
    g <- step_survival(cf, x$time[cases], left_limit = TRUE)
    if (any(!is.finite(g) | g <= 0)) stop("Nonpositive censoring survival at an event.")
    w <- 1 / g
    r_state <- km3(x[state])[["risk"]]; r_history <- km3(x[history])[["risk"]]
    capture_state <- sum(w * state[cases]) / sum(w)
    capture_history <- sum(w * history[cases]) / sum(w)
    f <- mean(state)
    c(state_selected_risk_3y = r_state, history_selected_risk_3y = r_history, selected_risk_gain = r_history - r_state,
      state_event_capture_3y = capture_state, history_event_capture_3y = capture_history, event_capture_gain = capture_history - capture_state,
      state_events_per_1000_candidates = 1000 * f * r_state, history_events_per_1000_candidates = 1000 * f * r_history,
      events_gain_per_1000_candidates = 1000 * f * (r_history - r_state))
  }
  screens <- list(); members <- list(); reclassified <- list(); screen_boot <- list()
  for (k in seq_len(nrow(expected))) {
    cohort_name <- expected$cohort[k]; landmark_name <- expected$landmark[k]
    x <- copy(oof[which(oof[["cohort"]] == cohort_name & oof[["landmark"]] == landmark_name)])
    setorder(x, id)
    point <- screen_metrics(x)
    err <- rep("", screening_boot_reps)
    boot <- t(vapply(seq_len(screening_boot_reps), function(b) {
      ix <- sample.int(nrow(x), nrow(x), replace = TRUE)
      tryCatch(screen_metrics(x[ix]), error = function(e) {err[b] <<- conditionMessage(e); rep(NA_real_, length(point))})
    }, numeric(length(point))))
    colnames(boot) <- names(point)
    state <- select_capacity(x$state_risk_oof, x$id); history <- select_capacity(x$history_risk_oof, x$id)
    members[[k]] <- copy(x)[, `:=`(selected_state = state, selected_history = history, capacity_fraction = screening_fraction)]
    screens[[k]] <- rbindlist(lapply(seq_along(point), function(j) {
      limits <- ci(boot[, j])
      data.table(cohort = cohort_name, landmark = landmark_name, analysis_role = if (landmark_name == "MRI2") "Primary_capacity_scenario" else "Secondary_reassessment_scenario",
        participants = nrow(x), selected_per_rule = sum(state), actual_selected_fraction = mean(state), events_total = sum(x$event),
        events_by_3y = sum(x$event == 1L & x$time <= horizon), metric = names(point)[j], estimate = point[j], CI_low = limits[1L], CI_high = limits[2L],
        bootstrap_requested = screening_boot_reps, bootstrap_valid = sum(is.finite(boot[, j])), CI_status = if (all(is.finite(limits))) "AVAILABLE" else "BOOTSTRAP_STABILITY_REVIEW")
    }))
    group <- fifelse(state & history, "Selected_by_both", fifelse(state, "Selected_by_state_only", fifelse(history, "Selected_by_history_only", "Selected_by_neither")))
    reclassified[[k]] <- copy(x)[, group := group][, .(participants = .N, observed_events_by_3y = sum(event == 1L & time <= horizon),
      known_3y_status = sum(time >= horizon | (event == 1L & time <= horizon))), by = .(cohort, landmark, group)]
    screen_boot[[k]] <- cbind(data.table(cohort = cohort_name, landmark = landmark_name, replicate = seq_len(screening_boot_reps)), as.data.table(boot))
    if (any(nzchar(err))) boot_errors[[length(boot_errors) + 1L]] <- data.table(cohort = cohort_name, comparison = paste("capacity", landmark_name), replicate = which(nzchar(err)), error = err[nzchar(err)])
  }
  save_tsv(rbindlist(screens), "Table6_paired_fixed_20pct_screening.tsv")
  save_tsv(rbindlist(reclassified), "Table7_screening_overlap_and_observed_events.tsv")
  save_tsv(rbindlist(members), "Frozen_OOF_screening_membership.tsv", "source_data")
  save_tsv(rbindlist(screen_boot), "Paired_screening_bootstrap_draws.tsv", "source_data")
  save_tsv(fread(inputs[["increment66"]]), "Table8_existing_Step66_AUC_Brier_net_benefit.tsv")

  cat("Step 70: same-day acquisition repeats at the frozen OASIS MRI anchors\n")
  muse <- fread(inputs[["muse"]])
  if (anyDuplicated(muse$Code)) stop("Duplicate MUSE Code.")
  muse[, `:=`(id = paste0("OAS", sub(".*OAS-[0-9]-([0-9]{5})-.*", "\\1", as.character(Code))),
    mri_day = as.integer(as.Date(EXAMDATE) - as.Date("2000-01-01")), hippo = number(Left_Hippocampus) + number(Right_Hippocampus), icv = number(ICV))]
  muse <- muse[id %in% risk[cohort == "OASIS-3", id]]
  digits <- sub(".*-", "", muse$Code)
  muse[, run := number(substr(digits, nchar(digits) - 2L, nchar(digits) - 1L))]
  coded_day <- number(substr(digits, 1L, nchar(digits) - 3L))
  muse[, valid_code_day := is.finite(coded_day) & coded_day == mri_day & grepl("^[0-9]+0$", digits)]
  meta <- fread(inputs[["mr_metadata"]])
  signature_fields <- c("DeviceSerialNumber", "ManufacturersModelName", "MagneticFieldStrength", "SoftwareVersions", "SeriesDescription", "RepetitionTime", "EchoTime", "FlipAngle", "SliceThickness")
  if (!all(c("subject_id", "label", "filename", signature_fields) %in% names(meta))) stop("Missing MRI metadata matching fields.")
  meta <- meta[grepl("_T1w\\.json$", filename) & subject_id %in% muse$id]
  meta[, `:=`(id = as.character(subject_id), mri_day = number(sub(".*_d([0-9]+).*", "\\1", label)),
    run = fifelse(grepl("_run-[0-9]+_", filename), number(sub(".*_run-([0-9]+)_.*", "\\1", filename)), 0))]
  for (field_name in signature_fields) set(meta, j = field_name, value = trimws(as.character(meta[[field_name]])))
  meta[, metadata_complete := rowSums(vapply(.SD, function(v) !is.na(v) & nzchar(v) & v != "NA", logical(nrow(meta)))) == length(signature_fields), .SDcols = signature_fields]
  meta[, signature := do.call(paste, c(.SD, sep = "|")), .SDcols = signature_fields]
  metadata_match <- meta[, .(metadata_rows = .N, distinct_signatures = uniqueN(signature), metadata_complete = all(metadata_complete),
    distinct_filenames = uniqueN(filename), signature = signature[1L], filename = as.character(filename[1L])), by = .(id, mri_day, run)]
  metadata_match[, exact_metadata := metadata_complete & distinct_signatures == 1L & distinct_filenames == 1L]
  scans <- merge(muse[, .(id, Code, mri_day, run, valid_code_day, hippo, icv)], metadata_match, by = c("id", "mri_day", "run"), all.x = TRUE)
  anchors <- fread(inputs[["scan_identity"]])
  scanner <- fread(inputs[["scanner"]])
  anchor_files <- rbindlist(list(scanner[, .(id = OASISID, scan = "MRI1", expected_filename = MRI1_filename)],
    scanner[, .(id = OASISID, scan = "MRI2", expected_filename = MRI2_filename)]))
  anchors <- merge(anchors, anchor_files, by = c("id", "scan"), all.x = TRUE)
  if (nrow(anchors) != 612L || anyDuplicated(anchors, by = c("id", "scan"))) stop("Frozen OASIS anchors changed.")
  repeat_pairs <- rbindlist(lapply(seq_len(nrow(anchors)), function(j) {
    anchor <- anchors[j]
    a <- scans[Code == anchor$Code]
    if (nrow(a) != 1L || a$mri_day != anchor$anchor_day) stop("Anchor scan identity mismatch: ", anchor$Code)
    candidates <- scans[id == anchor$id & mri_day == anchor$anchor_day & Code != anchor$Code]
    anchor_resolved <- isTRUE(a$valid_code_day) && isTRUE(a$exact_metadata) && isTRUE(a$filename == anchor$expected_filename)
    status <- if (!anchor_resolved) "ANCHOR_METADATA_UNRESOLVED" else "NO_MATCHED_ALTERNATE"
    b <- candidates[valid_code_day %in% TRUE & exact_metadata %in% TRUE & signature == a$signature]
    setorder(b, run, Code)
    if (nrow(b)) {
      b <- b[1L]
      status <- if (all(is.finite(c(a$hippo, a$icv, b$hippo, b$icv))) && min(c(a$hippo, a$icv, b$hippo, b$icv)) > 0 && anchor_resolved) "MATCHED_REPEAT" else "INVALID_VOLUME_OR_ANCHOR_METADATA"
    }
    valid <- identical(status, "MATCHED_REPEAT")
    data.table(id = anchor$id, scan = anchor$scan, mri_day = anchor$anchor_day, frozen_Code = anchor$Code,
      alternate_Code = if (nrow(b)) b$Code else NA_character_, same_day_alternates = nrow(candidates), status,
      frozen_hippo = a$hippo, frozen_ICV = a$icv, alternate_hippo = if (valid) b$hippo else NA_real_, alternate_ICV = if (valid) b$icv else NA_real_,
      raw_hippo_log_difference_pct = if (valid) 100 * log(b$hippo / a$hippo) else NA_real_,
      normalized_hippo_log_difference_pct = if (valid) 100 * log((b$hippo / b$icv) / (a$hippo / a$icv)) else NA_real_)
  }))
  save_tsv(scans, "OASIS_repeat_scan_metadata_matching.tsv", "audit")
  save_tsv(repeat_pairs, "OASIS_frozen_anchor_repeat_pairs.tsv", "source_data")
  save_tsv(repeat_pairs[, .(anchor_scans = .N, participants = uniqueN(id)), by = .(scan, status)], "Table9_repeat_acquisition_matching_counts.tsv")
  repeat_stat <- function(diff) {
    if (length(diff) < 2L) return(setNames(rep(NA_real_, 6L), c("mean_difference", "SD_difference", "lower_agreement_limit", "upper_agreement_limit", "within_acquisition_SD", "median_absolute_difference")))
    s <- sd(diff); m <- mean(diff)
    c(mean_difference = m, SD_difference = s, lower_agreement_limit = m - 1.96 * s,
      upper_agreement_limit = m + 1.96 * s, within_acquisition_SD = s / sqrt(2), median_absolute_difference = median(abs(diff)))
  }
  repeat_summaries <- list()
  for (field_name in c("raw_hippo_log_difference_pct", "normalized_hippo_log_difference_pct")) {
    z <- repeat_pairs[status == "MATCHED_REPEAT" & is.finite(get(field_name))]
    ids <- unique(z$id)
    point <- repeat_stat(z[[field_name]])
    if (length(ids) >= 2L) {
      blocks <- split(seq_len(nrow(z)), z$id)
      boot <- t(replicate(model_boot_reps, {
        sample_ids <- sample(ids, length(ids), replace = TRUE)
        ix <- unlist(blocks[sample_ids], use.names = FALSE)
        repeat_stat(z[[field_name]][ix])
      }))
    } else boot <- matrix(NA_real_, nrow = model_boot_reps, ncol = length(point))
    repeat_summaries[[field_name]] <- rbindlist(lapply(seq_along(point), function(j) {
      limits <- ci(boot[, j])
      data.table(measure = field_name, units = "100_times_log_ratio_percent", participants = length(ids), anchor_pairs = nrow(z), metric = names(point)[j],
        estimate = point[j], CI_low = limits[1L], CI_high = limits[2L], bootstrap_unit = "participant_preserving_both_anchors")
    }))
  }
  save_tsv(rbindlist(repeat_summaries), "Table10_same_day_volume_repeatability.tsv")
  both <- if (any(repeat_pairs$status == "MATCHED_REPEAT"))
    dcast(repeat_pairs[status == "MATCHED_REPEAT"], id ~ scan, value.var = c("alternate_hippo", "alternate_ICV")) else data.table()
  history_agreement <- data.table(status = "NO_PARTICIPANTS_WITH_MATCHED_REPEATS_AT_BOTH_ANCHORS")
  if (all(c("alternate_hippo_MRI1", "alternate_hippo_MRI2", "alternate_ICV_MRI2") %in% names(both))) {
    both <- both[is.finite(alternate_hippo_MRI1) & is.finite(alternate_hippo_MRI2) & is.finite(alternate_ICV_MRI2)]
    repeated <- merge(d[cohort == "OASIS-3"], both, by = "id")
    repeated[, `:=`(alternate_loss_pct_year = 100 * log(alternate_hippo_MRI1 / alternate_hippo_MRI2) / interval_years,
      alternate_hippo2_pct_ICV = 100 * alternate_hippo_MRI2 / alternate_ICV_MRI2)]
    save_tsv(repeated, "OASIS_two_anchor_alternate_history.tsv", "source_data")
    if (nrow(repeated) >= 3L) {
      agreement <- function(z) {
        a <- z$hippo_loss_pct_year; b <- z$alternate_loss_pct_year
        c(repeat_stat(b - a), Pearson_r = if (sd(a) > 0 && sd(b) > 0) cor(a, b) else NA_real_,
          concordance_correlation = 2 * cov(a, b) / (var(a) + var(b) + (mean(a) - mean(b))^2))
      }
      point <- agreement(repeated)
      boot <- t(replicate(model_boot_reps, agreement(repeated[sample.int(nrow(repeated), nrow(repeated), replace = TRUE)])))
      history_agreement <- rbindlist(lapply(seq_along(point), function(j) {
        limits <- ci(boot[, j])
        data.table(status = "TWO_ALTERNATE_ANCHORS", participants = nrow(repeated), events_total = sum(repeated$event), metric = names(point)[j],
          estimate = point[j], CI_low = limits[1L], CI_high = limits[2L],
          units = if (names(point)[j] %in% c("Pearson_r", "concordance_correlation")) "dimensionless" else "percentage_points_per_year")
      }))
      xr <- repeated[apoe_e4 %in% c(0, 1)]
      alternative <- copy(xr)
      alternative[, `:=`(hippo_loss_pct_year = alternate_loss_pct_year, hippo2_pct_ICV = alternate_hippo2_pct_ICV)]
      rr <- tryCatch(pair_models(xr, "OASIS-3", "Frozen_pair_vs_two_alternate_acquisitions", apoe_formula, apoe_formula, alternate = alternative), error = function(e) e)
      if (inherits(rr, "error")) {
        save_tsv(data.table(status = "MODEL_NOT_ESTIMABLE", reason = conditionMessage(rr), participants = nrow(xr), events_total = sum(xr$event)), "Table12_repeat_history_clinical_comparison.tsv")
      } else {
        save_tsv(rr$summary, "Table12_repeat_history_clinical_comparison.tsv")
        save_tsv(rr$coefficients, "Table13_repeat_history_Cox_coefficients.tsv")
      }
    }
  }
  save_tsv(history_agreement, "Table11_two_anchor_history_agreement.tsv")
  if (!file.exists(file.path(out, "tables/Table12_repeat_history_clinical_comparison.tsv"))) save_tsv(data.table(status = "NO_ESTIMABLE_TWO_ANCHOR_REPEAT_SUBSET"), "Table12_repeat_history_clinical_comparison.tsv")
  if (length(boot_errors)) save_tsv(rbindlist(boot_errors), "Bootstrap_failures.tsv", "audit")
  capture.output(sessionInfo(), file = file.path(out, "audit/R_session_info.txt"))
  writeLines(c("STEP_70_COMPLETED", normalizePath(out, winslash = "/", mustWork = TRUE)), file.path(out, "70_COMPLETED.txt"), useBytes = TRUE)
  cat("Step 70 complete: ", out, "\n", sep = "")
  print(fread(file.path(out, "tables/Table2_same_sample_standardized_risk_comparisons.tsv"))[metric == "risk_difference",
    .(cohort, comparison, participants, events_total, base, extended, extended_CI_low, extended_CI_high, CI_status)])
  print(fread(file.path(out, "tables/Table6_paired_fixed_20pct_screening.tsv"))[metric %in% c("selected_risk_gain", "event_capture_gain", "events_gain_per_1000_candidates"),
    .(cohort, landmark, selected_per_rule, metric, estimate, CI_low, CI_high, CI_status)])
})
