local({
  packages <- c("data.table", "lme4", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))
  library(data.table)

  root <- .study$output
  paths <- c(
    clinical = file.path(root, "cn_progression_62/run_frozen/tables/Table1_person_level_clinical_course.tsv"),
    later_mri = file.path(root, "cn_progression_62/run_frozen/tables/Table2_MRI2_and_later_structural_scans.tsv"),
    cognition = file.path(root, "cognitive_trajectory_63/run_frozen/tables/Table2_person_visit_cognitive_trajectories.tsv"),
    quartiles = file.path(root, "absolute_risk_65/run_frozen/tables/Table1_frozen_prior_loss_quartiles.tsv"),
    mayo_risk = file.path(root, "mcsa_event_43/run_frozen/tables/Table2_person_level_frozen_risk_sets.tsv"),
    mayo_zip = file.path(.study$data,"mayo","Demographics_and_Assessments.zip"),
    oasis_risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    oasis_visits = file.path(.study$output,"44_OASIS3_full_clinical_link_results/07_D1_all_visit_stages.csv"),
    oasis_scanner = file.path(root, "oasis_event_40/run_frozen/source_data/Table6_all_312_exact_scan_crosswalk.tsv")
  )
  if (any(!file.exists(paths)))
    stop("Missing frozen input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "aligned_longitudinal_68",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  save_tsv(data.table(input = names(paths), file = unname(paths),
                      md5 = unname(tools::md5sum(paths))), "Input_manifest.tsv", "audit")

  clinical <- fread(paths[["clinical"]])
  if (nrow(clinical[cohort == "Mayo"]) != 925L ||
      nrow(clinical[cohort == "OASIS-3"]) != 312L ||
      anyDuplicated(clinical, by = c("cohort", "id")))
    stop("Step 62 clinical population changed.")
  eligible <- clinical[prior_MCI_record == 0L]
  if (nrow(eligible[cohort == "Mayo"]) != 925L ||
      nrow(eligible[cohort == "OASIS-3"]) != 306L ||
      eligible[cohort == "Mayo", sum(impairment_event)] != 177L ||
      eligible[cohort == "OASIS-3", sum(impairment_event)] != 59L)
    stop("The final first-impairment risk sets changed.")
  quartiles <- fread(paths[["quartiles"]])
  check_q <- eligible[, .(loss_Q1 = as.numeric(quantile(hippo_loss_pct_year, 0.25)),
                          loss_Q3 = as.numeric(quantile(hippo_loss_pct_year, 0.75))), by = cohort]
  check_q <- merge(check_q, quartiles, by = "cohort", suffixes = c("_now", "_65"))
  if (nrow(check_q) != 2L ||
      any(abs(check_q$loss_Q1_now - check_q$loss_Q1_65) > 1e-8) ||
      any(abs(check_q$loss_Q3_now - check_q$loss_Q3_65) > 1e-8))
    stop("Step 65 prior-loss quartiles changed.")

  roles <- c("Primary_Mayo_global_z", "Mayo_memory_z",
             "OASIS_administered_MMSE")
  cognition <- fread(paths[["cognition"]])[role %in% roles]
  cognition <- merge(cognition,
                     eligible[, .(cohort, id, impairment_event,
                                  impairment_time_years)],
                     by = c("cohort", "id"))
  if (cognition[role == roles[1L], uniqueN(id)] != 860L ||
      cognition[role == roles[2L], uniqueN(id)] != 914L ||
      cognition[role == roles[3L], uniqueN(id)] != 305L ||
      anyDuplicated(cognition, by = c("cohort", "role", "id", "time_years")))
    stop("Step 63 records do not align with the final risk sets.")
  cognition[, window := "All_post_MRI2_5y"]
  pre_event <- copy(cognition[
    (impairment_event == 1L &
       time_years < impairment_time_years - 1 / 365.25) |
      (impairment_event == 0L &
         time_years <= impairment_time_years + 1 / 365.25)])
  pre_event[, window := "Before_first_impairment"]
  cognitive_sets <- rbindlist(list(cognition, pre_event))
  coverage <- cognitive_sets[, .(visits = .N,
                                 event = impairment_event[1L],
                                 last_year = max(time_years)),
                             by = .(cohort, role, window, id)][,
    .(participants = .N, future_visits = sum(visits),
      with_2plus_future_visits = sum(visits >= 2L),
      event_participants = sum(event),
      median_last_year = median(last_year)),
    by = .(cohort, role, window)]
  save_tsv(coverage, "Table1_aligned_cognitive_followup.tsv")

  fit_cognition <- function(x) {
    cohort_name <- x$cohort[1L]
    role_name <- x$role[1L]
    window_name <- x$window[1L]
    d <- as.data.frame(x)
    d$id <- factor(d$id)
    d$sex <- factor(d$sex)
    fit <- lme4::lmer(
      score ~ time_years * (hippo_loss_pct_year + hippo2_pct_ICV +
        baseline_score + age_MRI2 + sex + education + baseline_CDRSB) +
        (time_years | id), data = d, REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 200000)))
    term <- "time_years:hippo_loss_pct_year"
    beta <- unname(lme4::fixef(fit)[term])
    se <- sqrt(as.matrix(vcov(fit))[term, term])
    q <- quartiles[cohort == cohort_name]
    gap <- q$loss_Q3 - q$loss_Q1
    estimate <- data.table(
      cohort = cohort_name, role = role_name, window = window_name,
      participants = nlevels(d$id), observations = nrow(d),
      beta_per_1pct_year = beta, CI_low = beta - 1.96 * se,
      CI_high = beta + 1.96 * se,
      Wald_p = 2 * pnorm(-abs(beta / se)),
      loss_Q1 = q$loss_Q1, loss_Q3 = q$loss_Q3,
      Q3_minus_Q1_score_change_at_3y = 3 * gap * beta,
      contrast_CI_low = 3 * gap * (beta - 1.96 * se),
      contrast_CI_high = 3 * gap * (beta + 1.96 * se),
      singular = lme4::isSingular(fit))
    curve <- data.table(cohort = cohort_name, role = role_name,
                        window = window_name, years_after_MRI2 = 0:5)
    curve[, `:=`(Q3_minus_Q1_score_change = years_after_MRI2 * gap * beta,
                 CI_low = years_after_MRI2 * gap * (beta - 1.96 * se),
                 CI_high = years_after_MRI2 * gap * (beta + 1.96 * se))]
    list(estimate = estimate, curve = curve)
  }
  groups <- unique(cognitive_sets[, .(cohort, role, window)])
  fits <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i]
    fit_cognition(cognitive_sets[cohort == g$cohort & role == g$role &
                                   window == g$window])
  })
  cognitive_models <- rbindlist(lapply(fits, `[[`, "estimate"))
  save_tsv(cognitive_models, "Table2_aligned_cognitive_trajectories.tsv")
  save_tsv(rbindlist(lapply(fits, `[[`, "curve")),
           "Figure68_cognitive_trajectory_source.tsv", "source_data")

  later_mri <- fread(paths[["later_mri"]])
  later_mri <- merge(later_mri,
                     eligible[, .(cohort, id, impairment_event,
                                  impairment_time_years)],
                     by = c("cohort", "id"))
  later_mri <- later_mri[time_years > 0 &
    ((impairment_event == 1L &
        time_years < impairment_time_years - 1 / 365.25) |
       (impairment_event == 0L &
          time_years <= impairment_time_years + 1 / 365.25))]
  setorder(later_mri, cohort, id, time_years)
  repeat_ids <- later_mri[, .N, by = .(cohort, id)][N >= 2L,
                          .(cohort, id)]
  two_scans <- later_mri[repeat_ids, on = .(cohort, id)]
  setorder(two_scans, cohort, id, time_years)
  scan_coverage <- later_mri[, .(later_scans = .N,
                                 event = impairment_event[1L]),
                             by = .(cohort, id)][,
    .(participants_with_later_scan = .N,
      participants_with_2plus_later_scans = sum(later_scans >= 2L),
      later_scans = sum(later_scans),
      event_participants = sum(event)), by = cohort]
  save_tsv(scan_coverage, "Table3_post_MRI2_scan_coverage_before_first_event.tsv")

  fit_structure <- function(x, cohort_name) {
    d <- as.data.frame(x)
    d$id <- factor(d$id)
    d$sex <- factor(d$sex)
    fit <- lme4::lmer(
      log_hippo ~ time_years * (hippo_loss_pct_year + age_MRI2 +
        sex + education + baseline_CDRSB + hippo2_pct_ICV) + (1 | id),
      data = d, REML = FALSE)
    term <- "time_years:hippo_loss_pct_year"
    beta <- unname(lme4::fixef(fit)[term])
    se <- sqrt(as.matrix(vcov(fit))[term, term])
    data.table(cohort = cohort_name, participants = nlevels(d$id),
               later_scans = nrow(d),
               beta_per_1pct_year = beta, CI_low = beta - 1.96 * se,
               CI_high = beta + 1.96 * se,
               Wald_p = 2 * pnorm(-abs(beta / se)),
               singular = lme4::isSingular(fit))
  }
  structure_models <- rbindlist(lapply(c("Mayo", "OASIS-3"), function(name)
    fit_structure(two_scans[cohort == name], name)))
  save_tsv(structure_models,
           "Table4_future_hippocampal_change_without_shared_MRI2.tsv")
  future_rates <- two_scans[, .(
    later_scan_count = .N,
    later_scan_interval_years = time_years[.N] - time_years[1L],
    independent_future_change_pct_year =
      (log_hippo[.N] - log_hippo[1L]) /
      (time_years[.N] - time_years[1L]),
    prior_loss_pct_year = hippo_loss_pct_year[1L]),
    by = .(cohort, id)]
  if (any(!is.finite(future_rates$independent_future_change_pct_year)))
    stop("A later MRI pair has zero interval or invalid volume.")
  save_tsv(future_rates, "Table5_independent_future_scan_pairs.tsv")
  future_summary <- future_rates[, {
    test <- suppressWarnings(cor.test(prior_loss_pct_year,
      independent_future_change_pct_year, method = "spearman", exact = FALSE))
    .(participants = .N,
      median_later_scan_interval_years = median(later_scan_interval_years),
      spearman_rho = unname(test$estimate), spearman_p = test$p.value)
  }, by = cohort]
  save_tsv(future_summary, "Table5b_independent_future_scan_summary.tsv")

  mayo_risk <- fread(paths[["mayo_risk"]])[stage == "CN",
    .(id, MRI2_date = as.Date(MRI2_date), MRI_gap_years,
      same_manufacturer = as.integer(same_manufacturer))]
  oasis_risk <- fread(paths[["oasis_risk"]])[baseline_stage == "CN",
    .(id = OASISID, MRI2_day, MRI_span_years =
        (MRI2_day - MRI1_day) / 365.25)]
  oasis_scanner <- fread(paths[["oasis_scanner"]])
  if (nrow(oasis_scanner) != 312L || anyDuplicated(oasis_scanner$OASISID))
    stop("Step 40 OASIS-3 scan crosswalk changed.")
  mayo_protocol <- merge(eligible[cohort == "Mayo"], mayo_risk, by = "id")
  oasis_protocol <- merge(eligible[cohort == "OASIS-3"], oasis_scanner[
    , .(id = OASISID, same_scanner, same_scanner_series,
        both_exact, serial_known)], by = "id")
  oasis_protocol <- merge(oasis_protocol,
                          oasis_risk[, .(id, MRI_span_years)], by = "id")
  if (nrow(mayo_protocol) != 925L || nrow(oasis_protocol) != 306L)
    stop("Pre-landmark protocol metadata do not match the clinical risk sets.")
  protocol_coverage <- rbindlist(list(
    mayo_protocol[, .(cohort = "Mayo", comparison =
                        "Same manufacturer across MRI1/MRI2",
                      participants = .N, consistent = sum(same_manufacturer == 1L),
                      consistent_events = sum(impairment_event[same_manufacturer == 1L]),
                      median_mri_interval_years = median(MRI_gap_years))],
    oasis_protocol[, .(cohort = "OASIS-3", comparison =
                         "Same physical device and T1 series",
                       participants = .N, consistent = sum(same_scanner_series),
                       consistent_events = sum(impairment_event[same_scanner_series]),
                       median_mri_interval_years = median(MRI_span_years))]))
  save_tsv(protocol_coverage, "Table6_prelandmark_protocol_coverage.tsv")

  fit_event <- function(d, cohort_name, analysis) {
    x <- as.data.frame(d)
    x$sex <- factor(x$sex)
    fit <- survival::coxph(
      survival::Surv(impairment_time_years, impairment_event) ~
        age_MRI2 + sex + education + baseline_CDRSB +
        hippo2_pct_ICV + hippo_loss_pct_year,
      data = x, ties = "efron")
    term <- "hippo_loss_pct_year"
    beta <- unname(coef(fit)[term])
    se <- sqrt(vcov(fit)[term, term])
    data.table(cohort = cohort_name, analysis = analysis,
               participants = nrow(x), events = sum(x$impairment_event),
               HR_per_1pct_year = exp(beta),
               CI_low = exp(beta - 1.96 * se),
               CI_high = exp(beta + 1.96 * se),
               Wald_p = 2 * pnorm(-abs(beta / se)))
  }
  protocol_models <- rbindlist(list(
    fit_event(mayo_protocol, "Mayo", "All eligible"),
    fit_event(mayo_protocol[same_manufacturer == 1L], "Mayo",
              "Same manufacturer"),
    fit_event(oasis_protocol, "OASIS-3", "All eligible"),
    fit_event(oasis_protocol[same_scanner_series == TRUE], "OASIS-3",
              "Same device and T1 series")))
  reference_hr <- c(Mayo = 1.17390287870943, `OASIS-3` = 1.61323880367617)
  if (any(abs(protocol_models[analysis == "All eligible", HR_per_1pct_year] -
              reference_hr[protocol_models[analysis == "All eligible", cohort]]) > 1e-5))
    stop("The unrestricted Step 62 first-impairment model did not reproduce.")
  save_tsv(protocol_models, "Table7_protocol_restricted_first_event_models.tsv")

  mayo_raw <- as.data.table(read.csv(unz(paths[["mayo_zip"]],
                                     "MCSA_Data_12Nov2026.csv"),
                                   stringsAsFactors = FALSE, check.names = FALSE))
  mayo_raw <- mayo_raw[MCSA_ID %in% eligible[
    cohort == "Mayo" & impairment_event == 1L, id],
    .(id = MCSA_ID, visit_date = as.Date(VisitDate),
      stage_code = suppressWarnings(as.integer(FCOGIMPR)))]
  mayo_visits <- merge(mayo_raw, mayo_risk[, .(id, MRI2_date)], by = "id")
  mayo_visits[, day := as.numeric(visit_date - MRI2_date)]
  first_unknown <- mayo_visits[day > 0 & stage_code == 4L,
                               .(unknown_day = min(day)), by = id]
  mayo_visits <- merge(mayo_visits, first_unknown, by = "id", all.x = TRUE)
  mayo_visits <- mayo_visits[day > 0 &
    (is.na(unknown_day) | day < unknown_day) &
    stage_code %in% c(0L, 1L, 3L),
    .(cohort = "Mayo", id, day,
      stage = c("CN", "MCI", "Dementia")[match(stage_code,
                                               c(0L, 1L, 3L))])]
  oasis_visits <- merge(fread(paths[["oasis_visits"]])[
    OASISID %in% eligible[
      cohort == "OASIS-3" & impairment_event == 1L, id] &
      stage %in% c("CN", "MCI", "Dementia"),
    .(id = OASISID, day, stage)],
    oasis_risk[, .(id, MRI2_day)], by = "id")
  oasis_visits[, day := day - MRI2_day]
  oasis_visits <- oasis_visits[day > 0,
    .(cohort = "OASIS-3", id, day, stage)]
  visits <- unique(rbindlist(list(mayo_visits, oasis_visits)),
                   by = c("cohort", "id", "day", "stage"))
  if (anyDuplicated(visits, by = c("cohort", "id", "day")))
    stop("Conflicting clinical classifications occur on one visit date.")
  setkey(visits, cohort, id)
  event_people <- eligible[impairment_event == 1L,
    .(cohort, id, sequence, impairment_time_years)]
  event_gaps <- rbindlist(lapply(seq_len(nrow(event_people)), function(i) {
    person <- event_people[i]
    v <- visits[.(person$cohort, person$id)]
    setorder(v, day)
    impaired <- v[stage %in% c("MCI", "Dementia")]
    if (!nrow(impaired))
      stop("No first impairment visit: ", person$cohort, " / ", person$id)
    first <- impaired[1L]
    if (abs(first$day / 365.25 - person$impairment_time_years) > 1 / 365.25)
      stop("First impairment date changed: ", person$cohort, " / ", person$id)
    prior_cn <- v[stage == "CN" & day < first$day, day]
    last_cn <- if (length(prior_cn)) max(prior_cn) else 0
    data.table(cohort = person$cohort, id = person$id,
               first_event = first$stage,
               days_since_last_CN_or_MRI2 = first$day - last_cn,
               no_CN_visit_after_MRI2 = length(prior_cn) == 0L,
               years_from_MRI2_to_event = person$impairment_time_years)
  }))
  if (event_gaps[cohort == "Mayo" & first_event == "MCI", .N] != 173L ||
      event_gaps[cohort == "Mayo" & first_event == "Dementia", .N] != 4L ||
      event_gaps[cohort == "OASIS-3" & first_event == "MCI", .N] != 22L ||
      event_gaps[cohort == "OASIS-3" & first_event == "Dementia", .N] != 37L)
    stop("First-event composition changed.")
  save_tsv(event_gaps, "Table8_person_level_first_event_gap.tsv")
  gap_summary <- event_gaps[, .(
    events = .N,
    no_CN_visit_after_MRI2 = sum(no_CN_visit_after_MRI2),
    median_days_since_last_CN_or_MRI2 = median(days_since_last_CN_or_MRI2),
    Q1_days = as.numeric(quantile(days_since_last_CN_or_MRI2, 0.25)),
    Q3_days = as.numeric(quantile(days_since_last_CN_or_MRI2, 0.75)),
    median_years_from_MRI2_to_event = median(years_from_MRI2_to_event)),
    by = .(cohort, first_event)]
  save_tsv(gap_summary, "Table9_first_event_composition_and_visit_gap.tsv")

  writeLines(capture.output(sessionInfo()),
             file.path(out, "audit", "R_sessionInfo.txt"))
  writeLines("STATUS=STEP68_COMPLETE", file.path(out, "audit", "68_COMPLETED.txt"))
  message("Step 68 complete: ", out)
  print(coverage, row.names = FALSE)
  print(cognitive_models, row.names = FALSE)
  print(scan_coverage, row.names = FALSE)
  print(structure_models, row.names = FALSE)
  print(future_summary, row.names = FALSE)
  print(protocol_models, row.names = FALSE)
  print(gap_summary, row.names = FALSE)
})
