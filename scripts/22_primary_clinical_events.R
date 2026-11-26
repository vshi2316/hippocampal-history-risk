local({
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("lme4", quietly = TRUE))
    stop("Install data.table, survival, and lme4 before running this script.")
  library(data.table)

  root <- .study$output
  paths <- c(
    oasis_risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    oasis_visits = file.path(.study$output,"44_OASIS3_full_clinical_link_results/07_D1_all_visit_stages.csv"),
    oasis_scans = file.path(.study$output,"44_OASIS3_full_clinical_link_results/06_selected_MUSE_scans.csv"),
    mayo_risk = file.path(root, "mcsa_event_43/run_frozen/tables/Table2_person_level_frozen_risk_sets.tsv"),
    mayo_zip = file.path(.study$data,"mayo","Demographics_and_Assessments.zip")
  )
  if (any(!file.exists(paths)))
    stop("Missing frozen input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "cn_progression_62",
                   paste0("run_", "frozen"))
  for (d in c("tables", "audit"))
    dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  save_tsv(data.table(input = names(paths), file = unname(paths),
                      md5 = unname(tools::md5sum(paths))),
           "Input_manifest.tsv", "audit")

  first_day <- function(x) if (length(x)) min(x) else NA_real_
  summarize_course <- function(v, censor_day, expected_event) {
    v <- v[day > 0 & day <= censor_day & stage %in% c("CN", "MCI", "Dementia")]
    setorder(v, day)
    dementia_day <- first_day(v[stage == "Dementia", day])
    if (as.integer(!is.na(dementia_day)) != expected_event ||
        (expected_event == 1L && dementia_day != censor_day))
      stop("A dementia record disagrees with the frozen endpoint.")
    mci_day <- first_day(v[stage == "MCI" &
                             (is.na(dementia_day) | day < dementia_day), day])
    impairment_day <- first_day(v[stage %in% c("MCI", "Dementia"), day])
    reversal <- !is.na(mci_day) &&
      any(v$stage == "CN" & v$day > mci_day &
            (is.na(dementia_day) | v$day < dementia_day))
    sequence <- if (!is.na(mci_day) && !is.na(dementia_day)) {
      "MCI_then_dementia"
    } else if (!is.na(mci_day)) {
      "MCI_no_dementia_record"
    } else if (!is.na(dementia_day)) {
      "Dementia_without_prior_MCI_record"
    } else "No_impairment_record"
    data.table(
      impairment_event = as.integer(!is.na(impairment_day)),
      impairment_time_years = if (!is.na(impairment_day))
        impairment_day / 365.25 else censor_day / 365.25,
      first_MCI_years = mci_day / 365.25,
      first_dementia_years = dementia_day / 365.25,
      MCI_to_dementia_years = if (!is.na(mci_day) && !is.na(dementia_day))
        (dementia_day - mci_day) / 365.25 else NA_real_,
      MCI_then_CN_before_dementia = as.integer(reversal),
      sequence = sequence,
      post_MRI2_eligible_visits = nrow(v)
    )
  }

  oasis <- fread(paths[["oasis_risk"]])[baseline_stage == "CN"]
  oasis[, event := as.integer(as.logical(event))]
  if (nrow(oasis) != 312L || sum(oasis$event) != 52L ||
      anyDuplicated(oasis$OASISID)) stop("OASIS-3 frozen CN set changed.")
  ov <- fread(paths[["oasis_visits"]])
  if (!all(c("OASISID", "day", "stage") %in% names(ov)))
    stop("OASIS-3 clinical visit fields changed.")
  setkey(ov, OASISID)
  oasis_course <- rbindlist(lapply(seq_len(nrow(oasis)), function(i) {
    r <- oasis[i]
    v <- copy(ov[.(r$OASISID)])
    prior_mci <- as.integer(any(v$stage == "MCI" & v$day < r$MRI2_day,
                                na.rm = TRUE))
    v[, day := day - r$MRI2_day]
    course <- summarize_course(v, r$followup_days, r$event)
    if (r$event == 1L &&
        abs(course$first_dementia_years * 365.25 -
              (r$first_dementia_day - r$MRI2_day)) > 1e-6)
      stop("OASIS-3 first dementia date changed: ", r$OASISID)
    cbind(data.table(cohort = "OASIS-3", id = r$OASISID,
                     age_MRI2 = r$age_MRI2, sex = r$sex,
                     education = r$education,
                     baseline_CDRSB = r$baseline_CDRSB,
                     hippo2_pct_ICV = r$hippo2_pct_ICV,
                     hippo_loss_pct_year = r$hippo_loss_pct_year,
                     prior_MCI_record = prior_mci,
                     dementia_event = r$event,
                     dementia_time_years = r$followup_days / 365.25),
          course)
  }))

  mayo <- fread(paths[["mayo_risk"]])[stage == "CN"]
  mayo[, event := as.integer(event)]
  if (nrow(mayo) != 925L || sum(mayo$event) != 21L ||
      anyDuplicated(mayo$id)) stop("Mayo frozen CN set changed.")
  raw <- as.data.table(read.csv(unz(paths[["mayo_zip"]],
                                  "MCSA_Data_12Nov2026.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE))
  if (!all(c("MCSA_ID", "VisitDate", "FCOGIMPR") %in% names(raw)))
    stop("Mayo clinical visit fields changed.")
  raw[, `:=`(VisitDate = as.Date(VisitDate),
             FCOGIMPR = suppressWarnings(as.numeric(FCOGIMPR)))]
  setkey(raw, MCSA_ID)
  mayo_course <- rbindlist(lapply(seq_len(nrow(mayo)), function(i) {
    r <- mayo[i]
    d <- raw[.(r$id)]
    index_date <- as.Date(r$MRI2_date)
    unknown <- d[!is.na(VisitDate) & VisitDate > index_date &
                   FCOGIMPR == 4, VisitDate]
    first_unknown <- if (length(unknown)) min(unknown) else as.Date(NA)
    future <- d[!is.na(VisitDate) & VisitDate > index_date &
                  FCOGIMPR %in% c(0, 1, 3) &
                  (is.na(first_unknown) | VisitDate < first_unknown)]
    if (!nrow(future)) stop("Frozen Mayo participant has no follow-up: ", r$id)
    setorder(future, VisitDate)
    first_dementia <- future[FCOGIMPR == 3, VisitDate]
    event <- as.integer(length(first_dementia) > 0L)
    end_date <- if (event == 1L) first_dementia[1L] else max(future$VisitDate)
    censor_day <- as.numeric(end_date - index_date)
    if (event != r$event ||
        abs(censor_day / 365.25 - r$followup_years) > 1e-6)
      stop("Mayo dementia outcome or follow-up disagrees with frozen risk set: ", r$id)
    v <- future[, .(day = as.numeric(VisitDate - index_date),
                    stage = c("CN", "MCI", "Dementia")[match(FCOGIMPR,
                                                             c(0, 1, 3))])]
    course <- summarize_course(v, censor_day, event)
    cbind(data.table(cohort = "Mayo", id = r$id,
                     age_MRI2 = r$age_MRI2, sex = r$male,
                     education = r$education,
                     baseline_CDRSB = r$baseline_CDRSB,
                     hippo2_pct_ICV = r$hippo2_pct_ICV,
                     hippo_loss_pct_year = r$hippo_loss_pct_year,
                     prior_MCI_record = 0L,
                     dementia_event = event,
                     dementia_time_years = r$followup_years),
          course)
  }))
  people <- rbindlist(list(oasis_course, mayo_course))
  if (nrow(people) != 1237L || sum(people$dementia_event) != 73L ||
      any(!is.finite(people$impairment_time_years)) ||
      any(people$impairment_time_years <= 0))
    stop("Combined frozen risk set or impairment times changed.")
  save_tsv(people, "Table1_person_level_clinical_course.tsv")

  oasis_scans <- fread(paths[["oasis_scans"]])
  if (!all(c("OASISID", "mri_day", "Code", "ICV", "hippocampus") %in%
           names(oasis_scans))) stop("OASIS-3 selected scan fields changed.")
  setkey(oasis_scans, OASISID)
  os <- rbindlist(lapply(seq_len(nrow(oasis)), function(i) {
    r <- oasis[i]
    s <- copy(oasis_scans[.(r$OASISID)])
    baseline <- s[Code == r$MRI2_Code & mri_day == r$MRI2_day]
    if (nrow(baseline) != 1L ||
        abs(100 * baseline$hippocampus / baseline$ICV -
              r$hippo2_pct_ICV) > 1e-6)
      stop("OASIS-3 MRI2 measurement disagrees with frozen risk set: ", r$OASISID)
    s <- s[mri_day == r$MRI2_day & Code == r$MRI2_Code |
             mri_day >= r$MRI2_day + 180L &
             mri_day <= r$MRI2_day + 1827L]
    s <- s[is.finite(hippocampus) & hippocampus > 0 &
             is.finite(ICV) & ICV > 0]
    setorder(s, mri_day, Code)
    s <- unique(s, by = "mri_day")
    if (nrow(s) < 2L) return(NULL)
    data.table(cohort = "OASIS-3", id = r$OASISID,
               time_years = (s$mri_day - r$MRI2_day) / 365.25,
               log_hippo = 100 * log(s$hippocampus / s$ICV),
               post_dementia = r$event == 1L &
                 s$mri_day >= r$first_dementia_day,
               within_clinical_followup =
                 s$mri_day - r$MRI2_day <= r$followup_days,
               age_MRI2 = r$age_MRI2, sex = r$sex,
               education = r$education, baseline_CDRSB = r$baseline_CDRSB,
               hippo2_pct_ICV = r$hippo2_pct_ICV,
               hippo_loss_pct_year = r$hippo_loss_pct_year)
  }), fill = TRUE)

  raw[, `:=`(ImagingDate = as.Date(ImagingDate),
             HPVOL = suppressWarnings(as.numeric(HPVOL)),
             TIV = suppressWarnings(as.numeric(TIV)),
             calc_age_vis = suppressWarnings(as.numeric(calc_age_vis)),
             Male = suppressWarnings(as.numeric(Male)),
             EDUC = suppressWarnings(as.numeric(EDUC)),
             CDRSUM = suppressWarnings(as.numeric(CDRSUM)))]
  raw[, imaging_gap_days := as.numeric(ImagingDate - VisitDate)]
  ms <- rbindlist(lapply(seq_len(nrow(mayo)), function(i) {
    r <- mayo[i]
    s <- copy(raw[.(r$id)])
    index_date <- as.Date(r$MRI2_date)
    s <- s[!is.na(ImagingDate) & !is.na(VisitDate) &
             imaging_gap_days >= 0 & imaging_gap_days <= 90 &
             is.finite(HPVOL) & HPVOL > 0 & is.finite(TIV) & TIV > 0 &
             is.finite(calc_age_vis) & is.finite(Male) &
             is.finite(EDUC) & is.finite(CDRSUM) &
             (ImagingDate == index_date & FCOGIMPR == 0 |
                ImagingDate >= index_date + 180 &
                ImagingDate <= index_date + 1827)]
    setorder(s, ImagingDate, VisitDate, visit_num)
    baseline <- s[ImagingDate == index_date][1L]
    if (nrow(baseline) != 1L ||
        abs(100 * baseline$HPVOL / baseline$TIV - r$hippo2_pct_ICV) > 1e-6)
      stop("Mayo MRI2 measurement disagrees with frozen risk set: ", r$id)
    s <- rbindlist(list(baseline,
                        unique(s[ImagingDate > index_date], by = "ImagingDate")))
    if (nrow(s) < 2L) return(NULL)
    data.table(cohort = "Mayo", id = r$id,
               time_years = as.numeric(s$ImagingDate - index_date) / 365.25,
               log_hippo = 100 * log(s$HPVOL / s$TIV),
               post_dementia = r$event == 1L &
                 as.numeric(s$ImagingDate - index_date) / 365.25 >=
                 r$followup_years,
               within_clinical_followup =
                 as.numeric(s$ImagingDate - index_date) / 365.25 <=
                 r$followup_years + 1e-6,
               age_MRI2 = r$age_MRI2, sex = r$male,
               education = r$education, baseline_CDRSB = r$baseline_CDRSB,
               hippo2_pct_ICV = r$hippo2_pct_ICV,
               hippo_loss_pct_year = r$hippo_loss_pct_year)
  }), fill = TRUE)
  scans <- rbindlist(list(os, ms), fill = TRUE)
  save_tsv(scans, "Table2_MRI2_and_later_structural_scans.tsv")
  scan_flow <- scans[, .(participants = uniqueN(id), scans = .N,
                         later_scans = sum(time_years > 0),
                         participants_with_3plus_scans =
                           sum(tabulate(match(id, unique(id))) >= 3L),
                         post_dementia_scans = sum(post_dementia)),
                     by = cohort]
  save_tsv(scan_flow, "Table3_postlandmark_scan_coverage.tsv")

  sequence_counts <- people[, .(participants = .N,
                                median_first_impairment_years =
                                  if (all(is.na(first_MCI_years)) &&
                                      all(is.na(first_dementia_years))) NA_real_
                                  else median(impairment_time_years[
                                    impairment_event == 1L])),
                            by = .(cohort, prior_MCI_record, sequence)][
                              order(cohort, prior_MCI_record, sequence)]
  save_tsv(sequence_counts, "Table5_observed_stage_sequences.tsv")
  lead <- people[sequence == "MCI_then_dementia",
                 .(participants = .N,
                   median_MCI_to_dementia_years = median(MCI_to_dementia_years),
                   Q1_years = unname(quantile(MCI_to_dementia_years, 0.25)),
                   Q3_years = unname(quantile(MCI_to_dementia_years, 0.75)),
                   MCI_then_CN_before_dementia = sum(MCI_then_CN_before_dementia)),
                 by = cohort]
  save_tsv(lead, "Table6_recorded_MCI_to_dementia_intervals.tsv")

  fit_endpoint <- function(d, endpoint) {
    x <- as.data.frame(copy(d))
    if (endpoint == "First_impairment")
      x <- x[x$prior_MCI_record == 0L, , drop = FALSE]
    x$time <- if (endpoint == "First_impairment")
      x$impairment_time_years else x$dementia_time_years
    x$status <- if (endpoint == "First_impairment")
      x$impairment_event else x$dementia_event
    covars <- c("time", "status", "age_MRI2", "sex", "education",
                "baseline_CDRSB", "hippo2_pct_ICV", "hippo_loss_pct_year")
    if (any(!complete.cases(x[covars])))
      stop("Missing model fields in ", d$cohort[1L], ": ", endpoint)
    if (d$cohort[1L] == "OASIS-3") x$sex <- factor(x$sex)
    base <- survival::coxph(
      survival::Surv(time, status) ~ age_MRI2 + sex + education +
        baseline_CDRSB + hippo2_pct_ICV, data = x, ties = "efron")
    plus <- survival::coxph(
      survival::Surv(time, status) ~ age_MRI2 + sex + education +
        baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year,
      data = x, ties = "efron")
    beta <- unname(coef(plus)["hippo_loss_pct_year"])
    se <- sqrt(vcov(plus)["hippo_loss_pct_year", "hippo_loss_pct_year"])
    lr <- 2 * (as.numeric(logLik(plus)) - as.numeric(logLik(base)))
    ph_p <- survival::cox.zph(plus)$table["hippo_loss_pct_year", "p"]
    data.table(cohort = d$cohort[1L], endpoint = endpoint,
               participants = nrow(x), events = sum(x$status),
               HR_per_1pct_year = exp(beta), CI_low = exp(beta - 1.96 * se),
               CI_high = exp(beta + 1.96 * se),
               Wald_p = 2 * pnorm(-abs(beta / se)),
               LR_chisq = lr, LR_p = pchisq(lr, 1, lower.tail = FALSE),
               PH_p_prior_loss = ph_p,
               apparent_C_base = unname(summary(base)$concordance[1L]),
               apparent_C_plus_loss = unname(summary(plus)$concordance[1L]))
  }
  models <- rbindlist(lapply(split(people, people$cohort), function(d)
    rbindlist(lapply(c("First_impairment", "Dementia"),
                     function(endpoint) fit_endpoint(d, endpoint)))))
  check_hr <- models[endpoint == "Dementia"]
  if (abs(check_hr[cohort == "OASIS-3", HR_per_1pct_year] - 1.543) > 0.005 ||
      abs(check_hr[cohort == "Mayo", HR_per_1pct_year] - 1.584) > 0.005)
    stop("The dementia Cox estimate no longer reproduces the manuscript.")
  save_tsv(models, "Table7_first_impairment_and_dementia_models.tsv")

  horizon <- rbindlist(lapply(c("impairment", "dementia"), function(endpoint)
    rbindlist(lapply(c(3, 5, 8), function(years) {
      time <- paste0(endpoint, "_time_years")
      event <- paste0(endpoint, "_event")
      tab <- people[, .(participants = .N,
                        known_status = sum(get(time) >= years |
                          (get(event) == 1L & get(time) <= years)),
                        events_by_horizon = sum(get(event) == 1L &
                                                  get(time) <= years)),
                    by = cohort]
      tab[, `:=`(endpoint = endpoint, horizon_years = years)]
      tab
    }))))
  save_tsv(horizon, "Table8_followup_horizon_information.tsv")

writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
})
