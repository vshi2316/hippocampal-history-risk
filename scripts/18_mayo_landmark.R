local({
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Install the survival package before running this script.")

  zip_path <- file.path(.study$data,"mayo","Demographics_and_Assessments.zip")
  oasis_path <- paste0(
    file.path(.study$output,"oasis_event_38/"),
    "run_frozen/tables/",
    "Table4_primary_increment_and_frozen_reproduction.tsv"
  )
  output_root <- file.path(.study$output,"mcsa_event_43")
  output_dir <- file.path(output_root, paste0("run_", "frozen"))
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

  save_tsv <- function(x, name) {
    write.table(x, file.path(output_dir, "tables", name), sep = "\t",
                row.names = FALSE, quote = FALSE, na = "")
  }
  valid <- function(x) is.finite(x) & !is.na(x)

  message("Reading MCSA longitudinal assessments and imaging summaries")
  raw <- read.csv(unz(zip_path, "MCSA_Data_12Nov2026.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
  dict <- read.csv(unz(zip_path, "MCSA_Data_Dictionary_12Nov2026.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("MCSA_ID", "VisitDate", "ImagingDate", "FCOGIMPR",
              "HPVOL", "TIV", "calc_age_vis", "Male", "EDUC", "CDRSUM",
              "MRIMANUFACTURER")
  if (!all(needed %in% names(raw))) stop("MCSA input is missing a required field.")
  raw$VisitDate <- as.Date(raw$VisitDate)
  raw$ImagingDate <- as.Date(raw$ImagingDate)
  for (v in c("FCOGIMPR", "HPVOL", "TIV", "calc_age_vis", "Male",
              "EDUC", "CDRSUM")) raw[[v]] <- as.numeric(raw[[v]])
  if (!all(c("HPVOL", "TIV", "FCOGIMPR") %in% dict$VARIABLE_NAME))
    stop("The MCSA dictionary does not match the analysis fields.")

  people <- split(raw, raw$MCSA_ID, drop = TRUE)
  message("Constructing MRI2 landmark risk sets")
  build_risk_set <- function(stage) {
    steps <- c(all_people = length(people), two_eligible_scans = 0L,
               no_prior_impairment_or_dementia = 0L,
               later_clinical_visit = 0L)
    records <- vector("list", length(people))
    k <- 0L
    for (d in people) {
      gap <- as.numeric(d$ImagingDate - d$VisitDate)
      ok <- d$FCOGIMPR == stage & valid(d$HPVOL) & d$HPVOL > 0 &
        valid(d$TIV) & d$TIV > 0 & valid(d$calc_age_vis) &
        valid(d$Male) & valid(d$EDUC) & valid(d$CDRSUM) &
        !is.na(d$VisitDate) & !is.na(d$ImagingDate) &
        !is.na(gap) & gap >= 0 & gap <= 90
      scans <- d[which(ok), , drop = FALSE]
      if (nrow(scans) < 2L) next
      scans <- scans[order(scans$ImagingDate, scans$VisitDate,
                           scans$visit_num), , drop = FALSE]
      chosen <- NULL
      for (j in seq.int(2L, nrow(scans))) {
        intervals <- as.numeric(scans$ImagingDate[j] -
                                  scans$ImagingDate[seq_len(j - 1L)])
        earlier <- which(intervals >= 365 & intervals <= 1827)
        if (length(earlier)) {
          i <- max(earlier)
          chosen <- list(first = scans[i, , drop = FALSE],
                         second = scans[j, , drop = FALSE],
                         days = intervals[i])
          break
        }
      }
      if (is.null(chosen)) next
      steps["two_eligible_scans"] <- steps["two_eligible_scans"] + 1L
      first <- chosen$first
      second <- chosen$second
      index_date <- second$ImagingDate
      previous <- d[!is.na(d$VisitDate) & d$VisitDate < index_date,
                    , drop = FALSE]
      excluded_codes <- if (stage == 0) c(1, 3, 4) else c(3, 4)
      if (any(previous$FCOGIMPR %in% excluded_codes, na.rm = TRUE)) next
      steps["no_prior_impairment_or_dementia"] <-
        steps["no_prior_impairment_or_dementia"] + 1L

      unknown <- d$VisitDate[!is.na(d$VisitDate) &
                               d$VisitDate > index_date & d$FCOGIMPR %in% 4]
      first_unknown <- if (length(unknown)) min(unknown) else as.Date(NA)
      future <- d[!is.na(d$VisitDate) & d$VisitDate > index_date &
                    d$FCOGIMPR %in% c(0, 1, 3) &
                    (is.na(first_unknown) | d$VisitDate < first_unknown),
                  , drop = FALSE]
      if (!nrow(future)) next
      steps["later_clinical_visit"] <- steps["later_clinical_visit"] + 1L
      future <- future[order(future$VisitDate), , drop = FALSE]
      dementia <- which(future$FCOGIMPR == 3)
      event <- as.integer(length(dementia) > 0L)
      endpoint <- if (event) future$VisitDate[dementia[1L]] else
        max(future$VisitDate)
      followup <- as.numeric(endpoint - index_date) / 365.25
      k <- k + 1L
      records[[k]] <- data.frame(
        id = second$MCSA_ID, stage = if (stage == 0) "CN" else "MCI",
        MRI1_date = first$ImagingDate, MRI2_date = index_date,
        MRI_gap_years = chosen$days / 365.25,
        visit_MRI2_gap_days = as.numeric(index_date - second$VisitDate),
        age_MRI2 = second$calc_age_vis,
        male = second$Male, education = second$EDUC,
        baseline_CDRSB = second$CDRSUM,
        hippo2_pct_ICV = 100 * second$HPVOL / second$TIV,
        hippo_loss_pct_year = 100 * log(first$HPVOL / second$HPVOL) /
          (chosen$days / 365.25),
        same_manufacturer = as.integer(first$MRIMANUFACTURER ==
                                         second$MRIMANUFACTURER),
        followup_years = followup, event = event,
        stringsAsFactors = FALSE
      )
    }
    risk <- if (k) do.call(rbind, records[seq_len(k)]) else data.frame()
    list(risk = risk,
         flow = data.frame(stage = if (stage == 0) "CN" else "MCI",
                           gate = names(steps), participants = as.integer(steps)))
  }

  cn <- build_risk_set(0)
  risk <- cn$risk
  save_tsv(cn$flow, "Table1_landmark_sample_flow.tsv")
  save_tsv(risk, "Table2_person_level_frozen_risk_sets.tsv")
  summary_tab <- do.call(rbind, lapply(split(risk, risk$stage), function(d) {
    data.frame(stage = d$stage[1L], participants = nrow(d),
               events = sum(d$event),
               same_manufacturer = sum(d$same_manufacturer),
               median_MRI_gap_years = median(d$MRI_gap_years),
               median_followup_years = median(d$followup_years))
  }))
  save_tsv(summary_tab, "Table3_clinical_and_amyloid_availability.tsv")

  fit_stage <- function(d) {
    m0 <- survival::coxph(
      survival::Surv(followup_years, event) ~ age_MRI2 + male + education +
        baseline_CDRSB + hippo2_pct_ICV, data = d, ties = "efron")
    m1 <- survival::coxph(
      survival::Surv(followup_years, event) ~ age_MRI2 + male + education +
        baseline_CDRSB + hippo2_pct_ICV + hippo_loss_pct_year,
      data = d, ties = "efron")
    coef <- summary(m1)$coefficients["hippo_loss_pct_year", ]
    ci <- exp(stats::confint(m1)["hippo_loss_pct_year", ])
    lr <- 2 * (as.numeric(logLik(m1)) - as.numeric(logLik(m0)))
    data.frame(cohort = "MCSA", stage = d$stage[1L], participants = nrow(d),
               events = sum(d$event),
               HR_per_1pct_year = exp(coef["coef"]), CI_low = ci[1L],
               CI_high = ci[2L], Wald_p = coef["Pr(>|z|)"],
               LR_chisq = lr, LR_p = pchisq(lr, 1, lower.tail = FALSE),
               apparent_C_base = summary(m0)$concordance[1L],
               apparent_C_plus_loss = summary(m1)$concordance[1L],
               stringsAsFactors = FALSE)
  }
  message("Estimating the frozen OASIS-style association in MCSA")
  results <- do.call(rbind, lapply(split(risk, risk$stage), fit_stage))
  save_tsv(results, "Table4_MCSA_event_associations.tsv")

  oasis <- read.delim(oasis_path, stringsAsFactors = FALSE)
  comparator <- data.frame(cohort = "OASIS-3", stage = "CN",
                           participants = oasis$participants[1L],
                           events = oasis$events[1L],
                           HR_per_1pct_year = oasis$HR_per_1pct_year[1L],
                           CI_low = oasis$CI_low[1L],
                           CI_high = oasis$CI_high[1L],
                           Wald_p = oasis$Wald_p[1L])
  comparison <- rbind(comparator, results[names(comparator)])
  save_tsv(comparison, "Table5_OASIS_MCSA_event_comparison.tsv")

  writeLines(capture.output(sessionInfo()),
             file.path(output_dir, "R_sessionInfo.txt"))
  message("Analysis complete: ", output_dir)
  print(summary_tab, row.names = FALSE)
  print(results, row.names = FALSE)
})
