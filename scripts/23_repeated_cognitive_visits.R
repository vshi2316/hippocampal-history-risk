local({
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("lme4", quietly = TRUE))
    stop("Install data.table and lme4 before running this script.")
  library(data.table)

  root <- .study$output
  paths <- c(
    oasis_risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    oasis_cognition = file.path(.study$output,"44_OASIS3_full_clinical_link_results/09_B4_CDR_by_visit.csv"),
    mayo_risk = file.path(root, "mcsa_event_43/run_frozen/tables/Table2_person_level_frozen_risk_sets.tsv"),
    mayo_zip = file.path(.study$data,"mayo","Demographics_and_Assessments.zip")
  )
  if (any(!file.exists(paths)))
    stop("Missing input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "cognitive_trajectory_63",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  save_tsv(data.table(input = names(paths), file = unname(paths),
                      md5 = unname(tools::md5sum(paths))),
           "Input_manifest.tsv", "audit")

  oasis <- fread(paths[["oasis_risk"]])[baseline_stage == "CN"]
  mayo <- fread(paths[["mayo_risk"]])[stage == "CN"]
  if (nrow(oasis) != 312L || sum(as.logical(oasis$event)) != 52L ||
      nrow(mayo) != 925L || sum(mayo$event) != 21L ||
      anyDuplicated(oasis$OASISID) || anyDuplicated(mayo$id))
    stop("The frozen MRI2 risk sets changed.")

  oasis_risk <- oasis[, .(id = OASISID, index = MRI2_day,
                          event_end = followup_days,
                          age_MRI2, sex = factor(sex), education,
                          baseline_CDRSB, hippo2_pct_ICV,
                          hippo_loss_pct_year)]
  ov <- fread(paths[["oasis_cognition"]])
  if (!all(c("OASISID", "day", "mmse") %in% names(ov)))
    stop("OASIS-3 MMSE fields changed.")
  ov <- ov[, .(id = OASISID, visit = day,
               mmse = suppressWarnings(as.numeric(mmse)))]
  ov <- merge(ov, oasis_risk, by = "id")
  ov[, `:=`(time_days = visit - index, cohort = "OASIS-3")]
  ov[mmse < 0 | mmse > 30, mmse := NA_real_]

  mayo[, `:=`(index = as.Date(MRI2_date),
              event_end = followup_years * 365.25,
              sex = factor(male))]
  mayo_risk <- mayo[, .(id, index, event_end, age_MRI2, sex,
                        education, baseline_CDRSB, hippo2_pct_ICV,
                        hippo_loss_pct_year)]
  raw <- as.data.table(read.csv(unz(paths[["mayo_zip"]],
                                  "MCSA_Data_12Nov2026.csv"),
                                stringsAsFactors = FALSE, check.names = FALSE))
  needed <- c("MCSA_ID", "VisitDate", "FCOGIMPR", "pzglobal",
              "pzmemory", "pzlanguage", "pzattention", "pzvisualsp",
              "faq_total", "MMSEcalc")
  if (!all(needed %in% names(raw)))
    stop("Mayo cognitive assessment fields changed.")
  raw[, VisitDate := as.Date(VisitDate)]
  for (name in setdiff(needed, c("MCSA_ID", "VisitDate")))
    set(raw, j = name,
        value = suppressWarnings(as.numeric(raw[[name]])))
  raw[, memory_specificity := pzmemory - rowMeans(
    as.matrix(raw[, .(pzlanguage, pzattention, pzvisualsp)]))]
  unknown <- merge(raw[FCOGIMPR == 4 & !is.na(VisitDate),
                       .(id = MCSA_ID, VisitDate)],
                   mayo_risk[, .(id, index)], by = "id")
  unknown <- unknown[VisitDate > index,
                     .(first_unknown = min(VisitDate)), by = id]
  score_names <- c("pzglobal", "pzmemory", "memory_specificity",
                   "faq_total", "MMSEcalc")
  mv <- raw[!is.na(VisitDate), c("MCSA_ID", "VisitDate", score_names),
            with = FALSE]
  setnames(mv, c("MCSA_ID", "VisitDate"), c("id", "visit"))
  mv <- mv[, lapply(.SD, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) median(x) else NA_real_
  }), by = .(id, visit), .SDcols = score_names]
  mv <- merge(mv, mayo_risk, by = "id")
  mv <- merge(mv, unknown, by = "id", all.x = TRUE)
  mv <- mv[is.na(first_unknown) | visit < first_unknown]
  mv[, `:=`(time_days = as.numeric(visit - index), cohort = "Mayo")]
  mv[faq_total < 0 | faq_total > 30, faq_total := NA_real_]
  mv[MMSEcalc < 0 | MMSEcalc > 30, MMSEcalc := NA_real_]
  mv[, log1p_faq := log1p(faq_total)]

  make_long <- function(d, measure, baseline_window_days, role,
                        until_event = FALSE) {
    x <- copy(d[is.finite(get(measure))])
    baseline <- x[time_days <= 0 & time_days >= -baseline_window_days]
    setorder(baseline, id, time_days)
    baseline <- baseline[, .SD[.N], by = id,
                         .SDcols = c("time_days", measure)]
    baseline <- baseline[, .(id, baseline_score = get(measure),
                             baseline_gap_days = -time_days)]
    future <- x[time_days > 0 & time_days <= 1827]
    if (until_event) future <- future[time_days <= event_end + 1e-6]
    future <- merge(future, baseline, by = "id")
    if (!nrow(future)) stop("No matched future scores: ", role)
    setorder(future, id, time_days)
    future <- unique(future, by = c("id", "time_days"))
    future[, .(cohort, role, id, time_years = time_days / 365.25,
               score = get(measure), baseline_score, baseline_gap_days,
               age_MRI2, sex, education, baseline_CDRSB,
               hippo2_pct_ICV, hippo_loss_pct_year)]
  }
  specs <- list(
    list(d = mv, field = "pzglobal", window = 90L,
         role = "Primary_Mayo_global_z"),
    list(d = mv, field = "pzmemory", window = 90L,
         role = "Mayo_memory_z"),


    list(d = ov, field = "mmse", window = 365L,
         role = "OASIS_administered_MMSE")
  )
  long <- rbindlist(lapply(specs, function(s)
    make_long(s$d, s$field, s$window, s$role)), fill = TRUE)
  save_tsv(long, "Table2_person_visit_cognitive_trajectories.tsv")
  coverage <- long[, .(future_visits = .N,
                       last_year = max(time_years)),
                   by = .(cohort, role, id)][,
                     .(participants = .N,
                       with_2plus_future_visits = sum(future_visits >= 2L),
                       with_3plus_future_visits = sum(future_visits >= 3L),
                       future_observations = sum(future_visits),
                       median_last_year = median(last_year)),
                     by = .(cohort, role)]
  save_tsv(coverage, "Table1_cognitive_followup_coverage.tsv")

writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
})
