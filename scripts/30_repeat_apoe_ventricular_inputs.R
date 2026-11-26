# Prepare clinical covariates and matched scans for the control analyses.
local({
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("Install required R package: data.table")
  library(data.table)

  root <- .study$output
  oasis_root <- .study$output
  cache <- file.path(root, "oasis_tau_mri_gate_27/run_frozen/audit/input_csv_cache")
  paths <- c(
    clinical = file.path(root, "cn_progression_62/run_frozen/tables/Table1_person_level_clinical_course.tsv"),
    mayo_risk = file.path(root, "mcsa_event_43/run_frozen/tables/Table2_person_level_frozen_risk_sets.tsv"),
    mayo_zip = file.path(.study$data,"mayo","Demographics_and_Assessments.zip"),
    oasis_risk = file.path(oasis_root, "45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    oasis_cdr = file.path(oasis_root, "44_OASIS3_full_clinical_link_results/09_B4_CDR_by_visit.csv"),
    oof = file.path(root, "history_increment_66/run_frozen/tables/Table2_person_level_out_of_fold_risk.tsv"),
    increment = file.path(root, "history_increment_66/run_frozen/tables/Table1_paired_out_of_fold_increment.tsv"),
    oasis_demo_zip = file.path(.study$data,"oasis3/download.zip"),
    oasis_muse_zip = file.path(.study$data,"oasis3/ADOPIC_Analysis.zip"),
    oasis_muse = file.path(cache, "MUSE_volume_OASIS_04Nov2026.csv"),
    oasis_mri_meta = file.path(cache, "OASIS3_MR_scan_summary.csv"),
    oasis_crosswalk = file.path(root, "oasis_event_40/run_frozen/source_data/Table6_all_312_exact_scan_crosswalk.tsv")
  )
  required <- c("clinical", "mayo_risk", "mayo_zip", "oasis_risk", "oasis_cdr", "oof", "increment")
  absent <- required[!file.exists(paths[required])]
  if (length(absent)) stop("Missing frozen input: ", paste(paths[absent], collapse = "; "))
  out <- file.path(root, "clinical_field_audit_69",
                   paste0("run_", "frozen"))
  if (dir.exists(out)) stop("Output directory already exists: ", out)
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA", quote = TRUE)
  number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  one_number <- function(x) {
    y <- unique(number(x)[is.finite(number(x))])
    if (length(y) == 1L) y else NA_real_
  }
  read_zip <- function(path, member, optional = FALSE) {
    if (!file.exists(path)) {
      if (optional) return(NULL)
      stop("Missing archive: ", path)
    }
    members <- unzip(path, list = TRUE)$Name
    exact <- members[basename(members) == member]
    if (length(exact) != 1L) {
      if (optional && !length(exact)) return(NULL)
      stop("Archive member is missing or ambiguous: ", member)
    }
    scratch <- tempfile("step69_zip_")
    dir.create(scratch)
    extracted <- file.path(scratch, basename(exact))
    on.exit({ unlink(extracted); unlink(scratch, recursive = FALSE) })
    unzip(path, files = exact, exdir = scratch, junkpaths = TRUE)
    fread(extracted, check.names = FALSE, na.strings = c("", "NA", "NaN", "NULL"))
  }
  manifest <- data.table(input = names(paths), file = unname(paths), exists = file.exists(paths))
  manifest[, md5 := NA_character_]
  manifest[exists == TRUE, md5 := unname(tools::md5sum(file))]
  save_tsv(manifest, "Input_manifest.tsv", "audit")
  clinical <- fread(paths[["clinical"]])
  risk <- copy(clinical[prior_MCI_record == 0L])
  expected <- data.table(cohort = c("Mayo", "OASIS-3"), participants = c(925L, 306L), events = c(177L, 59L))
  actual <- risk[, .(participants = .N, events = sum(impairment_event)), by = cohort]
  check <- merge(expected, actual, by = "cohort", suffixes = c("_expected", "_actual"))
  if (nrow(check) != 2L || any(check$participants_expected != check$participants_actual) ||
      any(check$events_expected != check$events_actual) || anyDuplicated(risk, by = c("cohort", "id")))
    stop("Frozen first-impairment risk sets changed.")
  mayo_anchor <- fread(paths[["mayo_risk"]])[stage == "CN",
    .(id, MRI1_date = as.Date(MRI1_date), MRI2_date = as.Date(MRI2_date), MRI_gap_years,
      same_manufacturer = as.integer(same_manufacturer))]
  oasis_anchor <- fread(paths[["oasis_risk"]])[baseline_stage == "CN",
    .(id = OASISID, MRI1_day, MRI2_day, MRI1_Code, MRI2_Code)]
  mayo_anchor <- merge(risk[cohort == "Mayo", .(id)], mayo_anchor, by = "id")
  oasis_anchor <- merge(risk[cohort == "OASIS-3", .(id)], oasis_anchor, by = "id")
  if (nrow(mayo_anchor) != 925L || nrow(oasis_anchor) != 306L ||
      anyDuplicated(mayo_anchor$id) || anyDuplicated(oasis_anchor$id)) stop("MRI anchors changed.")

  mayo <- read_zip(paths[["mayo_zip"]], "MCSA_Data_12Nov2026.csv")
  mayo_dictionary <- read_zip(paths[["mayo_zip"]], "MCSA_Data_Dictionary_12Nov2026.csv")
  mayo[, `:=`(MCSA_ID = as.character(MCSA_ID), VisitDate = as.Date(VisitDate), ImagingDate = as.Date(ImagingDate))]
  save_tsv(mayo_dictionary, "Mayo_data_dictionary.tsv", "audit")
  mayo_long <- merge(mayo, mayo_anchor, by.x = "MCSA_ID", by.y = "id")
  setnames(mayo_long, "MCSA_ID", "id")
  mb <- mayo_long[is.finite(number(pzglobal)) & !is.na(VisitDate),
    .(baseline_score = median(number(pzglobal)), MRI2_date = MRI2_date[1L]), by = .(id, VisitDate)]
  mb[, gap_days := as.numeric(MRI2_date - VisitDate)]
  mb <- mb[gap_days >= 0 & gap_days <= 90]
  setorder(mb, id, gap_days, VisitDate)
  mb <- unique(mb, by = "id")[, .(cohort = "Mayo", id, baseline_score, baseline_gap_days = gap_days)]
  cdr <- fread(paths[["oasis_cdr"]])
  ob <- merge(cdr[is.finite(mmse) & mmse >= 0 & mmse <= 30,
                  .(id = OASISID, day, baseline_score = mmse)], oasis_anchor[, .(id, MRI2_day)], by = "id")
  ob[, gap_days := MRI2_day - day]
  ob <- ob[gap_days >= 0 & gap_days <= 365]
  setorder(ob, id, gap_days, day)
  ob <- unique(ob, by = "id")[, .(cohort = "OASIS-3", id, baseline_score, baseline_gap_days = gap_days)]
  baseline <- rbindlist(list(mb, ob))
  risk <- merge(risk, baseline, by = c("cohort", "id"), all.x = TRUE)
  risk[, cognition_complete := is.finite(baseline_score)]
  if (risk[cohort == "Mayo" & cognition_complete == TRUE, .N] != 880L ||
      risk[cohort == "OASIS-3" & cognition_complete == TRUE, .N] != 306L ||
      risk[cohort == "Mayo" & cognition_complete == TRUE, sum(impairment_event)] != 168L ||
      risk[cohort == "OASIS-3" & cognition_complete == TRUE, sum(impairment_event)] != 59L)
    stop("Frozen cognition-complete risk sets changed.")
  predictions <- fread(paths[["oof"]])
  oof_mri2 <- predictions[landmark == "MRI2"]
  eligible_keys <- risk[cognition_complete == TRUE, paste(cohort, id)]
  if (anyDuplicated(oof_mri2, by = c("cohort", "id")) ||
      !setequal(eligible_keys, oof_mri2[, paste(cohort, id)])) stop("Step 66 identities changed.")
  identity <- merge(risk[cognition_complete == TRUE], oof_mri2, by = c("cohort", "id"))
  if (any(abs(identity$impairment_time_years - identity$time) > 1e-8) ||
      any(identity$impairment_event != identity$event)) stop("Step 66 event times changed.")
  risk[, population := "Full_first_impairment"]
  populations <- rbindlist(list(risk, copy(risk[cognition_complete == TRUE])[
    , population := "Cognition_complete_MRI2"]))
  coverage_rows <- list()
  add_coverage <- function(cohort_name, module, item, ids, status, note) {
    coverage_rows[[length(coverage_rows) + 1L]] <<- populations[cohort == cohort_name,
      .(participants_total = .N, participants_available = sum(id %in% ids),
        events_total = sum(impairment_event), events_available = sum(impairment_event[id %in% ids]),
        events_by_3y_available = sum(impairment_event[id %in% ids] == 1L &
                                      impairment_time_years[id %in% ids] <= 3)), by = .(cohort, population)][
      , `:=`(module = module, item = item, status = status, note = note)]
    invisible(NULL)
  }
  save_tsv(populations[, .(participants = .N, events = sum(impairment_event),
                            events_by_3y = sum(impairment_event == 1L & impairment_time_years <= 3)),
                        by = .(cohort, population)], "Table2_frozen_risk_set_check.tsv")
  baseline_summary <- risk[cognition_complete == TRUE,
    .(participants = .N, median = median(baseline_score), Q1 = as.numeric(quantile(baseline_score, .25)),
      Q3 = as.numeric(quantile(baseline_score, .75)), minimum = min(baseline_score), maximum = max(baseline_score),
      at_MMSE_30 = if (cohort[1L] == "OASIS-3") sum(baseline_score == 30) else NA_integer_,
      median_days_before_MRI2 = median(baseline_gap_days)), by = cohort]
  baseline_summary[, score := ifelse(cohort == "Mayo", "Global cognitive z score", "MMSE")]
  save_tsv(baseline_summary, "Table3_baseline_cognition_for_Table1.tsv")
  save_tsv(risk, "Frozen_MRI2_risk_with_baseline_score.tsv", "source_data")
  save_tsv(fread(paths[["increment"]]), "Table12_existing_Step66_performance_and_net_benefit.tsv")

  # APOE uses documented carrier coding in Mayo and explicit allele pairs in OASIS-3.
  apoe_mayo <- mayo_long[!is.na(VisitDate) & VisitDate <= MRI2_date & number(Any_E4) %in% c(0, 1),
    .(apoe_e4 = one_number(Any_E4), distinct_values = uniqueN(number(Any_E4))), by = id]
  apoe_mayo[, cohort := "Mayo"]
  add_coverage("Mayo", "APOE", "Any_E4_before_MRI2", apoe_mayo[distinct_values == 1L, id],
               "READY_FOR_SAME_SAMPLE_COMPARISON", "0/1 coding verified in the Mayo dictionary; conflicting records excluded.")
  demo <- read_zip(paths[["oasis_demo_zip"]], "OASIS3_demographics.csv", optional = TRUE)
  apoe_oasis <- data.table(id = character(), apoe_e4 = numeric(), distinct_values = integer(), cohort = character())
  if (!is.null(demo)) {
    demo <- merge(demo, oasis_anchor[, .(id)], by.x = "OASISID", by.y = "id")
    setnames(demo, "OASISID", "id")
    if (anyDuplicated(demo$id)) stop("OASIS demographic ID is not unique.")
    demo[, allele_pair := gsub("[^0-9]", "", as.character(APOE))]
    demo[!grepl("^[234]{2}$", allele_pair), allele_pair := NA_character_]
    demo[, candidate_e4 := fifelse(is.na(allele_pair), NA_real_, as.numeric(grepl("4", allele_pair)))]
    apoe_oasis <- demo[is.finite(candidate_e4), .(apoe_e4 = one_number(candidate_e4),
      distinct_values = uniqueN(candidate_e4)), by = id][, cohort := "OASIS-3"]
    save_tsv(demo[, .N, by = .(APOE, allele_pair)], "OASIS_APOE_recorded_codes.tsv", "audit")
    add_coverage("OASIS-3", "APOE", "APOE_allele_pair", apoe_oasis[distinct_values == 1L, id],
                 "READY_FOR_SAME_SAMPLE_COMPARISON", "Static genotype; explicit two-allele codes only. No acquisition date inferred.")
  } else add_coverage("OASIS-3", "APOE", "APOE_allele_pair", character(), "SOURCE_MISSING", "Demographic archive/member not found.")
  save_tsv(rbindlist(list(apoe_mayo, apoe_oasis), fill = TRUE), "Table4_person_level_APOE_availability.tsv", "source_data")

  # Use exact MUSE scan codes to align the anatomical comparator.
  muse <- if (file.exists(paths[["oasis_muse"]])) fread(paths[["oasis_muse"]]) else
    read_zip(paths[["oasis_muse_zip"]], "MUSE_volume_OASIS_04Nov2026.csv", optional = TRUE)
  if (!is.null(muse)) {
    dict <- read_zip(paths[["oasis_muse_zip"]], "MUSE_volume_OASIS_Dictionary_04Nov2026.csv", optional = TRUE)
    if (!is.null(dict)) save_tsv(dict, "OASIS_MUSE_dictionary.tsv", "audit")
    if (anyDuplicated(muse$Code)) stop("MUSE scan Code is not unique.")
    muse[, `:=`(id = paste0("OAS", sub(".*OAS-[0-9]-([0-9]{5})-.*", "\\1", as.character(Code))),
                 mri_day = as.integer(as.Date(EXAMDATE) - as.Date("2000-01-01")))]
    requested <- rbindlist(list(oasis_anchor[, .(id, scan = "MRI1", Code = MRI1_Code, anchor_day = MRI1_day)],
                                oasis_anchor[, .(id, scan = "MRI2", Code = MRI2_Code, anchor_day = MRI2_day)]))
    selected <- merge(requested, muse, by = c("id", "Code"), all.x = TRUE)
    if (any(selected$mri_day != selected$anchor_day, na.rm = TRUE)) stop("MUSE Code/day differs from the frozen anchor.")
    selected[, lateral_ventricles := number(Left_Lateral_Ventricle) + number(Right_Lateral_Ventricle)]
    selected[, hippocampus := number(Left_Hippocampus) + number(Right_Hippocampus)]
    selected[, matched := is.finite(mri_day)]
    hippo_check <- dcast(selected, id ~ scan, value.var = c("hippocampus", "ICV"))
    hippo_check <- merge(hippo_check, oasis_anchor[, .(id, interval_years = (MRI2_day - MRI1_day) / 365.25)], by = "id")
    hippo_check <- merge(hippo_check, risk[cohort == "OASIS-3", .(id, hippo2_pct_ICV, hippo_loss_pct_year)], by = "id")
    hippo_check[, `:=`(reconstructed_state = 100 * hippocampus_MRI2 / ICV_MRI2,
                       reconstructed_loss = 100 * log(hippocampus_MRI1 / hippocampus_MRI2) / interval_years)]
    if (hippo_check[is.finite(reconstructed_state), any(abs(reconstructed_state - hippo2_pct_ICV) > 1e-6)] ||
        hippo_check[is.finite(reconstructed_loss), any(abs(reconstructed_loss - hippo_loss_pct_year) > 1e-6)])
      stop("Exact MUSE pairs do not reconstruct the frozen hippocampal exposure.")
    save_tsv(hippo_check, "OASIS_MUSE_reconstructs_frozen_exposure.tsv", "audit")
    anatomical <- c("lateral_ventricles")
    pair_tables <- lapply(anatomical, function(field) {
      d <- selected[, .(id, scan, value = number(get(field)), ICV = number(ICV), matched)]
      wide <- dcast(d, id ~ scan, value.var = c("value", "ICV", "matched"))
      wide <- merge(wide, oasis_anchor[, .(id, interval_years = (MRI2_day - MRI1_day) / 365.25)], by = "id")
      wide[, usable_pair := matched_MRI1 %in% TRUE & matched_MRI2 %in% TRUE &
             is.finite(value_MRI1) & value_MRI1 > 0 & is.finite(value_MRI2) & value_MRI2 > 0 &
             is.finite(ICV_MRI1) & ICV_MRI1 > 0 & is.finite(ICV_MRI2) & ICV_MRI2 > 0 & interval_years > 0]
      wide[, `:=`(cohort = "OASIS-3", field = field,
                   annual_log_change_pct = NA_real_, annual_log_loss_pct = NA_real_)]
      wide[usable_pair == TRUE, annual_log_change_pct := 100 * log(value_MRI2 / value_MRI1) / interval_years]
      wide[, annual_log_loss_pct := -annual_log_change_pct]
      add_coverage("OASIS-3", "Anatomical_comparator", field, wide[usable_pair == TRUE, id],
                   "EXACT_SCAN_PAIR_AVAILABLE", "Values in cubic mm. Positive change means expansion; positive loss means shrinkage. Comparator not yet selected or modeled.")
      wide
    })
    save_tsv(rbindlist(pair_tables, fill = TRUE), "Table6_exact_scan_anatomical_pairs.tsv", "source_data")
    save_tsv(selected[, .(id, Code, scan, anchor_day, mri_day, matched, hippocampus, ICV)],
             "OASIS_selected_scan_identity.tsv", "audit")
    rep <- muse[id %in% oasis_anchor$id & is.finite(mri_day),
      .(scans = .N, distinct_codes = uniqueN(Code),
        codes = paste(sort(as.character(Code)), collapse = " | "),
        positive_hippocampal_scans = sum(is.finite(number(Left_Hippocampus) + number(Right_Hippocampus)) &
          number(Left_Hippocampus) + number(Right_Hippocampus) > 0)), by = .(id, mri_day)][distinct_codes >= 2L]
    anchor_days <- unique(requested[, .(id, mri_day = anchor_day)])
    rep[, is_frozen_anchor_day := paste(id, mri_day) %in% anchor_days[, paste(id, mri_day)]]
    save_tsv(rep, "Table7_OASIS_same_day_repeat_acquisitions.tsv", "source_data")
    add_coverage("OASIS-3", "Measurement", "Same_day_repeat_at_MRI1_or_MRI2", rep[is_frozen_anchor_day == TRUE, id],
                 "REPEAT_ACQUISITIONS_REQUIRE_METADATA_MATCH", "Repeat acquisitions are not repeat segmentation of one image. Match device/sequence before precision analysis.")
  } else {
    add_coverage("OASIS-3", "Anatomical_comparator", "MUSE", character(), "SOURCE_MISSING", "Neither pinned cache nor source archive was found.")
  }

  mmse <- merge(cdr[is.finite(mmse) & mmse >= 0 & mmse <= 30, .(id = OASISID, day, mmse)],
                oasis_anchor[, .(id, MRI2_day)], by = "id")
  mmse <- merge(mmse, risk[cohort == "OASIS-3", .(id, impairment_event, impairment_time_years)], by = "id")
  mmse[, time_years := (day - MRI2_day) / 365.25]
  mmse <- unique(mmse[time_years > 0 & time_years <= 5 &
    ((impairment_event == 1L & time_years < impairment_time_years - 1 / 365.25) |
       (impairment_event == 0L & time_years <= impairment_time_years + 1 / 365.25))], by = c("id", "day", "mmse"))
  mmse_summary <- mmse[, .(participants = uniqueN(id), observations = .N,
    observations_at_30 = sum(mmse == 30), proportion_at_30 = mean(mmse == 30),
    median = median(mmse), Q1 = as.numeric(quantile(mmse, .25)), Q3 = as.numeric(quantile(mmse, .75)))]
  save_tsv(mmse_summary, "Table9_OASIS_MMSE_ceiling_before_first_event.tsv")

  intervals <- rbindlist(list(mayo_anchor[, .(cohort = "Mayo", id, interval_years = MRI_gap_years)],
                             oasis_anchor[, .(cohort = "OASIS-3", id, interval_years = (MRI2_day - MRI1_day) / 365.25)]))
  intervals <- merge(intervals, risk[, .(cohort, id, impairment_event, impairment_time_years)], by = c("cohort", "id"))
  if (any(!is.finite(intervals$interval_years) | intervals$interval_years <= 0)) stop("Invalid frozen MRI interval.")
  intervals[, interval_bin := cut(interval_years, breaks = c(0, 1, 2, 3, 5, Inf), right = FALSE,
                                  labels = c("Below_1y", "1_to_below_2y", "2_to_below_3y", "3_to_below_5y", "5y_or_more"))]
  save_tsv(intervals[, .(participants = .N, events = sum(impairment_event),
    events_by_3y = sum(impairment_event == 1L & impairment_time_years <= 3),
    interval_median = median(interval_years)), by = .(cohort, interval_bin)], "Table11_MRI_interval_coverage_no_effect_test.tsv")
  save_tsv(intervals, "MRI_interval_in_final_risk_set.tsv", "source_data")
  add_coverage("Mayo", "Measurement", "MRI_manufacturer", mayo_anchor[!is.na(same_manufacturer), id],
               "MANUFACTURER_METADATA_ONLY", "Same manufacturer is not same device or repeated segmentation reliability.")
  for (cohort_name in c("Mayo", "OASIS-3"))
    add_coverage(cohort_name, "Measurement", "Scan_specific_segmentation_QC_and_processing_version", character(),
                 "NOT_FOUND_IN_SCOPED_VOLUME_TABLE", "Dataset update stamps and scanner software fields do not identify segmentation pass/fail or processing versions.")
  if (file.exists(paths[["oasis_crosswalk"]])) {
    scanner <- merge(fread(paths[["oasis_crosswalk"]]), oasis_anchor[, .(OASISID = id)], by = "OASISID")
    save_tsv(scanner, "OASIS_frozen_scan_metadata.tsv", "source_data")
    add_coverage("OASIS-3", "Measurement", "MRI_device_and_sequence", scanner[both_exact %in% TRUE, OASISID],
                 "METADATA_AVAILABLE", "Segmentation QC and pipeline version require their own evidence; scanner metadata do not establish them.")
  }
  readiness <- rbindlist(coverage_rows, fill = TRUE)
  setcolorder(readiness, c("cohort", "population", "module", "item", "status",
                           "participants_total", "participants_available", "events_total",
                           "events_available", "events_by_3y_available", "note"))
  save_tsv(readiness, "Table1_analysis_readiness_and_overlap.tsv")
  capture.output(sessionInfo(), file = file.path(out, "audit", "R_session_info.txt"))
  writeLines(c("STEP_69_COMPLETED", normalizePath(out, winslash = "/", mustWork = TRUE)),
             file.path(out, "69_COMPLETED.txt"), useBytes = TRUE)
  cat("Step 69 complete: ", out, "\n", sep = "")
  print(readiness[, .(cohort, population, module, item, status, participants_available, events_available)])
  cat("\nMain summary: tables/Table1_analysis_readiness_and_overlap.tsv\n")
})
