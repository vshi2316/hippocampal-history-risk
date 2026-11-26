required <- c("data.table")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages))
  stop("Install R packages: ", paste(missing_packages, collapse = ", "))

main <- function() {
  code_dir <- .study$output
  data_dir <- file.path(.study$data,"oasis3")
  out <- file.path(code_dir, "44_OASIS3_full_clinical_link_results")
  archives <- c(muse = file.path(data_dir, "ADOPIC_Analysis.zip"),
                clinical = file.path(data_dir, "Neuropsychological.zip"),
                demographic = file.path(data_dir, "download.zip"),
                imaging = file.path(data_dir, "Imaging.zip"))
  if (any(!file.exists(archives)))
    stop("Missing archive: ", paste(archives[!file.exists(archives)], collapse = "; "))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  writeLines("STATUS=STEP44_INCOMPLETE", file.path(out, "00_status.txt"))
  save_csv <- function(x, name)
    data.table::fwrite(x, file.path(out, name), na = "")
  num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  int <- function(x) suppressWarnings(as.integer(as.character(x)))
  read_member <- function(zip, member) {
    destination <- tempfile("oasis_member_")
    dir.create(destination)
    on.exit(unlink(destination, recursive = TRUE), add = TRUE)
    target <- utils::unzip(zip, files = member, exdir = destination)
    if (length(target) != 1L || !file.exists(target) || file.info(target)$size == 0)
      stop("Cannot extract archive member: ", member)
    data.table::as.data.table(utils::read.csv(
      target, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8"))
  }
  check_columns <- function(x, wanted, source_name) {
    absent <- setdiff(wanted, names(x))
    if (length(absent)) stop(source_name, " lacks: ", paste(absent, collapse = ", "))
  }
  info <- file.info(archives)
  save_csv(data.frame(archive = names(archives), path = unname(archives),
                      bytes = info$size, modified = format(info$mtime,
                                                           "%Y-%m-%d %H:%M:%S")),
           "01_input_manifest.csv")

  muse <- read_member(archives[["muse"]], "MUSE_volume_OASIS_04Nov2026.csv")
  demo <- read_member(archives[["demographic"]], "OASIS3_demographics.csv")
  d1 <- read_member(archives[["clinical"]], "OASIS3_UDSd1_diagnoses.csv")
  b4 <- read_member(archives[["clinical"]], "OASIS3_UDSb4_cdr.csv")
  mr_raw <- read_member(archives[["imaging"]], "OASIS3_MR_scan_summary.csv")
  check_columns(muse, c("SUBJECTID", "Code", "EXAMDATE", "ICV",
                        "Left_Hippocampus", "Right_Hippocampus"), "MUSE")
  check_columns(demo, c("OASISID", "AgeatEntry", "GENDER", "EDUC"), "demographics")
  check_columns(d1, c("OASISID", "days_to_visit", "DEMENTED", "NORMCOG"), "D1")
  check_columns(b4, c("OASISID", "days_to_visit", "CDRSUM", "MMSE"), "B4")
  check_columns(mr_raw, c("subject_id", "label", "filename"), "MR summary")
  if (anyDuplicated(demo$OASISID)) stop("Duplicate demographic subject IDs")
  demo[, `:=`(AgeatEntry = num(AgeatEntry), EDUC = num(EDUC),
              GENDER = as.character(GENDER))]

  mr_pattern <- "^OAS[0-9]{5}_MR_d([0-9]+)$"
  t1 <- unique(mr_raw[grepl(mr_pattern, label) &
                        grepl("_T1w\\.json$", filename), .(
    OASISID = as.character(subject_id),
    mri_day = int(sub(mr_pattern, "\\1", label)),
    run = int(ifelse(grepl("_run-[0-9]+_T1w\\.json$", filename),
                     sub("^.*_run-([0-9]+)_T1w\\.json$", "\\1", filename),
                     "0")))])
  data.table::setorder(t1, OASISID, mri_day, run)
  save_csv(t1, "03_T1_metadata_inventory.csv")

  muse[, OASISID := sub("^OAS-0-([0-9]{5})-.*$", "OAS\\1", Code)]
  muse[, code_day := int(sub(
    "^OAS-0-[0-9]{5}-m[0-9]+-([0-9]{4})[0-9]{3}$", "\\1", Code))]
  run_pattern <- "^OAS-0-[0-9]{5}-m[0-9]+-[0-9]{4}([0-9]{2})0$"
  muse[, run := int(ifelse(grepl(run_pattern, Code),
                           sub(run_pattern, "\\1", Code), NA_character_))]
  muse[, date_day := as.integer(as.Date(EXAMDATE) - as.Date("2000-01-01"))]
  muse[, subject_match := int(SUBJECTID) ==
         int(sub("^OAS", "", OASISID)) - 30000L]
  muse[, day_match := !is.na(code_day) & !is.na(date_day) &
         code_day == date_day]
  muse[, mri_day := code_day]
  muse[, has_T1_meta := paste(OASISID, mri_day, run, sep = "|") %in%
         paste(t1$OASISID, t1$mri_day, t1$run, sep = "|")]
  muse[, `:=`(ICV = num(ICV), Left_Hippocampus = num(Left_Hippocampus),
              Right_Hippocampus = num(Right_Hippocampus))]
  muse[, hippocampus := Left_Hippocampus + Right_Hippocampus]
  muse[, positive := is.finite(ICV) & ICV > 0 &
         is.finite(Left_Hippocampus) & Left_Hippocampus > 0 &
         is.finite(Right_Hippocampus) & Right_Hippocampus > 0]
  muse[, eligible := !is.na(subject_match) & subject_match & day_match &
         has_T1_meta == TRUE & positive == TRUE]
  save_csv(muse[, .(OASISID, SUBJECTID, Code, EXAMDATE, mri_day, run,
                    subject_match, day_match, has_T1_meta, positive, eligible,
                    ICV, Left_Hippocampus, Right_Hippocampus)],
           "04_MUSE_QC_and_provenance.csv")
  scan <- muse[eligible == TRUE, .(OASISID, mri_day, run, Code,
                                    ICV, hippocampus)]
  data.table::setorder(scan, OASISID, mri_day, run, Code)
  duplicate_scans <- scan[duplicated(scan, by = c("OASISID", "mri_day")) |
                            duplicated(scan, by = c("OASISID", "mri_day"),
                                       fromLast = TRUE)]
  save_csv(duplicate_scans, "05_same_day_MUSE_repeats.csv")
  scan <- unique(scan, by = c("OASISID", "mri_day"))
  save_csv(scan, "06_selected_MUSE_scans.csv")

  mci_columns <- grep("^MCI", names(d1), value = TRUE)
  if (!length(mci_columns)) stop("D1 has no MCI flags")
  d1[, day := int(days_to_visit)]
  d1 <- d1[!is.na(day) & !is.na(OASISID) & nzchar(OASISID)]
  d1[, dementia_flag := as.character(DEMENTED) == "1"]
  d1[, normal_flag := as.character(NORMCOG) == "1"]
  d1[, mci_flag := rowSums(as.matrix(d1[, mci_columns, with = FALSE]) == "1",
                           na.rm = TRUE) > 0]
  d1[, other_impairment_flag := if ("IMPNOMCI" %in% names(d1))
         as.character(IMPNOMCI) == "1" else FALSE]
  dx <- d1[, .(dementia_flag = any(dementia_flag %in% TRUE),
               mci_flag = any(mci_flag %in% TRUE),
               normal_flag = any(normal_flag %in% TRUE),
               other_impairment_flag = any(other_impairment_flag %in% TRUE),
               source_rows = .N), by = .(OASISID, day)]
  dx[, conflict := as.integer(dementia_flag) + as.integer(mci_flag) +
         as.integer(normal_flag) > 1L]
  dx[, stage := ifelse(conflict, "Conflict",
                ifelse(dementia_flag, "Dementia",
                ifelse(mci_flag, "MCI",
                ifelse(normal_flag, "CN",
                ifelse(other_impairment_flag, "Other impairment", "Unknown")))))]
  data.table::setorder(dx, OASISID, day)
  save_csv(dx, "07_D1_all_visit_stages.csv")
  save_csv(dx[conflict == TRUE], "07b_D1_diagnostic_conflicts.csv")
  writeLines(paste(mci_columns, collapse = ","),
             file.path(out, "07c_MCI_flag_columns.txt"))

  b4[, `:=`(day = int(days_to_visit),
             cdrsb_raw = num(CDRSUM), mmse_raw = num(MMSE))]
  b4[, valid_cdr := !is.na(day) & is.finite(cdrsb_raw) &
         cdrsb_raw >= 0 & cdrsb_raw <= 18]
  b4[, valid_mmse := !is.na(day) & is.finite(mmse_raw) &
         mmse_raw >= 0 & mmse_raw <= 30]
  save_csv(b4[, .(OASISID, days_to_visit, day, cdrsb_raw, mmse_raw,
                  valid_cdr, valid_mmse)], "08_B4_clinical_measure_inventory.csv")
  cdr <- b4[valid_cdr == TRUE, .(
    cdrsb = as.numeric(stats::median(cdrsb_raw)),
    mmse = as.numeric(stats::median(mmse_raw[valid_mmse == TRUE], na.rm = TRUE)),
    source_rows = .N,
    cdr_within_day_range = as.numeric(diff(range(cdrsb_raw)))),
    by = .(OASISID, day)]
  data.table::setorder(cdr, OASISID, day)
  save_csv(cdr, "09_B4_CDR_by_visit.csv")
  build_landmarks <- function(two_mri = FALSE) {
    result <- vector("list", data.table::uniqueN(scan$OASISID))
    k <- 0L
    for (id in unique(scan$OASISID)) {
      s <- scan[OASISID == id]
      person <- demo[OASISID == id]
      if (!nrow(person)) next
      all_dx <- dx[OASISID == id]
      all_cdr <- cdr[OASISID == id]
      for (j in seq_len(nrow(s))) {
        second <- s[j]
        prior <- s[mri_day <= second$mri_day - 180L]
        if (two_mri && !nrow(prior)) next
        first <- if (two_mri) prior[nrow(prior)] else second
        pre_dx <- all_dx[day <= second$mri_day &
                           day >= second$mri_day - 365L]
        pre_cdr <- all_cdr[day <= second$mri_day &
                             day >= second$mri_day - 365L]
        if (!nrow(pre_dx) || !nrow(pre_cdr)) next
        baseline_dx <- pre_dx[nrow(pre_dx)]
        baseline_cdr <- pre_cdr[nrow(pre_cdr)]
        if (!baseline_dx$stage %in% c("CN", "MCI") ||
            any(all_dx$day <= second$mri_day & all_dx$dementia_flag)) next
        later_dx <- all_dx[day > second$mri_day &
                             stage %in% c("CN", "MCI", "Dementia")]
        later_event <- later_dx[stage == "Dementia"]
        later_change <- later_dx[stage %in% c("MCI", "Dementia")]
        loss <- if (two_mri) 100 * log(first$hippocampus /
                                         second$hippocampus) /
                 ((second$mri_day - first$mri_day) / 365.25) else NA_real_
        k <- k + 1L
        result[[k]] <- data.frame(
          OASISID = id, MRI1_day = if (two_mri) first$mri_day else NA_integer_,
          MRI2_day = second$mri_day,
          MRI1_Code = if (two_mri) first$Code else NA_character_,
          MRI2_Code = second$Code,
          MRI_span_days = if (two_mri) second$mri_day - first$mri_day else NA_integer_,
          hippo1_mm3 = if (two_mri) first$hippocampus else NA_real_,
          hippo2_mm3 = second$hippocampus,
          hippo2_pct_ICV = 100 * second$hippocampus / second$ICV,
          hippo_loss_pct_year = loss,
          age_MRI2 = person$AgeatEntry + second$mri_day / 365.25,
          sex = person$GENDER, education = person$EDUC,
          baseline_D1_day = baseline_dx$day,
          baseline_stage = baseline_dx$stage,
          baseline_CDR_day = baseline_cdr$day,
          baseline_CDRSB = baseline_cdr$cdrsb,
          post_D1_visits = nrow(later_dx),
          last_D1_day = if (nrow(later_dx)) max(later_dx$day) else NA_integer_,
          first_dementia_day = if (nrow(later_event))
            later_event$day[1L] else NA_integer_,
          first_CN_progression_day = if (baseline_dx$stage == "CN" &&
                                          nrow(later_change))
            later_change$day[1L] else NA_integer_)
        break
      }
    }
    if (!k) return(data.table::data.table())
    data.table::as.data.table(do.call(rbind, result[seq_len(k)]))
  }

  two <- build_landmarks(TRUE)
  if (!nrow(two)) stop("No eligible clinical landmarks; inspect source inventories")
  add_outcomes <- function(x) {
    x[, `:=`(event = !is.na(first_dementia_day),
              followup_days = ifelse(!is.na(first_dementia_day),
                                     first_dementia_day - MRI2_day,
                                     last_D1_day - MRI2_day),
              D1_gap_days = MRI2_day - baseline_D1_day,
              CDR_gap_days = MRI2_day - baseline_CDR_day,
              extreme_MRI_change = abs(hippo_loss_pct_year) > 10)]
    x[, two_year_known := !is.na(followup_days) & followup_days >= 730L |
         event & followup_days <= 730L]
    x[, two_year_event := event & followup_days <= 730L]
    x[, five_year_known := !is.na(followup_days) & followup_days >= 1826L |
         event & followup_days <= 1826L]
    x[, five_year_event := event & followup_days <= 1826L]
    x
  }
  two <- add_outcomes(two)
  save_csv(two, "12_two_MRI_clinical_landmarks.csv")
  risk <- two[!is.na(followup_days) & followup_days > 0 &
                is.finite(hippo_loss_pct_year) & is.finite(age_MRI2) &
                is.finite(baseline_CDRSB)]
  risk[, baseline_stage := factor(baseline_stage)]
  event_count <- sum(risk$event)
  save_csv(risk[, .(OASISID, baseline_stage, MRI2_day, hippo_loss_pct_year,
                    age_MRI2, baseline_CDRSB, followup_days, event)],
           "19b_event_risk_set.csv")

writeLines("STATUS=STEP44_COMPLETE",file.path(out,"00_status.txt"))
}
main()
