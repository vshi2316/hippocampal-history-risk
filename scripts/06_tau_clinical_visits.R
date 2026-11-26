local({
  required <- c("data.table")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages(library(data.table))

  root <- .study$output
  input <- c(
    scan16 = file.path(root, "tau_hippocampal_models/run_frozen/analysis_data/Frozen_scan_audit.tsv"),
    long11 = file.path(root, "tau_longitudinal_inputs/run_frozen/analysis_data/MRI_longitudinal_analysis.tsv"),
    effect16 = file.path(root, "tau_hippocampal_models/run_frozen/tables/Table2_tau_slopes_standardized_and_raw.tsv"),
    influence16 = file.path(root, "tau_influence/run_frozen/tables/Table2_Stage16_tau_by_time_influence.tsv"),
    adni = file.path(.study$data,"adni/Clinical/ADNIMERGE_03Jan2026.csv"),
    habs_mmse = file.path(.study$data,"habs_hd/mmse_12Aug2026.csv"),
    habs_cdr = file.path(.study$data,"habs_hd/consensuscdr_12Aug2026.csv"),
    habs_memory = file.path(.study$data,"habs_hd/totals_12Aug2026.csv")
  )
  if (any(!file.exists(input))) stop("Missing input: ", paste(names(input)[!file.exists(input)], collapse = ", "))
  out <- file.path(root, "clinical_intersection_26",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "analysis_data", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t",
           quote = TRUE, na = "NA", bom = TRUE)
  number <- function(x) suppressWarnings(as.numeric(as.character(x)))
  clean_id <- function(x) sub("\\.0+$", "", trimws(as.character(x)))
  month_date <- function(year, month) {
    y <- suppressWarnings(as.integer(year))
    m <- suppressWarnings(as.integer(month))
    z <- trimws(as.character(month))
    use <- which(is.na(m) & !is.na(z))
    if (length(use)) m[use] <- match(tolower(substr(z[use], 1L, 3L)),
                                   tolower(month.abb))
    ans <- as.Date(rep(NA_real_, length(y)), origin = "1970-01-01")
    use <- which(is.finite(y) & y >= 1900 & y <= 2100 &
                 is.finite(m) & m %in% 1:12)
    ans[use] <- as.Date(sprintf("%04d-%02d-15", y[use], m[use]))
    ans
  }
  scans <- fread(input[["scan16"]], encoding = "UTF-8")[
    scope %in% c("ADNI_primary", "HABS_primary"),
    .(cohort, id, mri_date, image_id, stable_pipeline)]
  scans[, `:=`(cohort = as.character(cohort), id = as.character(id),
               image_id = as.character(image_id), mri_date = as.Date(mri_date))]
  long <- fread(input[["long11"]], encoding = "UTF-8")[,
    .(cohort, id, mri_date, image_id, baseline_date, tau_landmark,
      log_hip_ratio, amyloid_pos, reference_partition)]
  long[, `:=`(cohort = as.character(cohort), id = as.character(id),
              image_id = as.character(image_id), mri_date = as.Date(mri_date),
              baseline_date = as.Date(baseline_date), tau_landmark = as.Date(tau_landmark))]
  matched <- merge(scans, long,
                   by = c("cohort", "id", "mri_date", "image_id"), all.x = TRUE)
  if (nrow(matched) != nrow(scans) || anyNA(matched$baseline_date) ||
      any(matched$amyloid_pos != 1L) ||
      any(matched$reference_partition == "training"))
    stop("The frozen Stage 16 scan population did not match Stage 11")
  setorder(matched, cohort, id, mri_date, image_id)
  mri1 <- matched[mri_date == baseline_date, .SD[1L], by = .(cohort, id)]
  mri2 <- matched[mri_date >= tau_landmark + 365L &
                    mri_date <= tau_landmark + 1826L &
                    mri_date >= baseline_date + 180L & is.finite(log_hip_ratio),
                  .SD[1L], by = .(cohort, id)]
  mri1 <- mri1[, .(cohort, id, mri1_date = mri_date,
                    mri1_image_id = image_id, tau_date = tau_landmark)]
  mri2 <- mri2[, .(cohort, id, mri2_date = mri_date,
                    mri2_image_id = image_id)]
  anchors <- merge(mri1, mri2, by = c("cohort", "id"), all.x = TRUE)
  anchors[, index_date := pmax(mri1_date, tau_date)]
  frozen_n <- anchors[, .(n_people = as.integer(.N)), by = cohort]
  save_tsv(anchors, "Frozen_Stage16_participant_MRI_anchors.tsv", "analysis_data")

  adni <- fread(input[["adni"]],
                select = c("RID", "EXAMDATE", "MMSE", "CDRSB", "LDELTOTAL"),
                encoding = "UTF-8")
  habs_mmse <- fread(input[["habs_mmse"]],
                     select = c("MED_ID", "VISIT_YEAR", "VISIT_MONTH", "MMSE_TOTAL"),
                     encoding = "UTF-8")
  habs_cdr <- fread(input[["habs_cdr"]],
                    select = c("MED_ID", "VISIT_YEAR", "VISIT_MONTH", "CDR_SUM"),
                    encoding = "UTF-8")
  habs_memory <- fread(input[["habs_memory"]],
                       select = c("MED_ID", "VISIT_YEAR", "VISIT_MONTH", "LM2_A_TOTAL"),
                       encoding = "UTF-8")
  aid <- paste0("ADNI_", clean_id(adni$RID))
  adate <- as.Date(adni$EXAMDATE)
  visits <- rbindlist(list(
    data.table(cohort = "ADNI", id = aid, domain = "MMSE", visit_date = adate,
               score = number(adni$MMSE)),
    data.table(cohort = "ADNI", id = aid, domain = "CDR-SB", visit_date = adate,
               score = number(adni$CDRSB)),
    data.table(cohort = "ADNI", id = aid, domain = "Delayed recall", visit_date = adate,
               score = number(adni$LDELTOTAL)),
    data.table(cohort = "HABS-HD", id = paste0("HABS_", clean_id(habs_mmse$MED_ID)),
               domain = "MMSE", visit_date = month_date(habs_mmse$VISIT_YEAR,
                                                        habs_mmse$VISIT_MONTH),
               score = number(habs_mmse$MMSE_TOTAL)),
    data.table(cohort = "HABS-HD", id = paste0("HABS_", clean_id(habs_cdr$MED_ID)),
               domain = "CDR-SB", visit_date = month_date(habs_cdr$VISIT_YEAR,
                                                          habs_cdr$VISIT_MONTH),
               score = number(habs_cdr$CDR_SUM)),
    data.table(cohort = "HABS-HD", id = paste0("HABS_", clean_id(habs_memory$MED_ID)),
               domain = "Delayed recall", visit_date = month_date(habs_memory$VISIT_YEAR,
                                                                   habs_memory$VISIT_MONTH),
               score = number(habs_memory$LM2_A_TOTAL))
  ))
  visits <- visits[!is.na(visit_date) & is.finite(score) &
    ((domain == "MMSE" & score >= 0 & score <= 30) |
     (domain == "CDR-SB" & score >= 0 & score <= 18) |
     (domain == "Delayed recall" & score >= 0 & score <= 25))]
  conflicts <- visits[, .(n_rows = as.integer(.N),
                          distinct_scores = as.integer(uniqueN(score))),
                      by = .(cohort, id, domain, visit_date)][distinct_scores > 1L]
  save_tsv(conflicts, "Conflicting_same_day_clinical_scores.tsv", "audit")
  if (nrow(conflicts)) visits <- visits[!conflicts,
                                      on = .(cohort, id, domain, visit_date)]
  setorder(visits, cohort, id, domain, visit_date)
  visits <- unique(visits, by = c("cohort", "id", "domain", "visit_date"))
  paired <- merge(anchors, visits, by = c("cohort", "id"), allow.cartesian = TRUE)
  paired[, `:=`(days_from_index = as.integer(visit_date - index_date),
                days_from_mri2 = as.integer(visit_date - mri2_date),
                target_index_730 = index_date + 730L,
                target_mri2_730 = mri2_date + 730L)]
  save_tsv(paired[, .(cohort, id, domain, mri1_date, tau_date, index_date,
                      mri2_date, visit_date, days_from_index, days_from_mri2,
                      score)], "All_valid_clinical_visits_on_frozen_timeline.tsv",
           "analysis_data")

  choose <- function(x, date_col, score_col, target_col) {
    if (!nrow(x)) return(NULL)
    x <- copy(x)
    x[, distance := abs(as.integer(visit_date - get(target_col)))]
    setorder(x, cohort, id, domain, distance, visit_date)
    x <- unique(x, by = c("cohort", "id", "domain"))[
      , .(cohort, id, domain, selected_date = visit_date, selected_score = score)]
    setnames(x, c("selected_date", "selected_score"), c(date_col, score_col))
    x
  }
  base <- choose(paired[abs(days_from_index) <= 90L],
                 "baseline_date_clinical", "baseline_score", "index_date")
  future <- choose(paired[days_from_index >= 365L & days_from_index <= 1826L],
                   "future_index_date", "future_index_score", "target_index_730")
  pre_mri2 <- choose(paired[days_from_mri2 >= -180L & days_from_mri2 <= 0L],
                     "pre_mri2_date", "pre_mri2_score", "mri2_date")
  post_mri2 <- choose(paired[days_from_mri2 >= 365L & days_from_mri2 <= 1096L],
                      "post_mri2_date", "post_mri2_score", "target_mri2_730")
  domains <- c("MMSE", "CDR-SB", "Delayed recall")
  person <- anchors[rep(seq_len(nrow(anchors)), each = length(domains))]
  person[, domain := rep(domains, times = nrow(anchors))]
  for (part in list(base, future, pre_mri2, post_mri2))
    if (!is.null(part)) person <- merge(person, part,
                                       by = c("cohort", "id", "domain"), all.x = TRUE)
  for (col in c("baseline_date_clinical", "future_index_date", "pre_mri2_date",
                "post_mri2_date"))
    if (!col %in% names(person)) person[, (col) := as.Date(NA)]
  for (col in c("baseline_score", "future_index_score", "pre_mri2_score",
                "post_mri2_score"))
    if (!col %in% names(person)) person[, (col) := NA_real_]
  person[, `:=`(has_mri2 = !is.na(mri2_date),
                baseline_index = !is.na(baseline_date_clinical),
                future_index = !is.na(future_index_date),
                baseline_mri2 = !is.na(pre_mri2_date),
                future_mri2 = !is.na(post_mri2_date))]
  person[, `:=`(index_change = future_index_score - baseline_score,
                mri2_change = post_mri2_score - pre_mri2_score)]
  save_tsv(person, "Table2_person_level_clinical_intersection.tsv", "analysis_data")

})
