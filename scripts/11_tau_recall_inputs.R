local({
# -*- coding: UTF-8 -*-
# Prepare tau PET with dated ADNI MMSE and HABS-HD delayed recall z-scores.

for (p in c("data.table", "readxl", "sandwich", "lmtest")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install R package: ", p)
}
library(data.table)

root <- .study$data
out <- file.path(file.path(.study$output,"tau_recall_inputs"),
                 paste0("run_", "frozen"))
if (dir.exists(out)) stop("Output directory already exists: ", out)
dir.create(out, recursive = TRUE)

write_tsv <- function(x, name) {
  fwrite(as.data.table(x), file.path(out, name), sep = "\t", na = "NA")
}
num <- function(x) {
  y <- suppressWarnings(as.numeric(x))
  y[y %in% c(-9999, -777777, -7777)] <- NA_real_
  y
}
midmonth <- function(year, month) {
  m <- suppressWarnings(as.integer(month))
  bad <- is.na(m)
  if (any(bad)) m[bad] <- match(tolower(substr(as.character(month[bad]), 1, 3)),
                               tolower(month.abb))
  as.Date(sprintf("%04d-%02d-15", as.integer(year), m))
}
mean_complete <- function(x, columns) {
  z <- as.matrix(x[, ..columns])
  storage.mode(z) <- "numeric"
  ans <- rowMeans(z, na.rm = TRUE)
  ans[rowSums(is.finite(z)) < length(columns) * 0.75] <- NA_real_
  ans[!is.finite(ans) | ans < 1 | ans > 5] <- NA_real_
  ans
}
nearest <- function(anchor, other, other_date, window_days) {
  x <- merge(anchor, other, by = "id", allow.cartesian = TRUE)
  x[, date_gap := as.integer(get(other_date) - tau_date)]
  x <- x[!is.na(date_gap) & abs(date_gap) <= window_days]
  x[, abs_gap := abs(date_gap)]
  setorder(x, id, abs_gap, date_gap)
  x <- x[, .SD[1L], by = id]
  x[, abs_gap := NULL]
  x
}
flow <- data.table(cohort = character(), step = character(), n = integer())
record <- function(cohort, step, x) {
  flow <<- rbind(flow, data.table(cohort, step, n = uniqueN(x$id)))
}

regions <- list(
  mtl = c("entorhinal", "parahippocampal", "inferiortemporal", "middletemporal"),
  pcc = c("posteriorcingulate", "isthmuscingulate"),
  lpar = c("inferiorparietal", "superiorparietal", "supramarginal")
)

read_adni_mri <- function() {
  base <- file.path(root, "adni/MRI")
  dict <- fread(file.path(base, "MRI_FreeSurfer_Data_Dictionary.csv"))
  dict <- dict[grepl("^Thickness Average of (Left|Right)", TEXT)]
  dict[, side := ifelse(grepl("^Thickness Average of Left", TEXT), "L", "R")]
  dict[, roi := tolower(sub("^Thickness Average of (Left|Right)", "", TEXT))]
  fields <- unique(unlist(lapply(regions, function(rois) {
    dict[roi %in% rois, FLDNAME]
  })))
  if (length(fields) != 18L) stop("ADNI dictionary does not map all 18 regional thickness columns")
  x <- fread(file.path(base, "MRI_FreeSurfer.csv"),
             select = c("RID", "EXAMDATE", "IMAGEUID", "OVERALLQC", "STATUS",
                        "FIELD_STRENGTH", "FSVER", fields))
  x <- x[!tolower(OVERALLQC) %in% c("fail", "hippocampus only")]
  x[, `:=`(id = as.character(RID), mri_date = as.Date(EXAMDATE),
           image_id = as.character(IMAGEUID),
           scanner = paste(FIELD_STRENGTH, FSVER, sep = "_"),
           mri_qc = fifelse(is.na(OVERALLQC) | OVERALLQC == "",
                            "Not provided", as.character(OVERALLQC)))]
  for (region in names(regions)) {
    columns <- dict[roi %in% regions[[region]], FLDNAME]
    x[, (region) := mean_complete(.SD, columns), .SDcols = columns]
  }
  x <- x[!is.na(mri_date) & complete.cases(x[, .(mtl, pcc, lpar)])]
  x <- x[, .(id, mri_date, image_id, mri_qc, scanner, mtl, pcc, lpar)]
  unique(x, by = c("id", "mri_date", "image_id"))
}

read_habs_mri <- function() {
  path <- file.path(root, "habs_hd/HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx")
  x <- rbindlist(lapply(1:4, function(i) {
    as.data.table(readxl::read_xlsx(path, sheet = paste("Inclusions - Scan", i)))
  }), fill = TRUE)
  x[, `:=`(id = as.character(Med_ID), mri_date = midmonth(MRI_year, MRI_month),
           image_id = as.character(Time_point_MRI), age_mri = num(Age),
           scanner = as.character(Scanner_MRI), mri_qc = "Included")]
  for (region in names(regions)) {
    columns <- unlist(lapply(regions[[region]], function(roi) {
      paste0(c("L_", "R_"), roi, "_thickavg")
    }))
    if (!all(columns %in% names(x))) stop("HABS thickness columns missing: ", region)
    x[, (region) := mean_complete(.SD, columns), .SDcols = columns]
  }
  x <- x[!is.na(mri_date) & complete.cases(x[, .(mtl, pcc, lpar)])]
  x <- x[, .(id, mri_date, image_id, mri_qc, age_mri, scanner, mtl, pcc, lpar)]
  unique(x, by = c("id", "mri_date", "image_id"))
}

read_adni_tau <- function() {
  x <- fread(file.path(root, "adni/PET/PET_Tau_Berkeley.csv"),
             select = c("RID", "SCANDATE", "TRACER", "META_TEMPORAL_SUVR",
                        "CTX_POSTERIORCINGULATE_SUVR", "CTX_INFERIORPARIETAL_SUVR",
                        "CTX_SUPERIORPARIETAL_SUVR", "CTX_SUPRAMARGINAL_SUVR"))
  x <- x[TRACER == "FTP"]
  x[, `:=`(id = as.character(RID), tau_date = as.Date(SCANDATE),
           tau_mtl = num(META_TEMPORAL_SUVR),
           tau_pcc = num(CTX_POSTERIORCINGULATE_SUVR))]
  x[, tau_lpar := rowMeans(cbind(num(CTX_INFERIORPARIETAL_SUVR),
                                 num(CTX_SUPERIORPARIETAL_SUVR),
                                 num(CTX_SUPRAMARGINAL_SUVR)), na.rm = FALSE)]
  x <- x[!is.na(tau_date) & complete.cases(x[, .(tau_mtl, tau_pcc, tau_lpar)]) &
           tau_mtl > 0 & tau_pcc > 0 & tau_lpar > 0]
  unique(x[, .(id, tau_date, tau_mtl, tau_pcc, tau_lpar)], by = c("id", "tau_date"))
}

read_habs_tau <- function() {
  path <- file.path(root, "habs_hd/HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx")
  x <- rbindlist(lapply(1:2, function(i) {
    as.data.table(readxl::read_xlsx(path, sheet = paste0("Scan", i, " - Inclusion")))
  }), fill = TRUE)
  x[, `:=`(id = as.character(Med_ID),
           tau_date = midmonth(PI2620_PET_year, PI2620_PET_month),
           tau_mtl = num(Medial_Temporal_PI2620_SUVR),
           tau_pcc = num(Posterior_Cingulate_PI2620_SUVR),
           tau_lpar = num(Lateral_Parietal_PI2620_SUVR))]
  x <- x[!is.na(tau_date) & complete.cases(x[, .(tau_mtl, tau_pcc, tau_lpar)]) &
           tau_mtl > 0 & tau_pcc > 0 & tau_lpar > 0]
  unique(x[, .(id, tau_date, tau_mtl, tau_pcc, tau_lpar)], by = c("id", "tau_date"))
}

read_adni_amyloid <- function() {
  x <- fread(file.path(root, "adni/PET/PET_Amyloid_Berkeley.csv"),
             select = c("RID", "SCANDATE", "AMYLOID_STATUS", "CENTILOIDS", "qc_flag"))
  x[, `:=`(id = as.character(RID), ab_date = as.Date(SCANDATE),
           amyloid_pos = num(AMYLOID_STATUS), amyloid_cl = num(CENTILOIDS))]
  unique(x[!is.na(ab_date) & amyloid_pos %in% 0:1 & qc_flag >= 1,
           .(id, ab_date, amyloid_pos, amyloid_cl)], by = c("id", "ab_date"))
}

read_habs_amyloid <- function() {
  path <- file.path(root, "habs_hd/HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx")
  x <- rbindlist(lapply(1:2, function(i) {
    as.data.table(readxl::read_xlsx(path, sheet = paste("Scan", i, "- Inclusion")))
  }), fill = TRUE)
  x[, `:=`(id = as.character(Med_ID),
           ab_date = midmonth(FBB_PET_year, FBB_PET_month),
           amyloid_pos = num(AB_pos_FBB),
           amyloid_cl = num(Centiloid_Units_GlobalAmyloid))]
  unique(x[!is.na(ab_date) & amyloid_pos %in% 0:1,
           .(id, ab_date, amyloid_pos, amyloid_cl)], by = c("id", "ab_date"))
}

read_adni_clinical <- function() {
  x <- fread(file.path(root, "adni/Clinical/ADNIMERGE_03Jan2026.csv"),
             select = c("RID", "EXAMDATE", "DX", "AGE", "PTGENDER", "MMSE"))
  x[, `:=`(id = as.character(RID), visit_date = as.Date(EXAMDATE),
           diagnosis = as.character(DX), age = num(AGE),
           sex = as.character(PTGENDER), memory = num(MMSE))]
  x[, age := age + as.numeric(visit_date - min(visit_date, na.rm = TRUE)) / 365.25,
    by = id]
  unique(x[!is.na(visit_date),
           .(id, visit_date, diagnosis, age, sex, memory)],
         by = c("id", "visit_date"))
}

read_habs_clinical <- function() {
  base <- file.path(root, "habs_hd")
  dx <- fread(file.path(base, "dx_12Aug2026.csv"),
              select = c("MED_ID", "VISIT_ID", "VISIT_MONTH", "VISIT_YEAR", "CDX_COG"))
  cog <- fread(file.path(base, "zscores_12Aug2026.csv"),
               select = c("MED_ID", "VISIT_ID", "SEVLT_DR_ZSCORE"))
  demo <- fread(file.path(base, "interviewdemographics_12Aug2026.csv"),
                select = c("MED_ID", "VISIT_ID", "ID_GENDER"))
  x <- merge(dx, cog, by = c("MED_ID", "VISIT_ID"), all.x = TRUE)
  x <- merge(x, demo, by = c("MED_ID", "VISIT_ID"), all.x = TRUE)
  x[, `:=`(id = as.character(MED_ID),
           visit_date = midmonth(VISIT_YEAR, VISIT_MONTH),
           diagnosis = fifelse(num(CDX_COG) == 0, "CN",
                               fifelse(num(CDX_COG) == 1, "MCI", "Other")),
           sex = fifelse(num(ID_GENDER) %in% 0:1,
                         as.character(ID_GENDER), NA_character_),
           memory = num(SEVLT_DR_ZSCORE))]
  unique(x[!is.na(visit_date),
           .(id, visit_date, diagnosis, sex, memory)],
         by = c("id", "visit_date"))
}

make_cohort <- function(cohort, mri, tau, amyloid, clinical) {
  write_tsv(mri[, .N, by = mri_qc], paste0(cohort, "_MRI_QC_counts.tsv"))
  record(cohort, "valid_tau_PET", tau)
  pairs <- merge(tau, mri, by = "id", allow.cartesian = TRUE)
  pairs[, mri_tau_gap := as.integer(mri_date - tau_date)]
  pairs <- pairs[!is.na(mri_tau_gap) & abs(mri_tau_gap) <= 365L]
  pairs[, abs_gap := abs(mri_tau_gap)]
  setorder(pairs, id, tau_date, abs_gap)
  pairs <- pairs[, .SD[1L], by = .(id, tau_date)]
  setorder(pairs, id, tau_date)
  anchor <- pairs[, .SD[1L], by = id]
  anchor[, abs_gap := NULL]
  record(cohort, "tau_with_baseline_MRI_365d", anchor)
  if (!nrow(anchor)) stop(cohort, ": no tau/MRI pairs")

  ab <- nearest(anchor, amyloid, "ab_date", 730L)
  setnames(ab, "date_gap", "amyloid_tau_gap")
  record(cohort, "plus_amyloid_PET_730d", ab)
  if (!nrow(ab)) stop(cohort, ": no amyloid-linked pairs")

  dx <- nearest(ab, clinical[, setdiff(names(clinical), "memory"), with = FALSE],
                "visit_date", 365L)
  setnames(dx, "date_gap", "diagnosis_tau_gap")
  dx <- dx[diagnosis %in% c("CN", "MCI")]
  record(cohort, "nondemented_baseline_diagnosis", dx)
  if (!nrow(dx)) stop(cohort, ": no CN/MCI baseline pairs")
  setnames(dx, "visit_date", "diagnosis_date")
  if (!"age" %in% names(dx)) dx[, age := age_mri]
  dx[, tau_mean := rowMeans(cbind(tau_mtl, tau_pcc, tau_lpar))]
  dx[, tau_contrast := tau_mtl - (tau_pcc + tau_lpar) / 2]
  write_tsv(dx, paste0(cohort, "_01_baseline_tau_MRI_amyloid.tsv"))

  fu <- copy(mri)
  setnames(fu, c("mri_date", "image_id", "mri_qc", "scanner", "mtl", "pcc", "lpar"),
           c("follow_date", "follow_image_id", "follow_mri_qc", "follow_scanner",
             "mtl_follow", "pcc_follow", "lpar_follow"))
  fu <- fu[, .(id, follow_date, follow_image_id, follow_mri_qc, follow_scanner,
               mtl_follow, pcc_follow, lpar_follow)]
  longitudinal <- merge(dx, fu, by = "id", allow.cartesian = TRUE)
  longitudinal <- longitudinal[
    follow_date >= mri_date + 365L & follow_date >= tau_date + 365L &
      follow_date <= tau_date + 1826L]
  setorder(longitudinal, id, follow_date)
  longitudinal <- longitudinal[, .SD[1L], by = id]
  record(cohort, "plus_followup_MRI_1to5y_after_tau", longitudinal)
  record(cohort, "MRI_followup_amyloid_positive", longitudinal[amyloid_pos == 1])
  longitudinal[, `:=`(
    mri_interval_years = as.numeric(follow_date - mri_date) / 365.25,
    tau_to_follow_days = as.integer(follow_date - tau_date),
    scanner_changed = as.integer(scanner != follow_scanner))]
  for (region in names(regions)) {
    longitudinal[, (paste0("thin_", region)) :=
                   (get(region) - get(paste0(region, "_follow"))) / mri_interval_years]
  }
  longitudinal[, thin_contrast := thin_mtl - (thin_pcc + thin_lpar) / 2]
  write_tsv(longitudinal, paste0(cohort, "_02_tau_to_MRI_change.tsv"))

  cognitive <- clinical[!is.na(memory), .(id, memory_date = visit_date, memory)]
  baseline_cog <- nearest(dx, cognitive, "memory_date", 365L)
  setnames(baseline_cog, c("memory_date", "memory"),
           c("baseline_memory_date", "baseline_memory"))
  follow_cog <- copy(cognitive)
  setnames(follow_cog, c("memory_date", "memory"),
           c("follow_memory_date", "follow_memory"))
  trajectory <- merge(baseline_cog, follow_cog, by = "id", allow.cartesian = TRUE)
  trajectory <- trajectory[follow_memory_date >= tau_date + 365L &
                             follow_memory_date <= tau_date + 1826L &
                             follow_memory_date > baseline_memory_date]
  trajectory[, target_gap := abs(as.integer(follow_memory_date - (tau_date + 1096L)))]
  setorder(trajectory, id, target_gap)
  trajectory <- trajectory[, .SD[1L], by = id]
  trajectory[, memory_change_year :=
               (follow_memory - baseline_memory) /
               (as.numeric(follow_memory_date - baseline_memory_date) / 365.25)]
  record(cohort, "plus_future_memory_1to5y", trajectory)
  record(cohort, "memory_followup_amyloid_positive", trajectory[amyloid_pos == 1])
  write_tsv(trajectory, paste0(cohort, "_03_tau_to_future_memory.tsv"))
  list(longitudinal = longitudinal, cognition = trajectory)
}
adni <- make_cohort("ADNI", read_adni_mri(), read_adni_tau(),
                    read_adni_amyloid(), read_adni_clinical())
habs <- make_cohort("HABS_HD", read_habs_mri(), read_habs_tau(),
                    read_habs_amyloid(), read_habs_clinical())
write_tsv(flow, "00_sample_flow.tsv")

})
