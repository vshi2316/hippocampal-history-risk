# Cortical tau and age-associated trajectories in ADNI and HABS-HD.

local({
  packages <- c("data.table", "readxl", "lme4", "lmerTest")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages(library(data.table))

  root <- .study$output
  files <- c(
    frozen = file.path(root, "tau_hippocampal_models/run_frozen/analysis_data/Frozen_scan_audit.tsv"),
    longitudinal = file.path(root, "tau_longitudinal_inputs/run_frozen/analysis_data/MRI_longitudinal_analysis.tsv"),
    hippocampus = file.path(root, "tau_hippocampal_models/run_frozen/tables/Table2_tau_slopes_standardized_and_raw.tsv"),
    adni_mri = file.path(.study$data,"adni/MRI/MRI_FreeSurfer.csv"),
    adni_dict = file.path(.study$data,"adni/MRI/MRI_FreeSurfer_Data_Dictionary.csv"),
    habs_mri = file.path(.study$data,"habs_hd/HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx")
  )
  if (any(!file.exists(files))) stop("Missing input: ", paste(files[!file.exists(files)], collapse = ", "))
  out <- file.path(root, "roi68_spatial_34", paste0("run_", "frozen"))
  for (subdir in c("tables", "source_data", "audit"))
    dir.create(file.path(out, subdir), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, subdir = "tables")
    fwrite(as.data.table(x), file.path(out, subdir, name), sep = "\t", quote = TRUE, na = "NA", bom = TRUE)
  log_line <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
  numeric_value <- function(x) suppressWarnings(as.numeric(as.character(x)))
  month_date <- function(year, month) {
    y <- as.integer(numeric_value(year))
    m <- as.integer(numeric_value(month))
    missing <- which(is.na(m) & !is.na(month))
    if (length(missing)) m[missing] <- match(tolower(substr(as.character(month[missing]), 1L, 3L)),
                                            tolower(month.abb))
    as.Date(ifelse(is.finite(y) & is.finite(m) & m >= 1L & m <= 12L,
                   sprintf("%04d-%02d-15", y, m), NA_character_))
  }
  save_tsv(data.table(input = names(files), path = unname(files),
                      md5 = unname(tools::md5sum(files))), "Input_manifest.tsv", "audit")
  log_line("Mapping the 68 Desikan-Killiany cortical thickness fields")

  dict <- fread(files[["adni_dict"]], encoding = "UTF-8")
  map <- dict[grepl("^Thickness Average of (Left|Right)", TEXT) &
              grepl("^ST[0-9]+TA$", FLDNAME), .(adni_field = FLDNAME, label = TEXT)]
  map[, `:=`(side = ifelse(grepl("^Thickness Average of Left", label), "L", "R"),
             roi = tolower(sub("^Thickness Average of (Left|Right)", "", label)))]
  map[, habs_field := paste0(side, "_", roi, "_thickavg")]
  medial <- c("entorhinal", "parahippocampal")
  other_temporal <- c("bankssts", "fusiform", "inferiortemporal", "middletemporal",
                      "superiortemporal", "temporalpole", "transversetemporal")
  map[, band := fifelse(roi %in% medial, "Medial temporal",
                        fifelse(roi %in% other_temporal, "Other temporal", "Non-temporal"))]
  setorder(map, roi, side)
  if (nrow(map) != 68L || anyDuplicated(map[, .(side, roi)]))
    stop("ADNI cortical atlas does not contain 68 unique bilateral regions")
  save_tsv(map, "Table1_68_region_field_and_band_dictionary.tsv")

  long11 <- fread(files[["longitudinal"]], encoding = "UTF-8")
  frozen <- fread(files[["frozen"]], encoding = "UTF-8")
  for (x in list(long11, frozen)) {
    set(x, j = "id", value = as.character(x$id))
    set(x, j = "image_id", value = as.character(x$image_id))
    set(x, j = "mri_date", value = as.Date(x$mri_date))
  }
  keys <- c("cohort", "id", "mri_date", "image_id")
  if (anyDuplicated(long11, by = keys)) stop("Stage 11 has repeated MRI keys")
  frozen <- frozen[scope %in% c("ADNI_primary", "HABS_primary")]
  if (nrow(frozen) != 404L || uniqueN(frozen[cohort == "ADNI", id]) != 94L ||
      uniqueN(frozen[cohort == "HABS-HD", id]) != 57L || anyDuplicated(frozen, by = keys))
    stop("Stage 16 primary scan set differs from the frozen 94/277 and 57/127 samples")
  frozen_key <- frozen[, .(cohort, id, mri_date, image_id,
                           frozen_time = time_years, frozen_tau = tau_z,
                           frozen_cortex = cortex_mm)]
  pos <- merge(frozen_key, long11, by = keys, all.x = TRUE, sort = FALSE)
  if (nrow(pos) != 404L || anyNA(pos$time_years) ||
      any(abs(pos$time_years - pos$frozen_time) > 1e-7) ||
      any(abs(pos$tau_z - pos$frozen_tau) > 1e-7) ||
      any(abs(pos$cortex_mm - pos$frozen_cortex) > 1e-5) ||
      any(pos$amyloid_pos != 1L)) stop("Stage 16 scans do not reproduce Stage 11 source data")
  pos[, group := "Apositive"]
  neg <- copy(long11[amyloid_pos == 0L & reference_partition != "training"])
  neg[, group := "Anegative"]
  scans <- rbindlist(list(pos, neg), use.names = TRUE, fill = TRUE)
  scans[, age10 := age_c]
  if (anyDuplicated(scans, by = c("group", keys))) stop("Repeated analysis MRI keys")
  save_tsv(scans[, .(participants = uniqueN(id), scans = .N,
                       baseline_scans = sum(time_years == 0),
                       future_scans = sum(time_years > 0)), by = .(cohort, group)],
           "Table2_frozen_and_age_reference_scan_counts.tsv")
  log_line("Reading ADNI FreeSurfer and HABS-HD inclusion-sheet source rows")

  adni_fields <- map$adni_field
  a <- fread(files[["adni_mri"]], select = c("RID", "EXAMDATE", "IMAGEUID",
             "FIELD_STRENGTH", "FSVER", adni_fields), encoding = "UTF-8")
  a[, `:=`(id = paste0("ADNI_", sub("\\.0+$", "", as.character(RID))),
            mri_date = as.Date(EXAMDATE), image_id = as.character(IMAGEUID),
            scanner = as.character(FIELD_STRENGTH),
            processing_version = as.character(FSVER),
            source_row = .I)]
  a <- a[id %in% scans[cohort == "ADNI", id]]
  a <- a[, c("id", "mri_date", "image_id", "scanner", "processing_version",
             "source_row", adni_fields), with = FALSE]
  adni_join <- c("id", "mri_date", "image_id", "scanner", "processing_version")
  if (anyDuplicated(a, by = adni_join)) stop("Repeated ADNI FreeSurfer source key")

  sheet_names <- paste("Inclusions - Scan", 1:4)
  h <- rbindlist(lapply(sheet_names, function(sheet) {
    x <- as.data.table(readxl::read_excel(files[["habs_mri"]], sheet = sheet,
          col_types = "text", na = c("", "NA", "N/A"), .name_repair = "unique_quiet"))
    x <- x[!is.na(Med_ID)]
    set(x, j = "source_sheet", value = rep(sheet, nrow(x)))
    set(x, j = "source_row", value = seq_len(nrow(x)) + 1L)
    x
  }), use.names = TRUE, fill = TRUE)
  habs_fields <- map$habs_field
  if (!all(habs_fields %in% names(h))) stop("HABS-HD inclusion sheets lack mapped ROI fields")
  h[, `:=`(id = paste0("HABS_", sub("\\.0+$", "", as.character(Med_ID))),
            mri_date = month_date(MRI_year, MRI_month),
            image_id = as.character(Time_point_MRI), scanner = as.character(Scanner_MRI))]
  h <- h[id %in% scans[cohort == "HABS-HD", id]]
  h <- h[, c("id", "mri_date", "image_id", "scanner", "source_sheet",
             "source_row", habs_fields), with = FALSE]
  habs_join <- c("id", "mri_date", "image_id", "scanner")
  if (anyDuplicated(h, by = habs_join)) stop("Repeated HABS-HD cortical source key")
  for (field in adni_fields) set(a, j = field, value = numeric_value(a[[field]]))
  for (field in habs_fields) set(h, j = field, value = numeric_value(h[[field]]))

  make_roi <- function(target_cohort, source, join, fields, source_map) {
    target <- copy(scans[scans$cohort == target_cohort])
    wide <- merge(target, source, by = join, all.x = TRUE, sort = FALSE)
    if (nrow(wide) != nrow(target)) stop(target_cohort, ": source join multiplied frozen scans")
    wide[, source_matched := !is.na(source_row)]
    if (!"source_sheet" %in% names(wide)) wide[, source_sheet := NA_character_]
    status <- wide[, .(cohort, group, id, mri_date, image_id, time_years,
                       source_matched, source_sheet, source_row, scanner, processing_version)]
    id_vars <- c("cohort", "group", "id", "mri_date", "image_id", "time_years",
                 "tau_z", "amyloid_z", "age10", "sex_female", "education_c",
                 "apoe4_carrier", "diagnosis", "scanner", "processing_version",
                 "source_matched")
    x <- melt(wide, id.vars = id_vars, measure.vars = fields,
              variable.name = "field", value.name = "thickness_mm", variable.factor = FALSE)
    x <- merge(x, source_map, by = "field", all.x = TRUE, sort = FALSE)
    list(roi = x, status = status)
  }
  adni_map <- map[, .(field = adni_field, roi, side, band)]
  habs_map <- map[, .(field = habs_field, roi, side, band)]
  adni <- make_roi("ADNI", a, adni_join, adni_fields, adni_map)
  habs <- make_roi("HABS-HD", h, habs_join, habs_fields, habs_map)
  roi <- rbindlist(list(adni$roi, habs$roi), use.names = TRUE)
  if (nrow(roi) != nrow(scans) * 68L) stop("Each scan must have exactly 68 regional rows")
  roi[, valid := source_matched == TRUE & is.finite(thickness_mm) &
                 thickness_mm >= 1 & thickness_mm <= 5]
  roi[, y := NA_real_]
  roi[valid == TRUE, y := 100 * log(thickness_mm)]
  status <- rbindlist(list(adni$status, habs$status), use.names = TRUE)
  save_tsv(status, "Scan_to_regional_source_match.tsv", "audit")
  scan_counts <- roi[, .(valid_roi = sum(valid == TRUE)),
                     by = .(cohort, group, id, mri_date, image_id)]
  scan_counts <- merge(status, scan_counts, by = c("cohort", "group", "id", "mri_date", "image_id"))
  save_tsv(scan_counts[, .(people = uniqueN(id), scans = .N,
                           source_matched_scans = sum(source_matched == TRUE),
                           complete_68_scans = sum(valid_roi == 68L),
                           median_valid_roi = as.numeric(median(valid_roi))), by = .(cohort, group)],
           "Table3_source_match_and_68_region_completeness.tsv")
  save_tsv(roi[, .(valid_scans = sum(valid == TRUE),
                   valid_people = uniqueN(id[valid == TRUE]),
                   longitudinal_people = uniqueN(id[valid == TRUE & time_years > 0 &
                       id %in% id[valid == TRUE & time_years == 0]])),
                 by = .(cohort, group, side, roi, band, field)],
           "Table4_per_region_measurement_coverage.tsv")
  mtl <- roi[roi %in% c("entorhinal", "parahippocampal", "inferiortemporal", "middletemporal"),
             .(source_mtl_mm = if (.N == 8L && all(valid == TRUE)) mean(thickness_mm) else NA_real_),
             by = .(cohort, group, id, mri_date, image_id)]
  mtl <- merge(mtl, scans[, .(cohort, group, id, mri_date, image_id, frozen_cortex)],
               by = c("cohort", "group", "id", "mri_date", "image_id"))
  mtl <- mtl[group == "Apositive"]
  mtl[, absolute_difference := abs(source_mtl_mm - frozen_cortex)]
  save_tsv(mtl, "Frozen_eight_field_cortical_reproduction.tsv", "audit")
  if (any(!is.finite(mtl$absolute_difference)) || any(mtl$absolute_difference > 1e-5))
    stop("Regional-source extraction does not reproduce the frozen eight-field cortex measure")
  log_line("Source alignment complete; fitting parcel and paired-band models")

  fit_one <- function(d, outcome, cohort, group, feature, feature_type,
                      random = "intercept") {
    covariates <- c("sex_female", "education_c", "apoe4_carrier", "diagnosis",
                    "scanner", "processing_version")
    needed <- c("id", "time_years", "tau_z", "amyloid_z", "age10", outcome, covariates)
    d <- copy(d[complete.cases(d[, ..needed])])
    eligible <- d[, .(n_scans = .N, baseline = any(time_years == 0)), by = id]
    d <- d[id %in% eligible[baseline == TRUE & n_scans >= 2L, id]]
    covariates <- covariates[vapply(covariates, function(v) uniqueN(d[[v]]) > 1L, logical(1))]
    empty <- data.table(cohort, group, feature, feature_type, random,
                        term = c("time_years:tau_z", "time_years:age10"),
                        estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
                        ci_high = NA_real_, p = NA_real_, n_people = uniqueN(d$id),
                        n_scans = nrow(d), singular = NA, status = "insufficient_data")
    if (uniqueN(d$id) < 20L || nrow(d) < 40L) return(empty)
    for (v in covariates) if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d, j = v, value = droplevels(factor(d[[v]])))
    random_term <- if (random == "slope") "(1 + time_years | id)" else "(1 | id)"
    formula <- as.formula(paste(outcome, "~",
      paste(c("time_years * tau_z", "time_years * age10", "time_years * amyloid_z",
              covariates, random_term), collapse = " + ")))
    fit <- tryCatch(suppressWarnings(lmerTest::lmer(formula, data = d, REML = FALSE,
                    control = lme4::lmerControl(optimizer = "bobyqa",
                      optCtrl = list(maxfun = 100000)))), error = identity)
    if (inherits(fit, "error")) {
      empty[, status := paste0("model_error: ", conditionMessage(fit))]
      return(empty)
    }
    coefficients <- summary(fit)$coefficients
    for (i in seq_len(nrow(empty))) {
      term <- empty$term[i]
      if (!term %in% rownames(coefficients)) next
      beta <- unname(coefficients[term, "Estimate"])
      se <- unname(coefficients[term, "Std. Error"])
      p <- if ("Pr(>|t|)" %in% colnames(coefficients))
        unname(coefficients[term, "Pr(>|t|)"]) else 2 * pnorm(-abs(beta / se))
      set(empty, i = i, j = c("estimate", "se", "ci_low", "ci_high", "p", "singular", "status"),
          value = list(beta, se, beta - 1.96 * se, beta + 1.96 * se, p,
                       lme4::isSingular(fit, tol = 1e-4), "fitted"))
    }
    empty
  }
  all_fits <- list()
  index <- 0L
  for (co in c("ADNI", "HABS-HD")) for (gr in c("Apositive", "Anegative")) {
    fields <- if (co == "ADNI") adni_fields else habs_fields
    for (selected_field in fields) {
      index <- index + 1L
      d <- roi[cohort == co & group == gr & field == selected_field & valid == TRUE]
      all_fits[[index]] <- fit_one(d, "y", co, gr, selected_field, "cortical_region")
      if (index %% 34L == 0L) log_line("Cortical regions modeled: ", index, "/272")
    }
  }
  effects <- rbindlist(all_fits)
  effects <- merge(effects, rbindlist(list(
      adni_map[, .(cohort = "ADNI", field, roi, side, band)],
      habs_map[, .(cohort = "HABS-HD", field, roi, side, band)])),
      by.x = c("cohort", "feature"), by.y = c("cohort", "field"), all.x = TRUE)
  effects[, q_bh := ifelse(is.finite(p), p.adjust(p, method = "BH"), NA_real_),
          by = .(cohort, group, term)]
  save_tsv(effects, "Table5_68_region_tau_and_age_effects.tsv")

  expected <- map[, .(expected_roi = .N), by = band]
  band <- roi[valid == TRUE, .(n_roi = .N, y = mean(y)),
              by = .(cohort, group, id, mri_date, image_id, band)]
  band <- merge(band, expected, by = "band")
  band <- band[n_roi == expected_roi]
  bandwide <- dcast(band, cohort + group + id + mri_date + image_id ~ band, value.var = "y")
  scan_covars <- scans[, .(cohort, group, id, mri_date, image_id, time_years,
                           stable_pipeline, tau_z, amyloid_z, age10, sex_female, education_c,
                           apoe4_carrier, diagnosis, scanner, processing_version)]
  bandwide <- merge(scan_covars, bandwide,
                    by = c("cohort", "group", "id", "mri_date", "image_id"), all.x = TRUE)
  save_tsv(bandwide[group == "Apositive"], "Paired_band_scan_data.tsv", "source_data")
  contrasts <- list(
    "Medial minus other temporal" = c("Medial temporal", "Other temporal"),
    "Medial minus non-temporal" = c("Medial temporal", "Non-temporal"),
    "Other temporal minus non-temporal" = c("Other temporal", "Non-temporal"))
  band_results <- list()
  index <- 0L
  for (co in c("ADNI", "HABS-HD")) for (gr in c("Apositive", "Anegative"))
    for (name in names(contrasts)) {
      index <- index + 1L
      pair <- contrasts[[name]]
      d <- copy(bandwide[bandwide$cohort == co & bandwide$group == gr])
      d[, y := get(pair[1L]) - get(pair[2L])]
      band_results[[index]] <- fit_one(d, "y", co, gr, name, "paired_band_contrast")
    }
  band_effects <- rbindlist(band_results)
  save_tsv(band_effects, "Table6_direct_paired_band_differences.tsv")

  sensitivity <- list()
  index <- 0L
  for (selected_field in adni_fields) {
    index <- index + 1L
    sensitivity[[index]] <- fit_one(roi[cohort == "ADNI" & group == "Apositive" &
                                      field == selected_field & valid == TRUE],
                                    "y", "ADNI", "Apositive", selected_field,
                                    "cortical_region", random = "slope")
  }
  for (name in names(contrasts)) {
    index <- index + 1L
    pair <- contrasts[[name]]
    d <- copy(bandwide[cohort == "ADNI" & group == "Apositive"])
    d[, y := get(pair[1L]) - get(pair[2L])]
    sensitivity[[index]] <- fit_one(d, "y", "ADNI", "Apositive", name,
                                    "paired_band_contrast", random = "slope")
  }
  save_tsv(rbindlist(sensitivity), "Table7_ADNI_random_slope_sensitivity.tsv")
  cn <- list()
  index <- 0L
  for (co in c("ADNI", "HABS-HD")) for (name in names(contrasts)) {
    index <- index + 1L
    pair <- contrasts[[name]]
    d <- copy(bandwide[cohort == co & group == "Anegative" & diagnosis == "CN"])
    d[, y := get(pair[1L]) - get(pair[2L])]
    cn[[index]] <- fit_one(d, "y", co, "Anegative_CN", name, "paired_band_contrast")
  }
  save_tsv(rbindlist(cn), "Table8_amyloid_negative_CN_age_reference_sensitivity.tsv")

  map_effects <- effects[(group == "Apositive" & term == "time_years:tau_z") |
                         (group == "Anegative" & term == "time_years:age10")]
  map_effects[, signature := fifelse(group == "Apositive", "Tau in Aβ+", "Age in Aβ−")]
  save_tsv(map_effects, "Figure34a_region_effect_source.tsv", "source_data")
  correlations <- rbindlist(lapply(unique(map_effects$signature), function(label) {
    z <- dcast(map_effects[signature == label], roi + side + band ~ cohort,
               value.var = "estimate")
    z <- z[is.finite(ADNI) & is.finite(`HABS-HD`)]
    data.table(signature = label, n_regions = nrow(z),
               pearson = if (nrow(z) > 2L) cor(z$ADNI, z$`HABS-HD`) else NA_real_,
               spearman = if (nrow(z) > 2L) cor(z$ADNI, z$`HABS-HD`, method = "spearman") else NA_real_)
  }))
  save_tsv(correlations, "Table9_cross_cohort_spatial_correlations_descriptive.tsv")
  within_cohort <- rbindlist(lapply(c("ADNI", "HABS-HD"), function(co) {
    tau <- effects[cohort == co & group == "Apositive" & term == "time_years:tau_z",
                   .(roi, side, tau = estimate)]
    age <- effects[cohort == co & group == "Anegative" & term == "time_years:age10",
                   .(roi, side, age = estimate)]
    x <- merge(tau, age, by = c("roi", "side"))[is.finite(tau) & is.finite(age)]
    data.table(cohort = co, n_regions = nrow(x),
               pearson = if (nrow(x) > 2L) cor(x$tau, x$age) else NA_real_,
               spearman = if (nrow(x) > 2L) cor(x$tau, x$age, method = "spearman") else NA_real_)
  }))
  save_tsv(within_cohort, "Table10_within_cohort_age_tau_map_correlations_descriptive.tsv")
  save_tsv(fread(files[["hippocampus"]], encoding = "UTF-8")[
             scope %in% c("ADNI_primary", "HABS_primary") & metric == "log_hip_ratio"],
           "Table11_frozen_whole_hippocampus_context.tsv")

  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  log_line("Analysis complete: ", out)
  print(scan_counts[, .(people = uniqueN(id), scans = .N,
                        source_matched = sum(source_matched == TRUE),
                        complete_68 = sum(valid_roi == 68L)), by = .(cohort, group)])
  print(band_effects[term %in% c("time_years:tau_z", "time_years:age10") &
                       feature == "Medial minus non-temporal",
                     .(cohort, group, term, estimate, ci_low, ci_high, p, n_people, n_scans)])
  print(correlations)
})
