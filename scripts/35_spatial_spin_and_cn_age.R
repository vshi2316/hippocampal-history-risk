local({
  needed <- c("data.table", "readxl", "lme4", "lmerTest",
              "freesurferformats")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse = ", "))
  library(data.table)
  root <- .study$output
  work <- .study$work
  atlas_dir <- file.path(work, "inputs", "ENIGMA_BigBrain")
  out <- file.path(root, "bigbrain_spin_60",
                   paste0("run_", "frozen"))
  for (d in c("tables", "audit"))
    dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  paths <- c(
    matched = file.path(work, "outputs", "normative_input_audit_59",
                        "Apositive_Tau_effects_matched_to_BigBrain.tsv"),
    reference = file.path(work, "outputs", "normative_input_audit_59",
                          "DK68_BigBrain_reference.tsv"),
    long = file.path(root, "tau_longitudinal_inputs/run_frozen/analysis_data",
                     "MRI_longitudinal_analysis.tsv"),
    dictionary = file.path(root, "roi68_spatial_34/run_frozen/tables",
                           "Table1_68_region_field_and_band_dictionary.tsv"),
    adni = file.path(.study$data,"adni/MRI/MRI_FreeSurfer.csv"),
    habs = file.path(.study$data,"habs_hd/HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx")
  )
  if (any(!file.exists(paths)))
    stop("Missing input: ", paste(paths[!file.exists(paths)], collapse = "; "))

  # Download only the public atlas sphere and annotation assets, at a fixed commit.
  commit <- "b08974b55243060cbc1fad12c87048037446e8f7"
  atlas_paths <- c(
    L_sphere = "matlab/shared/surfaces/fsa5_sphere_lh",
    R_sphere = "matlab/shared/surfaces/fsa5_sphere_rh",
    L_annot = "enigmatoolbox/permutation_testing/annot/fsa5_lh_aparc.annot",
    R_annot = "enigmatoolbox/permutation_testing/annot/fsa5_rh_aparc.annot"
  )
  dir.create(atlas_dir, recursive = TRUE, showWarnings = FALSE)
  atlas_files <- file.path(atlas_dir, basename(atlas_paths))
  names(atlas_files) <- names(atlas_paths)
  for (key in names(atlas_paths)) if (!file.exists(atlas_files[[key]])) {
    url <- paste0("https://raw.githubusercontent.com/MICA-MNI/ENIGMA/",
                  commit, "/", atlas_paths[[key]])
    utils::download.file(url, atlas_files[[key]], mode = "wb", quiet = TRUE)
  }
  write_tsv(data.table(file = c(paths, atlas_files),
                       md5 = unname(tools::md5sum(c(paths, atlas_files)))),
            "Input_manifest.tsv", "audit")

  reference <- fread(paths[["reference"]])
  tau <- fread(paths[["matched"]])
  setorder(reference, side, roi)
  setorder(tau, cohort, side, roi)
  stopifnot(nrow(reference) == 68L, nrow(tau) == 136L,
            all(reference$side == rep(c("L", "R"), each = 34L)))
  for (co in c("ADNI", "HABS-HD")) {
    d <- tau[cohort == co]
    if (nrow(d) != 68L || !identical(d[, .(side, roi)], reference[, .(side, roi)]) ||
        any(!is.finite(d$estimate))) stop("Stage 59 match changed: ", co)
  }

  clean_roi <- function(x) tolower(gsub("[^[:alnum:]]", "", x))
  centers <- list()
  for (hemisphere in c("L", "R")) {
    surf <- freesurferformats::read.fs.surface(
      atlas_files[[paste0(hemisphere, "_sphere")]], format = "bin")
    annot <- freesurferformats::read.fs.annot(atlas_files[[paste0(hemisphere, "_annot")]])
    xyz <- surf$vertices
    if (nrow(xyz) != length(annot$label_names)) stop("Atlas vertex count differs: ", hemisphere)
    roi_names <- reference[reference$side == hemisphere, roi]
    lab <- clean_roi(annot$label_names)
    if (anyDuplicated(clean_roi(roi_names))) stop("Duplicate DK names")
    z <- t(vapply(roi_names, function(region) {
      ix <- which(lab == clean_roi(region))
      if (!length(ix)) stop("Atlas annotation lacks: ", hemisphere, "_", region)
      v <- colMeans(xyz[ix, , drop = FALSE])
      v / sqrt(sum(v^2))
    }, numeric(3)))
    centers[[hemisphere]] <- z
  }
  distance_matrix <- function(a, b) {
    d2 <- outer(rowSums(a^2), rowSums(b^2), "+") - 2 * a %*% t(b)
    sqrt(pmax(d2, 0))
  }
  greedy_match <- function(d) {
    n <- nrow(d)
    ans <- integer(n)
    while (any(ans == 0L)) {
      available <- which(ans == 0L)
      row_min <- apply(d[available, , drop = FALSE], 1L, min)
      i <- available[which.max(row_min)]
      j <- which.min(d[i, ])
      ans[i] <- j
      d[, j] <- Inf
    }
    ans
  }
  rotation <- function() {
    x <- matrix(rnorm(9L), 3L, 3L)
    decomposition <- qr(x)
    q <- qr.Q(decomposition) %*% diag(sign(diag(qr.R(decomposition))))
    if (det(q) < 0) q[, 1L] <- -q[, 1L]
    q
  }
  set.seed(20261118)
  n_rot <- 2000L
  perm <- matrix(NA_integer_, 68L, n_rot)
  reflect <- diag(c(-1, 1, 1))
  k <- 0L
  while (k < n_rot) {
    rot_l <- rotation()
    rot_r <- reflect %*% rot_l %*% reflect
    ix_l <- greedy_match(distance_matrix(centers$L, centers$L %*% rot_l))
    ix_r <- greedy_match(distance_matrix(centers$R, centers$R %*% rot_r))
    ix <- c(ix_l, ix_r + 34L)
    if (all(ix == seq_len(68L))) next
    k <- k + 1L
    perm[, k] <- ix
  }
  saveRDS(perm, file.path(out, "audit", "Spin_permutations_2000.rds"))
  write_tsv(data.table(side = reference$side, roi = reference$roi,
                       x = c(centers$L[, 1], centers$R[, 1]),
                       y = c(centers$L[, 2], centers$R[, 2]),
                       z = c(centers$L[, 3], centers$R[, 3])),
            "DK68_spherical_centroids.tsv", "audit")

  # Rebuild the Stage 34 non-training Anegative scan set, restricted to baseline CN.
  scans <- fread(paths[["long"]])
  scans <- scans[amyloid_pos == 0L & reference_partition != "training" & diagnosis == "CN"]
  scans[, `:=`(id = as.character(id), image_id = as.character(image_id),
               mri_date = as.Date(mri_date), age10 = age_c)]
  dictionary <- fread(paths[["dictionary"]])
  adni_fields <- dictionary$adni_field
  habs_fields <- dictionary$habs_field
  a <- fread(paths[["adni"]], select = c("RID", "EXAMDATE", "IMAGEUID",
             "FIELD_STRENGTH", "FSVER", adni_fields))
  a[, `:=`(id = paste0("ADNI_", sub("\\.0+$", "", as.character(RID))),
            mri_date = as.Date(EXAMDATE), image_id = as.character(IMAGEUID),
            scanner = as.character(FIELD_STRENGTH),
            processing_version = as.character(FSVER), source_row = .I)]
  a <- a[id %in% scans[cohort == "ADNI", id]]
  a <- a[, c("id", "mri_date", "image_id", "scanner", "processing_version",
             "source_row", adni_fields), with = FALSE]
  if (anyDuplicated(a, by = c("id", "mri_date", "image_id", "scanner", "processing_version")))
    stop("Duplicate ADNI source scan key")
  month_date <- function(year, month) {
    y <- suppressWarnings(as.integer(year))
    m <- suppressWarnings(as.integer(month))
    missing <- which(is.na(m) & !is.na(month))
    if (length(missing)) m[missing] <- match(tolower(substr(month[missing], 1, 3)),
                                            tolower(month.abb))
    as.Date(ifelse(is.finite(y) & is.finite(m) & m >= 1 & m <= 12,
                   sprintf("%04d-%02d-15", y, m), NA_character_))
  }
  h <- rbindlist(lapply(paste("Inclusions - Scan", 1:4), function(sheet) {
    x <- as.data.table(readxl::read_excel(paths[["habs"]], sheet = sheet,
                      col_types = "text", na = c("", "NA", "N/A"),
                      .name_repair = "unique_quiet"))
    x <- x[!is.na(Med_ID)]
    x[, source_row := .I + 1L]
    x
  }), use.names = TRUE, fill = TRUE)
  h[, `:=`(id = paste0("HABS_", sub("\\.0+$", "", as.character(Med_ID))),
            mri_date = month_date(MRI_year, MRI_month),
            image_id = as.character(Time_point_MRI),
            scanner = as.character(Scanner_MRI))]
  h <- h[id %in% scans[cohort == "HABS-HD", id]]
  h <- h[, c("id", "mri_date", "image_id", "scanner", "source_row", habs_fields),
         with = FALSE]
  if (anyDuplicated(h, by = c("id", "mri_date", "image_id", "scanner")))
    stop("Duplicate HABS-HD source scan key")
  for (field in adni_fields) set(a, j = field, value = suppressWarnings(as.numeric(a[[field]])))
  for (field in habs_fields) set(h, j = field, value = suppressWarnings(as.numeric(h[[field]])))

  region_data <- function(co, source, join, fields, field_col) {
    target <- copy(scans[cohort == co])
    wide <- merge(target, source, by = join, all.x = TRUE, sort = FALSE)
    if (nrow(wide) != nrow(target)) stop("Source join multiplied scans: ", co)
    id_vars <- c("cohort", "id", "mri_date", "image_id", "time_years",
                 "tau_z", "amyloid_z", "age10", "sex_female", "education_c",
                 "apoe4_carrier", "scanner", "processing_version", "source_row")
    long <- melt(wide, id.vars = id_vars, measure.vars = fields,
                 variable.name = "field", value.name = "thickness_mm",
                 variable.factor = FALSE)
    map <- dictionary[, .(field = get(field_col), side, roi)]
    merge(long, map, by = "field", all.x = TRUE, sort = FALSE)
  }
  cn <- rbindlist(list(
    region_data("ADNI", a,
                c("id", "mri_date", "image_id", "scanner", "processing_version"),
                adni_fields, "adni_field"),
    region_data("HABS-HD", h,
                c("id", "mri_date", "image_id", "scanner"),
                habs_fields, "habs_field")
  ), use.names = TRUE, fill = TRUE)
  cn[, valid := !is.na(source_row) & is.finite(thickness_mm) &
       thickness_mm >= 1 & thickness_mm <= 5]
  cn[valid == TRUE, y := 100 * log(thickness_mm)]
  write_tsv(cn[, .(people = uniqueN(id), scans = uniqueN(paste(id, mri_date, image_id)),
                   matched_region_rows = sum(valid)), by = cohort],
            "CN_source_coverage.tsv", "audit")

  fit_age <- function(d) {
    covars <- c("sex_female", "education_c", "apoe4_carrier", "scanner",
                "processing_version")
    require <- c("id", "time_years", "tau_z", "amyloid_z", "age10", "y", covars)
    d <- copy(d[complete.cases(d[, ..require])])
    eligible <- d[, .(n_scans = .N, baseline = any(time_years == 0)), by = id]
    d <- d[id %in% eligible[n_scans >= 2L & baseline == TRUE, id]]
    answer <- data.table(n_people = uniqueN(d$id), n_scans = nrow(d),
                         estimate = NA_real_, se = NA_real_,
                         ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
                         singular = NA, status = "insufficient_data")
    if (answer$n_people < 20L || answer$n_scans < 40L) return(answer)
    covars <- covars[vapply(covars, function(x) uniqueN(d[[x]]) > 1L, logical(1))]
    for (v in covars) if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d, j = v, value = droplevels(factor(d[[v]])))
    form <- as.formula(paste("y ~ time_years * tau_z + time_years * age10 +",
      "time_years * amyloid_z +", paste(c(covars, "(1 | id)"), collapse = " + ")))
    fit <- tryCatch(suppressWarnings(lmerTest::lmer(form, d, REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 100000)))), error = identity)
    if (inherits(fit, "error")) {
      answer[, status := paste0("model_error: ", conditionMessage(fit))]
      return(answer)
    }
    z <- summary(fit)$coefficients
    term <- "time_years:age10"
    if (!term %in% rownames(z)) return(answer)
    answer[, `:=`(estimate = z[term, "Estimate"], se = z[term, "Std. Error"],
                  ci_low = z[term, "Estimate"] - 1.96 * z[term, "Std. Error"],
                  ci_high = z[term, "Estimate"] + 1.96 * z[term, "Std. Error"],
                  p = z[term, "Pr(>|t|)"], singular = lme4::isSingular(fit),
                  status = "fitted")]
    answer
  }
  keys <- unique(cn[, .(cohort, side, roi, field)])
  age <- rbindlist(lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i]
    d <- cn[cohort == key$cohort & field == key$field & valid == TRUE]
    cbind(key, fit_age(d))
  }), use.names = TRUE)
  setorder(age, cohort, side, roi)
  write_tsv(age, "CN_68_region_age_effects.tsv")
  if (nrow(age) != 136L || any(age$status != "fitted") ||
      any(!is.finite(age$estimate)))
    stop("CN regional fits incomplete; inspect CN_68_region_age_effects.tsv")

  maps <- rbindlist(list(
    tau[, .(cohort, side, roi, signature = "Tau in Abeta+",
            estimate, se, n_people, n_scans)],
    age[, .(cohort, side, roi, signature = "Age in Abeta- CN",
            estimate, se, n_people, n_scans)]
  ))
  maps <- merge(maps, reference, by = c("side", "roi"), all.x = TRUE, sort = FALSE)
  setorder(maps, cohort, signature, side, roi)
  write_tsv(maps, "Figure60_region_source_data.tsv")
  vec <- function(co, sig_name) {
    z <- maps[cohort == co & signature == sig_name]
    setorder(z, side, roi)
    if (nrow(z) != 68L || anyNA(z$bigbrain_gradient)) stop("Map mismatch")
    z
  }
  ref <- reference$bigbrain_gradient
  spin_cor <- function(y, type = "pearson") {
    obs <- cor(ref, y, method = type)
    null <- apply(perm, 2L, function(ix) cor(ref[ix], y, method = type))
    list(r = obs, p_spin = (1 + sum(abs(null) >= abs(obs))) / (length(null) + 1),
         null = null)
  }
  results <- list()
  nulls <- list()
  contrasts <- list()
  k <- 0L
  for (co in c("ADNI", "HABS-HD")) {
    tau_map <- vec(co, "Tau in Abeta+")
    age_map <- vec(co, "Age in Abeta- CN")
    for (sig in c("Tau in Abeta+", "Age in Abeta- CN")) {
      z <- if (sig == "Tau in Abeta+") tau_map else age_map
      k <- k + 1L
      test <- spin_cor(z$estimate)
      results[[k]] <- data.table(cohort = co, signature = sig, regions = 68L,
        people_min = min(z$n_people), people_max = max(z$n_people),
        scans_min = min(z$n_scans), scans_max = max(z$n_scans),
        r_pearson = test$r, r_spearman = cor(ref, z$estimate, method = "spearman"),
        p_spin_two_sided = test$p_spin)
      nulls[[k]] <- data.table(cohort = co, signature = sig,
                               rotation = seq_along(test$null), r_null = test$null)
    }
    obs_diff <- cor(ref, tau_map$estimate) - cor(ref, age_map$estimate)
    null_diff <- vapply(seq_len(n_rot), function(i)
      cor(ref[perm[, i]], tau_map$estimate) -
      cor(ref[perm[, i]], age_map$estimate), numeric(1))
    contrasts[[co]] <- data.table(cohort = co, delta_r_tau_minus_age = obs_diff,
      p_spin_two_sided = (1 + sum(abs(null_diff) >= abs(obs_diff))) / (n_rot + 1L))
  }
  results <- rbindlist(results)
  write_tsv(results, "Table1_BigBrain_spatial_tests.tsv")
  write_tsv(rbindlist(contrasts), "Table2_Tau_vs_CN_age_spatial_contrast.tsv")
  write_tsv(rbindlist(nulls), "Spin_null_distributions.tsv", "audit")
  adni_tau <- vec("ADNI", "Tau in Abeta+")$estimate
  habs_tau <- vec("HABS-HD", "Tau in Abeta+")$estimate
  observed_cross <- cor(adni_tau, habs_tau)
  cross_null <- apply(perm, 2L, function(ix) cor(adni_tau[ix], habs_tau))
  write_tsv(data.table(r_ADNI_HABS = observed_cross,
                       p_spin_two_sided = (1 + sum(abs(cross_null) >=
                                              abs(observed_cross))) / (n_rot + 1L)),
            "Table3_Tau_cross_cohort_spatial_alignment.tsv")

  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  cat("Output: ", out, "\n", sep = "")
  print(results)
  print(rbindlist(contrasts))
})
