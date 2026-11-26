local({
  packages <- c("data.table", "lme4", "lmerTest")
  for (pkg in packages)
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Install R package: ", pkg)
  suppressPackageStartupMessages(library(data.table))

  root <- .study$output
  inputs <- c(
    frozen = file.path(root, "tau_longitudinal_inputs/run_frozen/models/Longitudinal_models.rds"),
    ashs = file.path(.study$data,"adni/MRI/MRI_Hippocampus_ASHS.csv")
  )
  if (any(!file.exists(inputs)))
    stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  out <- file.path(root, "paired_segmentation",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name),
           sep = "\t", quote = TRUE, na = "NA", bom = TRUE)
  save_tsv(data.table(input = names(inputs), file = unname(inputs),
                      md5 = unname(tools::md5sum(inputs))),
           "Input_manifest.tsv", "audit")

  key <- "ADNI__Apositive_primary__cortex__tau_time"
  x <- copy(as.data.table(readRDS(inputs[["frozen"]])$MRI[[key]]$data))
  x[, `:=`(id = as.character(id), image_id = as.character(image_id),
           mri_date = as.Date(mri_date))]
  stopifnot(nrow(x) == 277L, uniqueN(x$id) == 94L,
            !anyDuplicated(x[, .(id, image_id)]))
  a <- fread(inputs[["ashs"]], encoding = "UTF-8",
    select = c("RID", "EXAMDATE", "IMAGEUID_T1", "IMAGEUID_T2",
               "LEFT_HIPP_VOL", "RIGHT_HIPP_VOL", "T2_QC_L", "T2_QC_R",
               "T2_QC_COMMENT"))
  a[, `:=`(id = paste0("ADNI_", as.character(RID)),
           image_id = as.character(IMAGEUID_T1),
           ashs_date = as.Date(EXAMDATE),
           t2_image_id = as.character(IMAGEUID_T2),
           qc_left = as.character(T2_QC_L),
           qc_right = as.character(T2_QC_R),
           qc_comment = as.character(T2_QC_COMMENT),
           ashs_mm3 = (as.numeric(LEFT_HIPP_VOL) + as.numeric(RIGHT_HIPP_VOL)) / 2)]
  a <- a[image_id %in% x$image_id,
         .(id, image_id, ashs_date, t2_image_id, qc_left, qc_right,
           qc_comment, ashs_mm3)]
  if (anyDuplicated(a[, .(id, image_id)])) stop("Duplicate ASHS T1 image ID")
  d <- merge(x, a, by = c("id", "image_id"), all.x = TRUE, sort = FALSE)
  stopifnot(nrow(d) == nrow(x))
  d <- d[!is.na(ashs_date) & ashs_date == mri_date &
           is.finite(ashs_mm3) & ashs_mm3 > 0]
  d <- d[id %in% d[time_years == 0, id]]
  setorder(d, id, mri_date)
  d[, `:=`(fs_pct = 100 * log_hip_ratio,
           ashs_pct = 100 * log(ashs_mm3 / baseline_icv_mm3),
           fs_mm3 = as.numeric(hip_mean_mm3))]
  d[, method_diff_pct := fs_pct - ashs_pct]
  stopifnot(nrow(d) == 239L, uniqueN(d$id) == 85L,
            !anyDuplicated(d[, .(id, mri_date)]),
            all(is.finite(d$fs_pct)), all(is.finite(d$ashs_pct)))

  ids3 <- d[, .N, by = id][N >= 3, id]
  three <- d[id %in% ids3]
  agreement <- three[, {
    fs_log <- log(fs_mm3)
    ashs_log <- log(ashs_mm3)
    t <- time_years
    .(fs_baseline = fs_mm3[1L] / baseline_icv_mm3[1L],
      ashs_baseline = ashs_mm3[1L] / baseline_icv_mm3[1L],
      fs_two = 100 * (fs_log[2L] - fs_log[1L]) / (t[2L] - t[1L]),
      ashs_two = 100 * (ashs_log[2L] - ashs_log[1L]) / (t[2L] - t[1L]),
      fs_multi = unname(coef(lm(fs_log ~ t))[2L]) * 100,
      ashs_multi = unname(coef(lm(ashs_log ~ t))[2L]) * 100,
      n_scans = .N, followup_years = max(t) - min(t))
  }, by = id]
  save_tsv(agreement, "Individual_same_person_measurement_agreement.tsv", "source_data")
  agreement_stats <- rbindlist(lapply(list(
    baseline = c("fs_baseline", "ashs_baseline"),
    first_two = c("fs_two", "ashs_two"),
    three_plus = c("fs_multi", "ashs_multi")), function(cols) {
    z <- agreement[complete.cases(agreement[, ..cols])]
    data.table(n_people = nrow(z), pearson_r = cor(z[[cols[1L]]], z[[cols[2L]]]),
               spearman_r = cor(z[[cols[1L]]], z[[cols[2L]]], method = "spearman"))
  }), idcol = "level")
  save_tsv(agreement_stats, "Table6_same_person_agreement_levels.tsv")

  # Codes remain categorical until their released QC definitions are verified.
  d[, `:=`(previous_image_id = shift(image_id),
           previous_t2_id = shift(t2_image_id),
           previous_date = shift(mri_date),
           interval_years = time_years - shift(time_years),
           fs_rate = (fs_pct - shift(fs_pct)) / (time_years - shift(time_years)),
           ashs_rate = (ashs_pct - shift(ashs_pct)) / (time_years - shift(time_years))),
    by = id]
  d[, `:=`(method_gap = abs(fs_rate - ashs_rate),
           ashs_jump15 = is.finite(ashs_rate) & abs(ashs_rate) > 15)]
  d[is.na(ashs_jump15), ashs_jump15 := FALSE]
  shock_ids <- unique(d[ashs_jump15 == TRUE, id])
  qc <- d[, .(scans = .N, people = uniqueN(id),
              transition_scans = sum(is.finite(ashs_rate)),
              ashs_jump15_scans = sum(ashs_jump15)),
          by = .(qc_left, qc_right)][order(qc_left, qc_right)]
  save_tsv(qc, "Table1_ASHS_QC_codes_unranked.tsv")
  q4 <- d[qc_left == "4" & qc_right == "4"]
  q4_ids <- q4[, .(scans = .N, baseline = any(time_years == 0)), by = id]
  save_tsv(data.table(qc_pair = "4/4, descriptive only", scans = nrow(q4),
                      people = uniqueN(q4$id),
                      baseline_plus_repeat = sum(q4_ids$baseline & q4_ids$scans >= 2L)),
           "QC_4_4_feasibility.tsv", "audit")
  transitions <- d[is.finite(ashs_rate),
    .(id, previous_image_id, image_id, previous_t2_id, t2_image_id,
      previous_date, mri_date, interval_years, fs_rate, ashs_rate,
      method_gap, ashs_jump15, qc_left, qc_right, qc_comment)]
  setorder(transitions, -method_gap)
  save_tsv(transitions, "Table2_all_scan_transitions.tsv")
  review_ids <- unique(c("ADNI_6634", "ADNI_6480", shock_ids))
  review <- d[id %in% review_ids,
    .(id, image_id, t2_image_id, mri_date, previous_image_id,
      previous_t2_id, previous_date, fs_mm3, ashs_mm3, fs_rate, ashs_rate,
      method_gap, qc_left, qc_right, qc_comment)]
  setorder(review, id, mri_date)
  save_tsv(review, "Targeted_T1_T2_segmentation_review_manifest.tsv", "audit")
  save_tsv(data.table(rule = "Any observed ASHS interval with absolute annual log change >15%",
                      id = shock_ids), "Exploratory_jump_flagged_people.tsv", "audit")

  covariates <- c("sex_female", "education_c", "apoe4_carrier",
                  "diagnosis", "scanner", "processing_version")
  covariates <- covariates[vapply(covariates, function(v)
    v %in% names(d) && uniqueN(d[[v]][!is.na(d[[v]])]) > 1L, logical(1))]
  needed <- c("id", "time_years", "tau_z", "amyloid_z", "age_c",
              "fs_pct", "ashs_pct", "method_diff_pct", covariates)
  if (!all(complete.cases(d[, ..needed]))) stop("Unexpected missing model value in frozen matched scans")
  for (v in covariates)
    if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d, j = v, value = factor(d[[v]]))
  save_tsv(data.table(variable = covariates,
                      levels = vapply(covariates, function(v)
                        paste(unique(as.character(d[[v]])), collapse = " | "), character(1))),
           "Common_adjustment_set.tsv", "audit")

  make_formula <- function(y, group, random, covars) {
    re <- if (random == "random_slope")
      paste0("(1 + time_years | ", group, ")") else paste0("(1 | ", group, ")")
    as.formula(paste(y, "~", paste(c("time_years*tau_z",
      "time_years*amyloid_z", "time_years*age_c", covars, "time_years", re),
      collapse = " + ")))
  }
  fit_one <- function(dat, y, scope, method, random, covars) {
    formula <- make_formula(y, "id", random, covars)
    warnings <- character()
    model <- tryCatch(withCallingHandlers(
      lmerTest::lmer(formula, data = dat, REML = FALSE,
        control = lme4::lmerControl(optimizer = "bobyqa",
                                    optCtrl = list(maxfun = 100000))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(model, "error"))
      return(data.table(scope, method, random, outcome = y,
        n_people = uniqueN(dat$id), n_scans = nrow(dat),
        estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
        p = NA_real_, singular = NA, status = conditionMessage(model),
        warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
        formula = paste(deparse(formula), collapse = " ")))
    tab <- summary(model)$coefficients
    term <- "time_years:tau_z"
    beta <- unname(tab[term, "Estimate"])
    se <- unname(tab[term, "Std. Error"])
    data.table(scope, method, random, outcome = y,
      n_people = uniqueN(dat$id), n_scans = nrow(dat),
      estimate = beta, se = se,
      ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se,
      p = unname(tab[term, "Pr(>|t|)"]),
      singular = lme4::isSingular(model, tol = 1e-4), status = "fitted",
      warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
      formula = paste(deparse(formula), collapse = " "))
  }

  scopes <- list(
    all_matched = d,
    excluding_6634 = d[id != "ADNI_6634"],
    excluding_6480 = d[id != "ADNI_6480"],
    excluding_both = d[!id %in% c("ADNI_6634", "ADNI_6480")],
    jump15_people_removed_exploratory = d[!id %in% shock_ids]
  )
  results <- rbindlist(lapply(names(scopes), function(scope) {
    dat <- scopes[[scope]]
    cat("24 | Fitting", scope, ":", uniqueN(dat$id), "people /", nrow(dat), "scans\n")
    rbindlist(lapply(c("random_intercept", "random_slope"), function(random)
      rbindlist(list(
        fit_one(dat, "fs_pct", scope, "FreeSurfer", random, covariates),
        fit_one(dat, "ashs_pct", scope, "ASHS", random, covariates),
        fit_one(dat, "method_diff_pct", scope, "FS minus ASHS", random, covariates)
      ))))
  }))
  save_tsv(results, "Table3_paired_method_models_and_influence.tsv")

  # This is a paired method sensitivity, so processing_version is removed in all arms.
  covars_no_version <- setdiff(covariates, "processing_version")
  no_version <- rbindlist(lapply(c("random_intercept", "random_slope"), function(random)
    rbindlist(list(
      fit_one(d, "fs_pct", "all_matched_no_version", "FreeSurfer", random, covars_no_version),
      fit_one(d, "ashs_pct", "all_matched_no_version", "ASHS", random, covars_no_version),
      fit_one(d, "method_diff_pct", "all_matched_no_version", "FS minus ASHS", random, covars_no_version)
    ))))
  save_tsv(no_version, "Table4_processing_version_sensitivity.tsv")


  # Participant bootstrap of the paired contrast retains all visits per sampled person.
  set.seed(20261112)
  people <- unique(d$id)
  boot_n <- 400L
  boot_formula <- make_formula("method_diff_pct", "boot_id", "random_intercept", covariates)
  bootstrap <- rbindlist(lapply(seq_len(boot_n), function(b) {
    if (b %% 50L == 0L) cat("24 | Participant bootstrap", b, "/", boot_n, "\n")
    drawn <- sample(people, length(people), replace = TRUE)
    z <- rbindlist(lapply(seq_along(drawn), function(k) {
      part <- copy(d[id == drawn[k]])
      part[, boot_id := paste0("rep", k)]
      part
    }))
    bw <- character()
    fit <- tryCatch(withCallingHandlers(
      lme4::lmer(boot_formula, data = z,
        REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa",
                            optCtrl = list(maxfun = 100000))),
      warning = function(w) {
        bw <<- c(bw, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    estimate <- if (inherits(fit, "error")) NA_real_ else
      unname(lme4::fixef(fit)["time_years:tau_z"])
    data.table(replicate = b, estimate,
               singular = if (inherits(fit, "error")) NA else
                 lme4::isSingular(fit, tol = 1e-4),
               status = if (inherits(fit, "error")) conditionMessage(fit) else "fitted",
               warnings = if (length(bw)) paste(unique(bw), collapse = " | ") else NA_character_)
  }))
  save_tsv(bootstrap, "Participant_bootstrap_paired_difference.tsv", "source_data")
  valid <- bootstrap[is.finite(estimate), estimate]
  if (length(valid) < 300L) stop("Too few successful bootstrap fits; inspect saved replicates")
  ci <- unname(quantile(valid, c(0.025, 0.975), names = FALSE))
  boot_summary <- data.table(method = "FS minus ASHS", statistic = "tau x time, %/year per tau SD",
    participants = length(people), scans = nrow(d), requested = boot_n,
    successful = length(valid), percentile_low = ci[1L], percentile_high = ci[2L],
    singular_fits = bootstrap[is.finite(estimate), sum(singular)],
    warning_fits = bootstrap[is.finite(estimate), sum(!is.na(warnings))],
    wald_estimate = results[scope == "all_matched" & method == "FS minus ASHS" &
                              random == "random_intercept", estimate],
    wald_p = results[scope == "all_matched" & method == "FS minus ASHS" &
                       random == "random_intercept", p])
  save_tsv(boot_summary, "Table5_participant_bootstrap_paired_difference.tsv")

  capture.output(sessionInfo(), file = file.path(out, "audit", "sessionInfo.txt"))
  print(results[, .(scope, method, random, estimate, ci_low, ci_high,
                    p, n_people, n_scans, singular, status)])
  print(boot_summary)
  message("Stage 24 complete: ", out)
})
