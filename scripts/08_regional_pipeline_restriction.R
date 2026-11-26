# Amyloid-positive cortical-band contrasts in full and scanner-stable samples.

local({
  packages <- c("data.table", "lme4", "lmerTest")
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  suppressPackageStartupMessages(library(data.table))

  root <- .study$output
  stage34_root <- file.path(root, "roi68_spatial_34")
  runs <- sort(list.files(stage34_root, pattern = "^run_frozen$",
                          full.names = TRUE))
  runs <- runs[file.exists(file.path(runs, "source_data", "Paired_band_scan_data.tsv"))]
  if (!length(runs)) stop("Run scripts/07_regional_tau_models.R first")
  stage34 <- tail(runs, 1L)
  scan_file <- file.path(stage34, "source_data", "Paired_band_scan_data.tsv")
  primary_file <- file.path(stage34, "tables", "Table6_direct_paired_band_differences.tsv")
  out <- file.path(root, "cortical_band_35",
                   paste0("run_", "frozen"))
  for (subdir in c("tables", "source_data", "audit"))
    dir.create(file.path(out, subdir), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, subdir = "tables")
    fwrite(as.data.table(x), file.path(out, subdir, name), sep = "\t",
           quote = TRUE, na = "NA", bom = TRUE)
  log_line <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                " | ", ..., "\n", sep = "")
  log_line("Reading corrected Stage 34: ", stage34)
  scan <- fread(scan_file, encoding = "UTF-8")
  previous <- fread(primary_file, encoding = "UTF-8")
  scan <- scan[group == "Apositive"]
  scan[, stable_pipeline := as.logical(stable_pipeline)]
  keys <- c("cohort", "id", "mri_date", "image_id")
  if (anyDuplicated(scan, by = keys)) stop("Repeated frozen MRI scan keys")
  counts <- scan[, .(people = uniqueN(id), scans = .N), by = cohort]
  if (nrow(counts) != 2L ||
      counts[cohort == "ADNI", people] != 94L ||
      counts[cohort == "ADNI", scans] != 277L ||
      counts[cohort == "HABS-HD", people] != 57L ||
      counts[cohort == "HABS-HD", scans] != 127L)
    stop("The Stage 34 band data do not reproduce the frozen Stage 16 scans")

  contrasts <- list(
    "Medial minus other temporal" = c("Medial temporal", "Other temporal"),
    "Medial minus non-temporal" = c("Medial temporal", "Non-temporal"),
    "Other temporal minus non-temporal" = c("Other temporal", "Non-temporal"))
  covariates <- c("sex_female", "education_c", "apoe4_carrier", "diagnosis",
                  "scanner", "processing_version")
  needed <- c("id", "time_years", "tau_z", "amyloid_z", "age10", "y", covariates)
  prepare <- function(cohort_name, contrast_name, stable) {
    pair <- contrasts[[contrast_name]]
    d <- copy(scan[scan$cohort == cohort_name])
    if (stable) d <- d[stable_pipeline == TRUE]
    d[, y := get(pair[1L]) - get(pair[2L])]
    d <- d[complete.cases(d[, ..needed])]
    eligible <- d[, .(n_scans = .N, baseline = any(time_years == 0)), by = id]
    d[id %in% eligible[baseline == TRUE & n_scans >= 2L, id]]
  }
  estimate_stable <- function(d, cohort_name, contrast_name) {
    row <- data.table(cohort = cohort_name, scope = "Stable", feature = contrast_name,
                      estimate = NA_real_, se = NA_real_, ci_low = NA_real_,
                      ci_high = NA_real_, p = NA_real_, n_people = uniqueN(d$id),
                      n_scans = nrow(d), singular = NA, status = "insufficient_data")
    if (row$n_people < 20L || row$n_scans < 40L) return(row)
    covs <- covariates[vapply(covariates, function(v)
      uniqueN(d[[v]]) > 1L, logical(1))]
    for (v in covs) if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d, j = v, value = droplevels(factor(d[[v]])))
    formula <- as.formula(paste("y ~", paste(c(
      "time_years * tau_z", "time_years * age10", "time_years * amyloid_z",
      covs, "(1 | id)"), collapse = " + ")))
    fit <- tryCatch(suppressWarnings(lmerTest::lmer(
      formula, data = d, REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 100000)))), error = identity)
    if (inherits(fit, "error")) {
      row[, status := paste0("model_error: ", conditionMessage(fit))]
      return(row)
    }
    cf <- summary(fit)$coefficients
    term <- "time_years:tau_z"
    if (!term %in% rownames(cf)) stop("Missing Tau-by-time coefficient")
    beta <- unname(cf[term, "Estimate"])
    beta_se <- unname(cf[term, "Std. Error"])
    beta_p <- unname(cf[term, "Pr(>|t|)"])
    if (!all(is.finite(c(beta, beta_se, beta_p))))
      stop("Non-finite fitted Tau estimate: ", cohort_name, ", ", contrast_name)
    row[, `:=`(estimate = beta, se = beta_se, ci_low = beta - 1.96 * beta_se,
               ci_high = beta + 1.96 * beta_se, p = beta_p,
               singular = lme4::isSingular(fit, tol = 1e-4), status = "fitted")]
    row
  }

  flow_rows <- list()
  result_rows <- list()
  index <- 0L
  for (co in c("ADNI", "HABS-HD")) for (name in names(contrasts)) {
    full <- prepare(co, name, stable = FALSE)
    stable <- prepare(co, name, stable = TRUE)
    old <- previous[cohort == co & group == "Apositive" &
                      feature == name & term == "time_years:tau_z"]
    if (nrow(old) != 1L || old$n_people != uniqueN(full$id) ||
        old$n_scans != nrow(full))
      stop("Stage 35 full-sample denominator differs from Stage 34: ", co, ", ", name)
    for (scope_name in c("Full", "Stable")) {
      d <- if (scope_name == "Full") full else stable
      index <- index + 1L
      flow_rows[[index]] <- data.table(
        cohort = co, feature = name, scope = scope_name,
        people = uniqueN(d$id), scans = nrow(d),
        baseline_scans = sum(d$time_years == 0),
        followup_scans = sum(d$time_years > 0))
    }
    result_rows[[length(result_rows) + 1L]] <- old[, .(
      cohort, scope = "Full", feature, estimate, se, ci_low, ci_high,
      p, n_people, n_scans, singular, status)]
    result_rows[[length(result_rows) + 1L]] <- estimate_stable(stable, co, name)
    log_line(co, " | ", name, " | stable ", uniqueN(stable$id),
             " people / ", nrow(stable), " scans")
  }
  flow <- rbindlist(flow_rows)
  results <- rbindlist(result_rows)
  results[, p_holm_three := p.adjust(p, method = "holm", n = 3L),
          by = .(cohort, scope)]
  setorder(flow, cohort, feature, scope)
  setorder(results, cohort, feature, scope)
  save_tsv(flow, "Table1_band_contrast_sample_flow.tsv")
  save_tsv(results, "Table2_full_and_stable_tau_effects.tsv")
  save_tsv(results, "Figure35_source_data.tsv", "source_data")

  save_tsv(data.table(input = c("Stage34_band_scans", "Stage34_primary_estimates"),
                      path = c(scan_file, primary_file),
                      md5 = unname(tools::md5sum(c(scan_file, primary_file)))),
           "Input_manifest.tsv", "audit")
  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "sessionInfo.txt"))
  log_line("Analysis complete: ", out)
  print(results[, .(cohort, feature, scope, n_people, n_scans,
                    estimate, ci_low, ci_high, p, p_holm_three, status)])
})
