local({
  needed <- c("data.table", "survival")
  absent <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install R packages: ", paste(absent, collapse = ", "))
  library(data.table)
  library(survival)

  root <- .study$output
  files <- c(
    risk = file.path(.study$output,"45_OASIS3_frozen_event_link_results/01_frozen_event_risk_set.csv"),
    parent = file.path(.study$output,"44_OASIS3_full_clinical_link_results/12_two_MRI_clinical_landmarks.csv"),
    muse = file.path(.study$output,"44_OASIS3_full_clinical_link_results/06_selected_MUSE_scans.csv"),
    metadata = file.path(root, "oasis_tau_mri_gate_27/run_frozen/audit/input_csv_cache/OASIS3_MR_scan_summary.csv"),
    reference = file.path(root, "oasis_event_38/run_frozen/tables/Table4_primary_increment_and_frozen_reproduction.tsv"),
    bootstrap = file.path(root, "oasis_bootstrap_39/run_frozen/tables/Table3_optimism_corrected_C_index.tsv")
  )
  missing_files <- files[!file.exists(files)]
  if (length(missing_files)) stop("Missing input: ", paste(missing_files, collapse = ", "))
  out <- file.path(root, "oasis_event_40", paste0("run_", "frozen"))
  for (subdir in c("tables", "source_data", "audit"))
    dir.create(file.path(out, subdir), recursive = TRUE, showWarnings = FALSE)
  write_tsv <- function(x, name, subdir = "tables")
    fwrite(as.data.table(x), file.path(out, subdir, name), sep = "\t", na = "NA")
  write_tsv(data.table(input = names(files), path = unname(files),
                       md5 = unname(tools::md5sum(files))), "Input_manifest.tsv", "audit")

  message("Reading the frozen OASIS clinical risk set")
  risk <- fread(files[["risk"]])
  risk[, event := as.logical(event)]
  cn <- copy(risk[baseline_stage == "CN"])
  if (nrow(cn) != 312L || sum(cn$event) != 52L || anyDuplicated(cn$OASISID))
    stop("The frozen CN risk set no longer has 312 people and 52 events")
  fields <- c("OASISID", "followup_days", "event", "age_MRI2", "sex", "education",
              "baseline_CDRSB", "hippo2_pct_ICV", "hippo_loss_pct_year")
  analysis <- as.data.frame(cn[complete.cases(cn[, ..fields])])
  analysis$sex <- factor(analysis$sex)
  analysis$followup_years <- analysis$followup_days / 365.25
  if (nrow(analysis) != 312L || sum(analysis$event) != 52L)
    stop("The original Cox model has a different complete-case sample")
  formula_base <- Surv(followup_years, event) ~ age_MRI2 + sex + education +
    baseline_CDRSB + hippo2_pct_ICV
  formula_full <- update(formula_base, . ~ . + hippo_loss_pct_year)
  fit_base <- coxph(formula_base, data = analysis, ties = "efron")
  fit_full <- coxph(formula_full, data = analysis, ties = "efron")
  extract <- function(fit, label, base_fit = NULL) {
    beta <- unname(coef(fit)["hippo_loss_pct_year"])
    se <- sqrt(vcov(fit)["hippo_loss_pct_year", "hippo_loss_pct_year"])
    data.table(model = label, participants = as.integer(fit$n),
               events = as.integer(fit$nevent), log_HR = beta,
               HR = exp(beta), CI_low = exp(beta - 1.96 * se),
               CI_high = exp(beta + 1.96 * se),
               Wald_p = 2 * pnorm(-abs(beta / se)),
               LR_chisq = if (is.null(base_fit)) NA_real_ else
                 2 * (fit$loglik[2] - base_fit$loglik[2]),
               apparent_delta_C = if (is.null(base_fit)) NA_real_ else
                 unname(summary(fit)$concordance[1] - summary(base_fit)$concordance[1]),
               status = "FITTED")
  }
  primary <- extract(fit_full, "Frozen primary CN model", fit_base)
  ref <- fread(files[["reference"]])
  if (nrow(ref) != 1L || abs(primary$HR - ref$HR_per_1pct_year) > 1e-6 ||
      abs(primary$LR_chisq - ref$LR_chisq) > 1e-6 ||
      abs(primary$apparent_delta_C - ref$apparent_delta_c_index) > 1e-6)
    stop("The Step 38 primary model was not reproduced")
  write_tsv(primary, "Table1_frozen_primary_reproduction.tsv")
  write_tsv(fread(files[["bootstrap"]]), "Table2_existing_bootstrap_discrimination.tsv")

  message("Calculating one-person-out influence on the fixed Cox model")
  loo <- rbindlist(lapply(seq_len(nrow(analysis)), function(i) {
    fit <- tryCatch(coxph(formula_full, data = analysis[-i, , drop = FALSE], ties = "efron"),
                    error = function(e) e)
    if (inherits(fit, "error"))
      return(data.table(OASISID = analysis$OASISID[i],
                        log_HR = NA_real_, HR = NA_real_, status = conditionMessage(fit)))
    beta <- unname(coef(fit)["hippo_loss_pct_year"])
    data.table(OASISID = analysis$OASISID[i], log_HR = beta,
               HR = exp(beta), status = "FITTED")
  }))
  loo[, log_HR_shift := log_HR - primary$log_HR]
  loo[, abs_log_HR_shift := abs(log_HR_shift)]
  setorderv(loo, "abs_log_HR_shift", order = -1L, na.last = TRUE)
  loo[, influence_rank := seq_len(.N)]
  write_tsv(loo, "Table3_all_312_leave_one_out_effects.tsv", "source_data")
  write_tsv(loo[1:10], "Table4_ten_largest_person_level_influences.tsv")

  message("Linking MUSE scans to the MRI acquisition metadata")
  muse <- fread(files[["muse"]], select = c("OASISID", "mri_day", "run", "Code"))
  meta <- fread(files[["metadata"]], select = c("subject_id", "label", "filename",
                 "acccession", "MagneticFieldStrength", "ManufacturersModelName",
                 "DeviceSerialNumber", "SeriesDescription"))
  meta <- meta[grepl("^OAS[0-9]{5}_MR_d[0-9]+$", label) &
                 grepl("_T1w\\.json$", filename)]
  meta[, `:=`(OASISID = subject_id,
              mri_day = as.integer(sub("^.*_MR_d([0-9]+)$", "\\1", label)),
              run = 0L)]
  has_run <- grepl("_run-[0-9]+_T1w\\.json$", meta$filename)
  meta[has_run, run := as.integer(sub("^.*_run-([0-9]+)_T1w\\.json$", "\\1", filename))]
  setnames(meta, c("acccession", "MagneticFieldStrength", "ManufacturersModelName",
                   "DeviceSerialNumber", "SeriesDescription"),
           c("accession", "field_strength", "scanner", "serial", "series"))
  meta[, serial := as.character(serial)]
  meta[!nzchar(trimws(serial)), serial := NA_character_]
  meta <- unique(meta[, .(OASISID, mri_day, run, filename, accession,
                          field_strength, scanner, serial, series)])
  meta_key <- meta[, .(candidate_files = .N,
                       filename = if (.N == 1L) filename[1] else NA_character_,
                       accession = if (.N == 1L) accession[1] else NA_character_,
                       field_strength = if (.N == 1L) as.character(field_strength[1]) else NA_character_,
                       scanner = if (.N == 1L) scanner[1] else NA_character_,
                       serial = if (.N == 1L) serial[1] else NA_character_,
                       series = if (.N == 1L) series[1] else NA_character_),
                   by = .(OASISID, mri_day, run)]
  if (anyDuplicated(muse$Code)) stop("MUSE scan Code is not unique")
  linked <- merge(muse, meta_key, by = c("OASISID", "mri_day", "run"), all.x = TRUE, sort = FALSE)
  linked[, metadata_status := fifelse(is.na(candidate_files), "NOT_FOUND",
                                      fifelse(candidate_files == 1L, "EXACT", "AMBIGUOUS"))]
  ambiguous_keys <- unique(linked[metadata_status == "AMBIGUOUS",
                                  .(OASISID, mri_day, run)])
  ambiguous_candidates <- merge(ambiguous_keys, meta,
                                by = c("OASISID", "mri_day", "run"), all.x = TRUE)
  write_tsv(ambiguous_candidates, "Table6b_ambiguous_T1_metadata_candidates.tsv", "source_data")
  m1 <- linked[, .(OASISID, MRI1_Code = Code, MRI1_metadata_status = metadata_status,
                   MRI1_scanner = scanner, MRI1_serial = serial, MRI1_series = series,
                   MRI1_field = field_strength,
                   MRI1_filename = filename, MRI1_accession = accession)]
  m2 <- linked[, .(OASISID, MRI2_Code = Code, MRI2_metadata_status = metadata_status,
                   MRI2_scanner = scanner, MRI2_serial = serial, MRI2_series = series,
                   MRI2_field = field_strength,
                   MRI2_filename = filename, MRI2_accession = accession)]
  joined <- merge(cn, m1, by = c("OASISID", "MRI1_Code"), all.x = TRUE, sort = FALSE)
  joined <- merge(joined, m2, by = c("OASISID", "MRI2_Code"), all.x = TRUE, sort = FALSE)
  if (nrow(joined) != 312L) stop("MUSE to MRI metadata join changed the risk set")
  for (phase in c("MRI1", "MRI2"))
    set(joined, which(is.na(joined[[paste0(phase, "_metadata_status")]])),
        paste0(phase, "_metadata_status"), "NOT_FOUND")
  joined[, both_exact := MRI1_metadata_status == "EXACT" & MRI2_metadata_status == "EXACT" &
           !is.na(MRI1_scanner) & !is.na(MRI2_scanner)]
  joined[, serial_known := both_exact & !is.na(MRI1_serial) & !is.na(MRI2_serial)]
  joined[, same_scanner := serial_known & MRI1_serial == MRI2_serial &
           MRI1_scanner == MRI2_scanner]
  joined[, same_scanner_series := same_scanner & !is.na(MRI1_series) &
           !is.na(MRI2_series) & MRI1_series == MRI2_series &
           !is.na(MRI1_field) & !is.na(MRI2_field) & MRI1_field == MRI2_field]
  joined[, scanner_changed := as.integer(serial_known & !same_scanner)]
  metadata_flow <- joined[, .(participants = .N, events = as.integer(sum(event))),
                          by = .(MRI1_metadata_status, MRI2_metadata_status, both_exact, serial_known,
                                 same_scanner, same_scanner_series)]
  write_tsv(metadata_flow, "Table5_scan_metadata_and_protocol_flow.tsv")
  write_tsv(joined[, .(OASISID, MRI1_Code, MRI2_Code, MRI1_day, MRI2_day,
                        MRI1_metadata_status, MRI2_metadata_status,
                        MRI1_scanner, MRI2_scanner, MRI1_serial, MRI2_serial,
                        MRI1_series, MRI2_series,
                        MRI1_field, MRI2_field, MRI1_accession, MRI2_accession,
                        MRI1_filename, MRI2_filename, both_exact, serial_known, same_scanner,
                        same_scanner_series, event)],
            "Table6_all_312_exact_scan_crosswalk.tsv", "source_data")

  sensitivity_fit <- function(dt, label, switch_adjusted = FALSE) {
    n <- as.integer(nrow(dt))
    ev <- as.integer(sum(dt$event))
    blank <- data.table(model = label, participants = n, events = ev,
                        log_HR = NA_real_, HR = NA_real_, CI_low = NA_real_,
                        CI_high = NA_real_, Wald_p = NA_real_, LR_chisq = NA_real_,
                        apparent_delta_C = NA_real_, status = "INSUFFICIENT_EVENTS")
    if (n < 30L || ev < 10L) return(blank)
    d <- as.data.frame(dt)
    d$sex <- factor(d$sex, levels = levels(analysis$sex))
    d$followup_years <- d$followup_days / 365.25
    f0 <- formula_base
    f1 <- formula_full
    if (switch_adjusted) {
      if (length(unique(d$scanner_changed)) < 2L) {
        blank$status <- "NO_SCANNER_TRANSITION_VARIATION"
        return(blank)
      }
      f0 <- update(f0, . ~ . + scanner_changed)
      f1 <- update(f1, . ~ . + scanner_changed)
    }
    fits <- tryCatch(list(base = coxph(f0, data = d, ties = "efron"),
                          full = coxph(f1, data = d, ties = "efron")),
                     error = function(e) e)
    if (inherits(fits, "error")) {
      blank$status <- paste0("MODEL_ERROR: ", conditionMessage(fits))
      return(blank)
    }
    extract(fits$full, label, fits$base)
  }
  s1 <- sensitivity_fit(joined[both_exact == TRUE], "Both scans: exact MRI metadata")
  s2 <- sensitivity_fit(joined[same_scanner == TRUE], "Same scanner at MRI1 and MRI2")
  s3 <- sensitivity_fit(joined[same_scanner_series == TRUE], "Same scanner and T1 series")
  s4 <- sensitivity_fit(joined[serial_known == TRUE], "Scanner change adjusted", TRUE)
  sensitivities <- rbindlist(list(primary, s1, s2, s3, s4))
  write_tsv(sensitivities, "Table7_frozen_model_protocol_sensitivities.tsv")

  message("Describing inclusion from the same two-MRI parent sample")
  parent <- fread(files[["parent"]])
  parent <- parent[baseline_stage == "CN"]
  parent[, in_event_risk_set := OASISID %in% cn$OASISID]
  if (!setequal(parent[in_event_risk_set == TRUE]$OASISID, cn$OASISID))
    stop("The CN event sample is not contained in the two-MRI parent")
  vars <- c("age_MRI2", "education", "baseline_CDRSB", "hippo2_pct_ICV",
            "hippo_loss_pct_year", "MRI_span_days", "D1_gap_days")
  selection <- rbindlist(lapply(vars, function(v) {
    parent[, .(variable = v, participants = .N, observed = as.integer(sum(is.finite(get(v)))),
               median = as.numeric(median(get(v), na.rm = TRUE)), q1 = as.numeric(quantile(get(v), .25, na.rm = TRUE)),
               q3 = as.numeric(quantile(get(v), .75, na.rm = TRUE))),
           by = .(in_event_risk_set)]
  }))
  write_tsv(selection, "Table8_CN_two_MRI_selection_descriptives.tsv")
  write_tsv(parent[, .(participants = .N, events_recorded = as.integer(sum(event, na.rm = TRUE)),
                        sex_1 = as.integer(sum(sex == 1, na.rm = TRUE)),
                        sex_2 = as.integer(sum(sex == 2, na.rm = TRUE))),
                   by = in_event_risk_set], "Table9_CN_two_MRI_selection_counts.tsv")


  decision <- data.table(
    item = c("Frozen participants", "First dementia records", "Both scans with exact metadata",
             "Both scanner serials known", "Same physical scanner", "Same scanner and T1 series", "Unresolved scan metadata",
             "Top-10 influence cases", "Person-level models fitted"),
    value = c(nrow(cn), sum(cn$event), sum(joined$both_exact), sum(joined$serial_known),
              sum(joined$same_scanner), sum(joined$same_scanner_series),
              sum(joined$MRI1_metadata_status != "EXACT" | joined$MRI2_metadata_status != "EXACT" |
                    is.na(joined$MRI1_scanner) | is.na(joined$MRI2_scanner) |
                    is.na(joined$MRI1_serial) | is.na(joined$MRI2_serial) |
                    is.na(joined$MRI1_series) | is.na(joined$MRI2_series)),
              10L, sum(loo$status == "FITTED")))
  write_tsv(decision, "Table11_audit_summary.tsv")
  writeLines(capture.output(sessionInfo()), file.path(out, "audit", "R_sessionInfo.txt"))
  writeLines("STATUS=STEP40_COMPLETE", file.path(out, "audit", "40_COMPLETED.txt"))
  message("Step 40 complete: ", out)
  print(primary)
  print(sensitivities)
  print(decision)
})
