local({
  root <- .study$work
  atlas_dir <- file.path(root, "inputs", "ENIGMA_BigBrain")
  out_dir <- file.path(root, "outputs", "normative_input_audit_59")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  effect_file <- paste0(
    file.path(.study$output,"roi68_spatial_34/"),
    "run_frozen/tables/Table5_68_region_tau_and_age_effects.tsv"
  )
  gradient_file <- file.path(atlas_dir, "bb_gradient_aparc.csv")
  label_file <- file.path(atlas_dir, "allepi_case-controls_CortThick.csv")
  required <- c(effect_file, gradient_file, label_file)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Missing input file(s):\n", paste(missing, collapse = "\n"))
  }

  gradient <- scan(gradient_file, what = numeric(), sep = ",", quiet = TRUE)
  labels <- read.csv(label_file, stringsAsFactors = FALSE)$Structure
  if (length(gradient) != 68L || length(labels) != 68L ||
      anyDuplicated(labels) || anyNA(gradient)) {
    stop("ENIGMA gradient and ordered DK labels must each contain 68 valid entries")
  }
  reference <- data.frame(
    side = substr(labels, 1L, 1L),
    roi = sub("^[LR]_", "", labels),
    bigbrain_gradient = gradient,
    stringsAsFactors = FALSE
  )
  if (!setequal(reference$side, c("L", "R")) ||
      any(table(reference$side) != 34L)) {
    stop("ENIGMA label order is not 34 left plus 34 right DK regions")
  }

  effects <- read.delim(effect_file, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("cohort", "group", "feature_type", "random", "term", "status",
              "side", "roi", "estimate", "se", "n_people", "n_scans")
  if (!all(needed %in% names(effects))) {
    stop("Stage 34 table is missing: ", paste(setdiff(needed, names(effects)), collapse = ", "))
  }
  effects$side <- trimws(effects$side)
  effects$roi <- trimws(effects$roi)
  cortical <- effects[
    effects$cohort %in% c("ADNI", "HABS-HD") &
      effects$feature_type == "cortical_region" &
      effects$random == "intercept" &
      effects$status == "fitted",
  ]
  target <- cortical[
    cortical$group == "Apositive" & cortical$term == "time_years:tau_z",
  ]
  if (anyDuplicated(target[c("cohort", "side", "roi")])) {
    stop("The frozen Apositive Tau-by-time table has duplicated cohort/side/ROI keys")
  }
  if (!identical(sort(unique(target$cohort)), c("ADNI", "HABS-HD")) ||
      any(table(target$cohort) != 68L)) {
    stop("Expected exactly 68 fitted Apositive Tau-by-time regions in each cohort")
  }

  matched <- merge(target, reference, by = c("side", "roi"), all.x = TRUE,
                   sort = FALSE)
  if (nrow(matched) != 136L || anyNA(matched$bigbrain_gradient) ||
      any(!is.finite(matched$estimate)) || any(!is.finite(matched$se)) ||
      !setequal(paste(reference$side, reference$roi),
                paste(target$side[target$cohort == "ADNI"],
                      target$roi[target$cohort == "ADNI"]))) {
    stop("The BigBrain gradient does not match all 68 frozen DK regions")
  }

  write.table(reference, file.path(out_dir, "DK68_BigBrain_reference.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(matched, file.path(out_dir, "Apositive_Tau_effects_matched_to_BigBrain.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

})
