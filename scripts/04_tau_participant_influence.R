# ADNI participant influence in the raw-unit tau models.
local({
suppressPackageStartupMessages(library(data.table))
root <- .study$output
out <- file.path(root,"tau_influence","run_frozen")
for(f in c("tables","audit")) dir.create(file.path(out,f),recursive=TRUE,showWarnings=FALSE)
save_tsv <- function(x,name,folder="tables") fwrite(as.data.table(x),file.path(out,folder,name),sep="\t",na="NA")
files <- c(model=file.path(root,"tau_longitudinal_inputs/run_frozen/models/Longitudinal_models.rds"),
 scan16=file.path(root,"tau_hippocampal_models/run_frozen/analysis_data/Frozen_scan_audit.tsv"),
 frozen16=file.path(root,"tau_hippocampal_models/run_frozen/tables/Table2_tau_slopes_standardized_and_raw.tsv"))
  frozen <- readRDS(files[["model"]])$MRI
  adni <- copy(as.data.table(frozen[["ADNI__Apositive_primary__cortex__tau_time"]]$data))
  adni[, `:=`(id = as.character(id), image_id = as.character(image_id))]
  if (uniqueN(adni$id) != 94L || nrow(adni) != 277L ||
      adni[id == "ADNI_6480", .N] != 4L)
    stop("The frozen Stage 16 ADNI model population changed")

  covariates16 <- c("sex_female", "education_c", "apoe4_carrier",
                    "diagnosis", "scanner", "processing_version")
  covariates16 <- covariates16[vapply(covariates16, function(v)
    v %in% names(adni) && uniqueN(adni[[v]][!is.na(adni[[v]])]) > 1L,
    logical(1))]
  fit16 <- function(x, scenario, outcome) {
    d <- copy(x)
    needed <- unique(c(outcome, "time_years", "tau_z", "amyloid_z",
                       "age_c", "id", covariates16))
    if (!all(needed %in% names(d)) || !all(complete.cases(d[, ..needed])))
      stop("Incomplete frozen model variables in ", scenario, " / ", outcome)
    d <- d[complete.cases(d[, ..covariates16])]
    for (v in covariates16) if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d, j = v, value = droplevels(factor(d[[v]])))
    formula <- as.formula(paste(outcome, "~",
      paste(c("time_years*tau_z", "time_years*amyloid_z", "time_years*age_c",
              covariates16, "time_years", "(1 + time_years | id)"), collapse = " + ")))
    warnings <- character()
    fit <- withCallingHandlers(lmerTest::lmer(formula, data = d, REML = FALSE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 100000))),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    tab <- summary(fit)$coefficients
    term <- "time_years:tau_z"
    beta <- unname(tab[term, "Estimate"])
    se <- unname(tab[term, "Std. Error"])
    p <- if ("Pr(>|t|)" %in% colnames(tab)) unname(tab[term, "Pr(>|t|)"]) else
      2 * pnorm(-abs(beta / se))
    scale <- if (outcome == "log_hip_ratio") 100 else 1
    data.table(scenario, metric = outcome, participants = uniqueN(d$id),
               scans = nrow(d), estimate = beta * scale,
               ci_low = (beta - 1.96 * se) * scale,
               ci_high = (beta + 1.96 * se) * scale, p,
               singular = lme4::isSingular(fit, tol = 1e-4),
               warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
               formula = paste(deparse(formula), collapse = " "))
  }
  datasets16 <- list(original = adni,
                     remove_ADNI_6480 = adni[id != "ADNI_6480"])
  outcomes <- c("cortex_mm","log_hip_ratio")
  effects16 <- rbindlist(lapply(names(datasets16), function(s)
    rbindlist(lapply(outcomes, function(y) fit16(datasets16[[s]], s, y)))))
  original16 <- fread(files[["frozen16"]], encoding = "UTF-8")
  original16 <- original16[scope == "ADNI_primary",
    .(metric, frozen_estimate = estimate, frozen_ci_low = ci_low,
      frozen_ci_high = ci_high, frozen_p = p, frozen_participants = participants,
      frozen_scans = scans)]
  reproduction16 <- merge(effects16[scenario == "original"], original16, by = "metric")
  save_tsv(reproduction16, "Stage16_exact_reproduction.tsv", "audit")
  if (nrow(reproduction16) != 2L ||
      any(abs(reproduction16$estimate - reproduction16$frozen_estimate) > 1e-8) ||
      any(abs(reproduction16$ci_low - reproduction16$frozen_ci_low) > 1e-8) ||
      any(reproduction16$participants != reproduction16$frozen_participants) ||
      any(reproduction16$scans != reproduction16$frozen_scans))
    stop("Stage 16 original model did not reproduce; no influence result can be interpreted")
  reference16 <- effects16[scenario == "original", .(metric, original_estimate = estimate)]
  effects16 <- merge(effects16, reference16, by = "metric")
  effects16[, change_from_original := estimate - original_estimate]
  setcolorder(effects16, c("scenario", "metric", "participants", "scans",
                           "estimate", "ci_low", "ci_high", "p", "original_estimate",
                           "change_from_original", "singular", "warnings", "formula"))
  save_tsv(effects16, "Table2_Stage16_tau_by_time_influence.tsv")

})
