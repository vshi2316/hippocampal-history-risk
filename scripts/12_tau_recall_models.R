local({
# Fit the reported temporal-versus-posterior tau contrast and cognitive change models.
suppressPackageStartupMessages({library(data.table); library(sandwich); library(lmtest)})
input <- file.path(.study$output, "tau_recall_inputs/run_frozen")
out <- file.path(.study$output, "tau_recall_models/run_frozen")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_cohort <- function(cohort) {
  d <- fread(file.path(input, paste0(cohort, "_03_tau_to_future_memory.tsv")),
             na.strings = c("NA", ""))
  d[, cohort := cohort]
  d
}
cognition <- rbindlist(lapply(c("ADNI", "HABS_HD"), read_cohort), fill = TRUE)
model <- function(x, cohort, outcome, scope, adjustment) {
  d <- copy(x)
  if (scope %in% c("amyloid_positive", "amyloid_positive_same_scanner")) {
    d <- d[amyloid_pos == 1]
  }
  if (scope == "amyloid_negative") d <- d[amyloid_pos == 0]
  needed <- c(outcome, "tau_contrast", "tau_mean", "amyloid_cl", "age", "sex",
              adjustment)
  d <- d[complete.cases(d[, ..needed])]
  if (scope == "amyloid_interaction" && uniqueN(d$amyloid_pos) < 2L) {
    return(data.table(cohort, outcome, scope, n = nrow(d),
                      term = "one_amyloid_group", estimate = NA_real_,
                      se = NA_real_, p = NA_real_))
  }
  if (nrow(d) < 30L) return(data.table(cohort, outcome, scope, n = nrow(d),
                                       term = "sample_below_30", estimate = NA_real_,
                                       se = NA_real_, p = NA_real_))
  d[, tau_contrast_z := as.numeric(scale(tau_contrast))]
  d[, tau_mean_z := as.numeric(scale(tau_mean))]
  d[, amyloid_cl_z := as.numeric(scale(amyloid_cl))]
  contrast_term <- if (scope == "amyloid_interaction")
    "tau_contrast_z * factor(amyloid_pos)" else "tau_contrast_z"
  f <- reformulate(c(contrast_term, "tau_mean_z", "amyloid_cl_z",
                     "age", "factor(sex)", adjustment), response = outcome)
  fit <- lm(f, data = d)
  tab <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  data.table(cohort, outcome, scope, n = nobs(fit), term = rownames(tab),
             estimate = tab[, 1], se = tab[, 2], lower95 = tab[, 1] - 1.96 * tab[, 2],
             upper95 = tab[, 1] + 1.96 * tab[, 2], p = tab[, 4])
}

results <- rbindlist(lapply(c("ADNI", "HABS_HD"), function(cohort_name) {
  model(cognition[cohort == cohort_name], cohort_name, "memory_change_year",
        "amyloid_positive", "baseline_memory")
}), fill = TRUE)
fwrite(results, file.path(out, "Future_memory_model.tsv"), sep = "\t", na = "NA")
})
