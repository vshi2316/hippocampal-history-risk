local({
  if (!requireNamespace("data.table", quietly = TRUE) ||
      !requireNamespace("lme4", quietly = TRUE) ||
      !requireNamespace("survival", quietly = TRUE))
    stop("Install data.table, lme4, and survival before running this script.")
  library(data.table)

  root <- .study$output
  paths <- c(
    cognitive_visits = file.path(root, "cognitive_trajectory_63/run_frozen/tables/Table2_person_visit_cognitive_trajectories.tsv"),
    clinical_course = file.path(root, "cn_progression_62/run_frozen/tables/Table1_person_level_clinical_course.tsv")
  )
  if (any(!file.exists(paths)))
    stop("Missing input: ", paste(paths[!file.exists(paths)], collapse = "; "))
  out <- file.path(root, "prediagnostic_64",
                   paste0("run_", "frozen"))
  for (folder in c("tables", "source_data", "audit"))
    dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
  save_tsv <- function(x, name, folder = "tables")
    fwrite(as.data.table(x), file.path(out, folder, name), sep = "\t", na = "NA")
  save_tsv(data.table(input = names(paths), file = unname(paths),
                      md5 = unname(tools::md5sum(paths))),
           "Input_manifest.tsv", "audit")

  clinical <- fread(paths[["clinical_course"]])
  if (nrow(clinical[cohort == "Mayo"]) != 925L ||
      nrow(clinical[cohort == "OASIS-3"]) != 312L ||
      sum(clinical[cohort == "Mayo", dementia_event]) != 21L ||
      sum(clinical[cohort == "OASIS-3", dementia_event]) != 52L)
    stop("The frozen clinical cohort changed.")
  visits <- fread(paths[["cognitive_visits"]])
  keys <- c("cohort", "id")
  clinical_fields <- c(keys, "prior_MCI_record", "impairment_event",
                       "impairment_time_years", "dementia_event")
  x <- merge(visits, clinical[, ..clinical_fields], by = keys)
  x <- x[!(cohort == "OASIS-3" & prior_MCI_record == 1L)]
  x <- x[(impairment_event == 1L &
            time_years < impairment_time_years) |
           (impairment_event == 0L &
              time_years <= impairment_time_years + 1e-6)]
  save_tsv(x, "Table2_person_visit_before_first_impairment.tsv")
  coverage <- x[, .(visits = .N, last_year = max(time_years)),
                by = .(cohort, role, id)][,
                  .(participants = .N,
                    participants_2plus = sum(visits >= 2L),
                    participants_3plus = sum(visits >= 3L),
                    observations = sum(visits),
                    median_last_year = median(last_year)),
                  by = .(cohort, role)]
  save_tsv(coverage, "Table1_prediagnostic_cognition_coverage.tsv")

  selected <- x[role %in% c("Primary_Mayo_global_z",
                            "OASIS_administered_MMSE")]
  counts <- selected[, .N, by = .(cohort, role, id)][N >= 2L]
  setkeyv(selected, c("cohort", "role", "id"))
  landmark <- selected[counts[, .(cohort, role, id)],
                       on = .(cohort, role, id), nomatch = 0L]
  setorder(landmark, cohort, role, id, time_years)
  landmark <- landmark[, .SD[2L], by = .(cohort, role, id)]
  landmark[, `:=`(risk_years = impairment_time_years - time_years,
                  age_landmark = age_MRI2 + time_years)]
  landmark <- landmark[risk_years > 0]
  save_tsv(landmark, "Table4_second_cognitive_visit_landmark.tsv")

  fit_event <- function(d) {
    d <- as.data.frame(copy(d))
    d$sex <- factor(d$sex)
    base <- survival::coxph(
      survival::Surv(risk_years, impairment_event) ~ age_landmark +
        sex + education + hippo2_pct_ICV + baseline_score + score +
        time_years, data = d, ties = "efron")
    plus <- survival::coxph(
      survival::Surv(risk_years, impairment_event) ~ age_landmark +
        sex + education + hippo2_pct_ICV + baseline_score + score +
        time_years + hippo_loss_pct_year,
      data = d, ties = "efron")
    term <- "hippo_loss_pct_year"
    beta <- unname(coef(plus)[term])
    se <- sqrt(vcov(plus)[term, term])
    lr <- max(0, 2 * (as.numeric(logLik(plus)) -
                       as.numeric(logLik(base))))
    data.table(cohort = d$cohort[1L], participants = nrow(d),
               later_first_impairment_events = sum(d$impairment_event),
               HR_per_1pct_year = exp(beta),
               CI_low = exp(beta - 1.96 * se),
               CI_high = exp(beta + 1.96 * se),
               Wald_p = 2 * pnorm(-abs(beta / se)),
               LR_p = pchisq(lr, 1, lower.tail = FALSE),
               apparent_C_base = unname(summary(base)$concordance[1L]),
               apparent_C_plus_loss = unname(summary(plus)$concordance[1L]))
  }
  event_models <- rbindlist(lapply(split(landmark, landmark$cohort),
                                  fit_event))
  save_tsv(event_models, "Table5_updated_landmark_first_impairment.tsv")

writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
})
