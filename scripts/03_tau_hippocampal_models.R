# Raw-unit Tau by time models, processing restrictions and observed changes.
local({
  required <- c("data.table", "lme4", "lmerTest")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop("Install R packages: ", paste(missing, collapse=", "))
  suppressPackageStartupMessages(library(data.table))

  cfg <- list(
    script=file.path(.study$root,"run_all.R"),
    stage11=file.path(.study$output,"tau_longitudinal_inputs/run_frozen"),
    output_root=file.path(.study$output,"tau_hippocampal_models"))
  model_path <- file.path(cfg$stage11,"models","Longitudinal_models.rds")
  for (p in c(model_path)) if (!file.exists(p)) stop("Missing input: ",p)
  out <- file.path(cfg$output_root,paste0("run_","frozen"))
  if (dir.exists(out)) stop("Output directory already exists: ",out)
  for (folder in c("tables","source_data","analysis_data","audit"))
    dir.create(file.path(out,folder),recursive=TRUE,showWarnings=FALSE)
  save_tsv <- function(x,name,folder="tables")
    fwrite(as.data.table(x),file.path(out,folder,name),sep="\t",quote=TRUE,na="NA",bom=TRUE)
  log_file <- file.path(out,"16_run_log.txt")
  log_line <- function(message) {
    line <- paste(format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"|",message)
    cat(line,"\n")
    cat(line,"\n",file=log_file,append=TRUE)
  }

  log_line("Reading frozen Stage 11 model data")
  frozen <- readRDS(model_path)$MRI
  key <- function(cohort,scope,tissue)
    paste(cohort,scope,tissue,"tau_time",sep="__")
  data_for <- function(cohort,scope) {
    x <- copy(as.data.table(frozen[[key(cohort,scope,"cortex")]]$data))
    x[, id := as.character(id)]
    x
  }
  adni <- data_for("ADNI","Apositive_primary")
  habs <- data_for("HABS-HD","Apositive_primary")
  stable <- data_for("HABS-HD","Apositive_stable_pipeline")
  stable_ids <- unique(stable$id)
  matched <- copy(habs[id %in% stable_ids])
  stopifnot(uniqueN(adni$id)==94L,uniqueN(habs$id)==57L,
            uniqueN(stable$id)==30L,setequal(stable_ids,unique(matched$id)),
            all(stable$stable_pipeline),
            !anyDuplicated(adni[,.(id,mri_date)]),
            !anyDuplicated(habs[,.(id,mri_date)]))
  scopes <- list(ADNI_primary=adni,HABS_primary=habs,
                 HABS_same30_all_scans=matched,HABS_same30_stable_scans=stable)
  flow <- rbindlist(lapply(names(scopes),function(s) {
    x <- scopes[[s]]
    data.table(scope=s,cohort=unique(x$cohort),participants=uniqueN(x$id),
      scans=nrow(x),baseline_scans=sum(x$time_years==0),
      future_scans=sum(x$time_years>0),
      future_scans_stable=sum(x$time_years>0 & x$stable_pipeline))
  }))
  stopifnot(flow[scope=="HABS_same30_stable_scans",future_scans]>0)
  save_tsv(flow,"Table1_frozen_sample_flow.tsv")

  # Every refit uses Stage 11's exact participant/scanner rows and adjustment set.
  fit_one <- function(x,scope,outcome) {
    d <- copy(x)
    needed <- c(outcome,"time_years","tau_z","amyloid_z","age_c","id")
    stopifnot(all(needed %in% names(d)),all(complete.cases(d[,..needed])))
    covariates <- c("sex_female","education_c","apoe4_carrier",
                    "diagnosis","scanner","processing_version")
    covariates <- covariates[vapply(covariates,function(v)
      v %in% names(d) && uniqueN(d[[v]][!is.na(d[[v]])])>1L,logical(1))]
    d <- d[complete.cases(d[,..covariates])]
    for (v in covariates) if (is.character(d[[v]]) || is.factor(d[[v]]))
      set(d,j=v,value=droplevels(factor(d[[v]])))
    random <- if (scope=="ADNI_primary") "(1 + time_years | id)" else "(1 | id)"
    rhs <- paste(c("time_years*tau_z","time_years*amyloid_z",
                   "time_years*age_c",covariates,"time_years",random),collapse=" + ")
    formula <- as.formula(paste(outcome,"~",rhs))
    warnings <- character()
    fit <- withCallingHandlers(lmerTest::lmer(formula,data=d,REML=FALSE,
      control=lme4::lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))),
      warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    tab <- summary(fit)$coefficients
    term <- "time_years:tau_z"
    if (!term %in% rownames(tab)) stop("Missing tau-by-time term: ",scope," / ",outcome)
    beta <- unname(tab[term,"Estimate"])
    se <- unname(tab[term,"Std. Error"])
    p <- if ("Pr(>|t|)" %in% colnames(tab)) unname(tab[term,"Pr(>|t|)"]) else
      2*pnorm(-abs(beta/se))
    unit_scale <- if (outcome=="log_hip_ratio") 100 else 1
    unit <- switch(outcome,C="injury SD/year per tau SD",
      H="injury SD/year per tau SD",cortex_mm="mm/year per tau SD",
      log_hip_ratio="100 x log(volume/ICV)/year per tau SD")
    data.table(scope=scope,cohort=unique(d$cohort),metric=outcome,
      estimate=beta*unit_scale,ci_low=(beta-1.96*se)*unit_scale,
      ci_high=(beta+1.96*se)*unit_scale,p=p,unit=unit,
      participants=uniqueN(d$id),scans=nrow(d),
      singular=lme4::isSingular(fit,tol=1e-4),
      warnings=if(length(warnings))paste(unique(warnings),collapse=" | ") else NA_character_,
      formula=paste(deparse(formula),collapse=" "))
  }

  log_line("Refitting standardized and raw-unit trajectories")
  outcomes <- c("cortex_mm","log_hip_ratio")
  effects <- rbindlist(lapply(names(scopes),function(s)
    rbindlist(lapply(outcomes,function(y)fit_one(scopes[[s]],s,y)))),fill=TRUE)
  stopifnot(nrow(effects)==8L,all(is.finite(effects$estimate)),
            all(is.finite(effects$ci_low)),all(is.finite(effects$ci_high)))
  save_tsv(effects,"Table2_tau_slopes_standardized_and_raw.tsv")

  log_line("Auditing scanner, processing version, and QC coverage")
  qc <- rbindlist(lapply(names(scopes),function(s) {
    x <- copy(scopes[[s]])
    x[, `:=`(scope=s,
             scan_phase=ifelse(time_years==0,"baseline","follow-up"))]
    x[,.(scans=.N,participants=uniqueN(id)),
      by=.(cohort,scope,scan_phase,qc_tier,scanner,
           processing_version,stable_pipeline)]
  }),fill=TRUE)
  save_tsv(qc,"Table3_scan_processing_and_QC.tsv")
  scan_audit <- rbindlist(lapply(names(scopes),function(s) {
    x <- copy(scopes[[s]])
    x[,.(scope=s,cohort,id,mri_date,image_id,time_years,
          baseline_qc,qc_tier,baseline_scanner,scanner,
          baseline_version,processing_version,stable_pipeline,
          cortex_mm,hip_mean_mm3,tau_z)]
  }),fill=TRUE)
  save_tsv(scan_audit,"Frozen_scan_audit.tsv","analysis_data")

  # Descriptive annualized raw changes use each person's latest included scan.
  annual <- rbindlist(lapply(names(scopes),function(s) {
    x <- scopes[[s]]
    base <- x[time_years==0,.(id,cortex0=cortex_mm,hippo0=hip_mean_mm3)]
    follow <- x[time_years>0][order(id,-time_years)][, .SD[1L],by=id]
    z <- merge(base,follow[,.(id,time_years,cortex_mm,hip_mean_mm3)],by="id")
    z[,.(scope=s,cohort=unique(x$cohort),id,time_years,
         cortex_mm_per_year=(cortex_mm-cortex0)/time_years,
         hippo_percent_per_year=100*(hip_mean_mm3/hippo0-1)/time_years)]
  }),fill=TRUE)
  stopifnot(nrow(annual)==sum(flow$participants),
            all(is.finite(annual$cortex_mm_per_year)),
            all(is.finite(annual$hippo_percent_per_year)))
  save_tsv(annual,"Annualized_observed_change.tsv","analysis_data")
writeLines("STATUS=STEP16_COMPLETE",file.path(out,"16_COMPLETED.txt"))
})
