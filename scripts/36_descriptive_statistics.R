# Summarize participant characteristics in the final analysis populations.
local({
library(data.table)
dest <- file.path(.study$output,"descriptive_statistics","run_frozen");dir.create(dest,recursive=TRUE,showWarnings=FALSE)
read_result <- function(folder,file) fread(file.path(.study$output,folder,"run_frozen",file))
put <- function(x,n) fwrite(as.data.table(x),file.path(dest,paste0(n,".tsv")),sep="\t",na="NA")
clinical <- read_result("cn_progression_62","tables/Table1_person_level_clinical_course.tsv")[prior_MCI_record==0]
cognition <- read_result("structural_history_70","source_data/MRI2_complete_model_inputs.tsv")[,.(cohort,id,baseline_score)]
clinical <- merge(clinical,cognition,by=c("cohort","id"),all.x=TRUE)
intervals <- read_result("clinical_field_audit_69","source_data/MRI_interval_in_final_risk_set.tsv")[,.(cohort,id,interval_years)]
apoe <- read_result("clinical_field_audit_69","source_data/Table4_person_level_APOE_availability.tsv")[,.(cohort,id,apoe_e4)]
stopifnot(!anyDuplicated(cognition,by=c("cohort","id")),
 !anyDuplicated(intervals,by=c("cohort","id")),!anyDuplicated(apoe,by=c("cohort","id")))
clinical <- merge(merge(clinical,intervals,by=c("cohort","id"),all.x=TRUE),apoe,by=c("cohort","id"),all.x=TRUE)
stopifnot(clinical[cohort=="Mayo",.N]==925L,clinical[cohort=="OASIS-3",.N]==306L)
fields <- c("age_MRI2","education","baseline_CDRSB","hippo2_pct_ICV","hippo_loss_pct_year","baseline_score","interval_years","impairment_time_years")
description <- rbindlist(lapply(fields,function(v)clinical[,.(variable=v,available_n=sum(is.finite(get(v))),mean=mean(get(v),na.rm=TRUE),SD=sd(get(v),na.rm=TRUE),median=as.numeric(median(get(v),na.rm=TRUE)),Q1=as.numeric(quantile(get(v),.25,na.rm=TRUE)),Q3=as.numeric(quantile(get(v),.75,na.rm=TRUE))),by=cohort]))
description[,unit:=fcase(variable%in%c("age_MRI2","education","interval_years","impairment_time_years"),"years",variable=="hippo2_pct_ICV","% ICV",variable=="hippo_loss_pct_year","%/year",variable=="baseline_score" & cohort=="Mayo","global z score",default="points")]
put(description,"clinical_continuous_characteristics")
put(clinical[,.(participants=.N,first_impairment_events=sum(impairment_event)),by=cohort],"clinical_participant_event_counts")
put(clinical[,.N,by=.(cohort,sex)],"clinical_sex_counts")
put(clinical[,.(participants=.N,genotyped=sum(apoe_e4%in%c(0,1)),carriers=sum(apoe_e4==1,na.rm=TRUE),carrier_percent=100*sum(apoe_e4==1,na.rm=TRUE)/sum(apoe_e4%in%c(0,1))),by=cohort],"clinical_apoe_counts")
put(read_result("prediagnostic_64","tables/Table5_updated_landmark_first_impairment.tsv")[,.(cohort,participants)],"reassessment_participant_counts")
tau <- fread(file.path(.study$output,"tau_hippocampal_models/run_frozen/analysis_data/Frozen_scan_audit.tsv"))[scope%in%c("ADNI_primary","HABS_primary")]
put(tau[,.(participants=uniqueN(id),scans=.N),by=cohort],"tau_participant_scan_counts")
base <- read_result("tau_longitudinal_inputs","analysis_data/Baseline_whole_hippocampus_scored.tsv")
base <- merge(base,unique(tau[,.(cohort,id)]),by=c("cohort","id"))
span <- tau[,.(structural_followup_years=max(time_years)),by=.(cohort,id)]
base <- merge(base,span,by=c("cohort","id"))
base[,hip_mean_pct_ICV:=100*hip_mean_mm3/icv_mm3]
put(rbindlist(lapply(c("age","sex_female","education","apoe4_carrier","baseline_mmse","hip_mean_pct_ICV","hip_mean_mm3","structural_followup_years","tau_raw","amyloid_raw"),function(v)base[,.(variable=v,available_n=sum(is.finite(get(v))),mean=mean(get(v),na.rm=TRUE),SD=sd(get(v),na.rm=TRUE),median=as.numeric(median(get(v),na.rm=TRUE)),Q1=as.numeric(quantile(get(v),.25,na.rm=TRUE)),Q3=as.numeric(quantile(get(v),.75,na.rm=TRUE))),by=cohort])),"tau_continuous_characteristics")
put(base[,.(participants=.N,female=sum(sex_female==1,na.rm=TRUE),female_percent=100*mean(sex_female==1,na.rm=TRUE),genotyped=sum(apoe4_carrier%in%c(0,1)),carriers=sum(apoe4_carrier==1,na.rm=TRUE),carrier_percent=100*sum(apoe4_carrier==1,na.rm=TRUE)/sum(apoe4_carrier%in%c(0,1))),by=cohort],"tau_categorical_characteristics")
})
