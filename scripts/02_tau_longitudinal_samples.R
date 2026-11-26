# Construct MRI/PET analysis samples and longitudinal scan tables.
local({
CONFIG <- list(
  script_dir = .study$output,
  stage10_dir = file.path(.study$output,"tau_baseline_inputs/run_frozen"),
  raw_root = .study$data,
  seed = 20261022L,
  reference_fraction = 0.70,
  minimum_reference_n = 40L,
  followup_years = 5,
  minimum_MRI_interval_years = 0.5,
  clinical_window_days = 90L,
  pet_window_days = 180L,
  maximum_hippo_CT_gap_days = 0L
)
required <- c("data.table","readxl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages: install.packages(c(",
  paste(sprintf('"%s"', missing), collapse = ", "), "))")
suppressPackageStartupMessages(library(data.table))

num <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", ".", "-1", "-4", "-9", "-99",
              "-999", "-9999", "-8888", "-7777", "-777777")] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", "", z, fixed = TRUE)))
}
txt <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  z
}
col_text <- function(d, name) if (name %in% names(d)) txt(d[[name]]) else rep(NA_character_, nrow(d))
col_num <- function(d, name) if (name %in% names(d)) num(d[[name]]) else rep(NA_real_, nrow(d))
col_id <- function(d, name) sub("\\.0+$", "", col_text(d, name))
date_any <- function(x) {
  if (inherits(x, "Date")) return(as.Date(as.numeric(x), origin = "1970-01-01"))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  z <- txt(x); ans <- as.Date(rep(NA_real_, length(z)), origin = "1970-01-01")
  for (fmt in c("%Y-%m-%d","%Y/%m/%d","%m/%d/%Y","%d/%m/%Y")) {
    ii <- which(is.na(ans) & !is.na(z))
    if (length(ii)) ans[ii] <- as.Date(z[ii], format = fmt)
  }
  v <- suppressWarnings(as.numeric(z))
  ii <- which(is.na(ans) & is.finite(v) & v > 10000 & v < 90000)
  if (length(ii)) ans[ii] <- as.Date(v[ii], origin = "1899-12-30")
  as.Date(as.numeric(ans), origin = "1970-01-01")
}
month_date <- function(y, m) {
  y <- num(y); text_m <- txt(m); m <- suppressWarnings(as.integer(text_m))
  ii <- which(is.na(m) & !is.na(text_m))
  if (length(ii)) m[ii] <- match(tolower(substr(text_m[ii], 1, 3)), tolower(month.abb))
  ans <- as.Date(rep(NA_real_, length(y)), origin = "1970-01-01")
  ii <- which(is.finite(y) & y >= 1900 & y <= 2100 & is.finite(m) & m >= 1 & m <= 12)
  if (length(ii)) ans[ii] <- as.Date(sprintf("%04d-%02d-15", as.integer(y[ii]), m[ii]))
  as.Date(as.numeric(ans), origin = "1970-01-01")
}
positive_log <- function(x) {
  ans <- rep(NA_real_,length(x))
  ii <- which(is.finite(x) & x>0)
  if(length(ii)) ans[ii] <- log(x[ii])
  ans
}
safe_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_sd <- function(x) if (sum(is.finite(x)) > 1L) sd(x[is.finite(x)]) else NA_real_
safe_median <- function(x) if (any(is.finite(x))) as.numeric(median(x[is.finite(x)])) else NA_real_
p_text <- function(x) ifelse(is.finite(x), format.pval(x, digits = 3, eps = 0.0001), "")
norm_dx <- function(x) {
  z <- toupper(txt(x)); z <- sub("^.*[[:space:]]TO[[:space:]]", "", z)
  out <- rep(NA_character_, length(z))
  out[!is.na(z) & grepl("CN|NORMAL|CONTROL|SMC",z)] <- "CN"
  out[!is.na(z) & grepl("MCI|EMCI|LMCI",z)] <- "MCI"
  out[!is.na(z) & grepl("DEMENT|ALZ|^AD$",z)] <- "Dementia"
  out
}
qc_failed <- function(x) !is.na(x) & grepl("fail|reject|exclude|poor|invalid|hippocampus only",
                                          tolower(x), perl = TRUE)
is_pass <- function(x) !is.na(x) & tolower(trimws(x)) %in% c("pass","passed","included")
read_csv <- function(p) {
  if (!file.exists(p)) stop("Missing input: ", p)
  as.data.table(fread(p, na.strings = c("", "NA", "N/A", "NULL"), encoding = "UTF-8"))
}
read_inclusions <- function(p) {
  if (!file.exists(p)) stop("Missing input: ", p)
  sheets <- readxl::excel_sheets(p)
  selected <- sheets[grepl("inclusion|included", sheets, ignore.case = TRUE) &
                       !grepl("exclusion|excluded|pending|fail", sheets, ignore.case = TRUE)]
  if (!length(selected)) stop("No inclusion sheets found: ", basename(p))
  rows <- lapply(selected, function(sh) {
    d <- as.data.table(readxl::read_excel(p, sheet = sh, col_types = "text",
      na = c("", "NA", "N/A", "NULL"), .name_repair = "unique_quiet"))
    set(d, j = "source_sheet", value = rep(sh, nrow(d)))
    d
  })
  rbindlist(rows, fill = TRUE, use.names = TRUE)
}
date_max <- function(...) {
  z <- do.call(cbind, lapply(list(...), as.numeric))
  ans <- apply(z, 1, function(v) if (any(is.finite(v))) max(v[is.finite(v)]) else NA_real_)
  as.Date(ans, origin = "1970-01-01")
}
deduplicate <- function(d, keys) {
  d <- copy(d)
  setorderv(d, c(keys, intersect("source_row", names(d))), na.last = TRUE)
  unique(d, by = keys)
}
cohort_z <- function(x) {
  s <- safe_sd(x)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - safe_mean(x)) / s
}

OUT <- file.path(.study$output,"tau_longitudinal_inputs","run_frozen")
for (f in c("tables","models","analysis_data","audit")) dir.create(file.path(OUT,f),recursive=TRUE,showWarnings=FALSE)
save_tsv <- function(d, name, folder = "tables") {
  path <- file.path(OUT, folder, name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(as.data.table(d), path, sep = "\t", quote = TRUE, na = "NA", bom = TRUE)
  invisible(path)
}
log_line <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(...))
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = file.path(OUT, "11_run_log.txt"), append = TRUE)
}
new_status <- function(model, cohort, scope, outcome, status, n = 0L, participants = 0L) {
  data.table(model = model, cohort = cohort, scope = scope, outcome = outcome,
    term = NA_character_, estimate = NA_real_, se = NA_real_, conf_low = NA_real_,
    conf_high = NA_real_, p = NA_real_, q = NA_real_, n = as.integer(n),
    participants = as.integer(participants), focal = FALSE, status = status)
}

load_inputs <- function() {
  root <- CONFIG$raw_root
  raw <- list(
    fs = read_csv(file.path(root, "adni/MRI/MRI_FreeSurfer.csv")),
    ashs = read_csv(file.path(root, "adni/MRI/MRI_Hippocampus_ASHS.csv")),
    adni_clin = read_csv(file.path(root, "adni/Clinical/ADNIMERGE_03Jan2026.csv")),
    adni_tau = read_csv(file.path(root, "adni/PET/PET_Tau_Berkeley.csv")),
    adni_amy = read_csv(file.path(root, "adni/PET/PET_Amyloid_Berkeley.csv")),
    hct = read_inclusions(file.path(root, "habs_hd/HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx")),
    hhip = read_inclusions(file.path(root, "habs_hd/HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx")),
    htau = read_inclusions(file.path(root, "habs_hd/HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx")),
    hamy = read_inclusions(file.path(root, "habs_hd/HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx")),
    hdx = read_csv(file.path(root, "habs_hd/dx_12Aug2026.csv")),
    hmmse = read_csv(file.path(root, "habs_hd/mmse_12Aug2026.csv"))
  )
  base <- rbindlist(lapply(c("adni","habs"), function(co) {
    d <- read_csv(file.path(CONFIG$stage10_dir, paste0(co, "_baseline_anchored_unscored.tsv")))
    d[, diagnosis := norm_dx(diagnosis)]
    for (v in intersect(c("baseline_date","mri_date"), names(d))) set(d,j=v,value=date_any(d[[v]]))
    d
  }), fill = TRUE, use.names = TRUE)
  stopifnot(!anyDuplicated(base$id))
  list(raw = raw, base = base)
}

make_MRI <- function(raw) {
  regions <- c("entorhinal","parahippocampal","inferiortemporal","middletemporal")
  adni_fields <- c("ST24TA","ST44TA","ST32TA","ST40TA","ST83TA","ST103TA","ST91TA","ST99TA")
  habs_fields <- c(paste0("L_",regions,"_thickavg"),paste0("R_",regions,"_thickavg"))
  stopifnot(all(c(adni_fields,"ST29SV","ST88SV","ST10CV") %in% names(raw$fs)),
            all(habs_fields %in% names(raw$hct)))
  f <- raw$fs
  mat <- sapply(adni_fields, function(v) num(f[[v]]))
  ct <- rowMeans(mat)
  ct[!apply(is.finite(mat) & mat >= 1 & mat <= 5,1,all)] <- NA_real_
  a <- data.table(id = paste0("ADNI_",col_id(f,"RID")), cohort = "ADNI",
    mri_date = date_any(f$EXAMDATE), image_id = col_text(f,"IMAGEUID"),
    source_row = seq_len(nrow(f)), cortex_mm = ct,
    hip_left_mm3 = col_num(f,"ST29SV"), hip_right_mm3 = col_num(f,"ST88SV"),
    icv_mm3 = col_num(f,"ST10CV"), scanner = col_text(f,"FIELD_STRENGTH"),
    processing_version = col_text(f,"FSVER"), overall_qc = col_text(f,"OVERALLQC"),
    cortex_qc = col_text(f,"TEMPQC"), hippo_qc = col_text(f,"HIPPOQC"),
    source_status = col_text(f,"STATUS"), pipeline = "FreeSurfer whole hippocampus",
    date_precision = "day", hip_CT_gap_days = 0L)
  a[, explicit_fail := qc_failed(overall_qc) | qc_failed(cortex_qc) | qc_failed(hippo_qc)]
  a[, qc_tier := fifelse(explicit_fail, "EXPLICIT_FAIL",
    fifelse(is_pass(cortex_qc) & is_pass(hippo_qc), "REPORTED_REGIONAL_PASS", "QC_NOT_REPORTED"))]
  c <- raw$hct; h <- raw$hhip
  cm <- sapply(habs_fields, function(v) num(c[[v]]))
  ct <- rowMeans(cm)
  ct[!apply(is.finite(cm) & cm >= 1 & cm <= 5,1,all)] <- NA_real_
  hc <- data.table(id = paste0("HABS_",col_id(c,"Med_ID")), cohort = "HABS-HD",
    mri_date = month_date(c$MRI_year,c$MRI_month), cortex_mm = ct,
    image_id = col_text(c,"Time_point_MRI"), source_row = seq_len(nrow(c)),
    scanner = col_text(c,"Scanner_MRI"), source_status = col_text(c,"source_sheet"))
  hh <- data.table(id = paste0("HABS_",col_id(h,"Med_ID")),
    mri_date = month_date(h$MRI_year,h$MRI_month),
    hip_left_mm3 = col_num(h,"L_hippocampus"), hip_right_mm3 = col_num(h,"R_hippocampus"),
    icv_mm3 = col_num(h,"ICV_HD"), hip_scanner = col_text(h,"Scanner_MRI"),
    processing_version = col_text(h,"Software"), hip_visit = col_text(h,"Time_point"),
    hip_source_row = seq_len(nrow(h)))
  hc <- deduplicate(hc,c("id","mri_date","scanner"))
  hh <- deduplicate(hh,c("id","mri_date","hip_scanner"))
  hd <- merge(hc,hh,by=c("id","mri_date"),all.x=TRUE,allow.cartesian=TRUE,sort=FALSE)
  hd[, scanner_agreement := is.na(scanner) | is.na(hip_scanner) | scanner == hip_scanner]
  hd[, matched_hippo := is.finite(hip_left_mm3) & is.finite(hip_right_mm3) & scanner_agreement]
  setorderv(hd,c("id","mri_date","matched_hippo","hip_source_row"),c(1,1,-1,1),na.last=TRUE)
  hd <- unique(hd,by=c("id","mri_date","scanner"))
  hd[matched_hippo == FALSE, c("hip_left_mm3","hip_right_mm3","icv_mm3") :=
       list(NA_real_,NA_real_,NA_real_)]
  hd[, `:=`(explicit_fail = FALSE, qc_tier = "INCLUSION_RELEASE",
    overall_qc = "Inclusion sheet", cortex_qc = "Inclusion sheet", hippo_qc = "Inclusion sheet",
    pipeline = "HippoDeep whole hippocampus", date_precision = "month", hip_CT_gap_days = 0L)]
  mri <- rbindlist(list(a,hd), fill=TRUE,use.names=TRUE)
  mri[, hip_mean_mm3 := (hip_left_mm3 + hip_right_mm3)/2]
  mri[, valid_numeric := is.finite(cortex_mm) & cortex_mm >= 1 & cortex_mm <= 5 &
    is.finite(hip_left_mm3) & hip_left_mm3 >= 500 & hip_left_mm3 <= 8000 &
    is.finite(hip_right_mm3) & hip_right_mm3 >= 500 & hip_right_mm3 <= 8000 &
    is.finite(icv_mm3) & icv_mm3 >= 800000 & icv_mm3 <= 2500000]
  mri[, eligible_MRI := valid_numeric & !explicit_fail & !is.na(mri_date)]
  # Selection uses QC, completeness and source order, never a clinical outcome.
  mri[,qc_rank:=fifelse(qc_tier%in%c("REPORTED_REGIONAL_PASS","INCLUSION_RELEASE"),2L,1L)]
  setorderv(mri,c("id","mri_date","eligible_MRI","qc_rank","source_row"),c(1,1,-1,-1,1),na.last=TRUE)
  mri <- unique(mri,by=c("id","mri_date"))
  map <- data.table(cohort=rep(c("ADNI","HABS-HD"),each=11),
    target=rep(c(paste0(rep(c("L_","R_"),each=4),rep(regions,2)),
                 "left_whole_hippocampus","right_whole_hippocampus","intracranial_volume"),2),
    field=c(adni_fields,"ST29SV","ST88SV","ST10CV",habs_fields,
            "L_hippocampus","R_hippocampus","ICV_HD"))
  list(mri=mri,field_map=map)
}

make_clinical <- function(raw) {
  a <- raw$adni_clin
  adx <- data.table(id=paste0("ADNI_",col_id(a,"RID")),cohort="ADNI",
    visit_date=date_any(a$EXAMDATE),diagnosis=norm_dx(a$DX),mmse=col_num(a,"MMSE"),
    source_row=seq_len(nrow(a)))
  h <- raw$hdx
  code <- col_num(h,"CDX_COG")
  hdx <- data.table(id=paste0("HABS_",col_id(h,"MED_ID")),cohort="HABS-HD",
    visit_date=month_date(h$VISIT_YEAR,h$VISIT_MONTH),
    diagnosis=fifelse(code==0,"CN",fifelse(code==1,"MCI",fifelse(code==2,"Dementia",NA_character_))),
    source_row=seq_len(nrow(h)))
  hm <- raw$hmmse
  hmmse <- data.table(id=paste0("HABS_",col_id(hm,"MED_ID")),cohort="HABS-HD",
    visit_date=month_date(hm$VISIT_YEAR,hm$VISIT_MONTH),mmse=col_num(hm,"MMSE_TOTAL"),
    source_row=seq_len(nrow(hm)))
  dx <- rbindlist(list(adx[,.(id,cohort,visit_date,diagnosis,source_row)],hdx),fill=TRUE)
  dx <- deduplicate(dx[!is.na(visit_date) & !is.na(diagnosis)],c("id","visit_date"))
  mmse <- rbindlist(list(adx[,.(id,cohort,visit_date,mmse,source_row)],hmmse))
  mmse <- deduplicate(mmse[!is.na(visit_date) & is.finite(mmse) & mmse>=0 & mmse<=30],
                      c("id","visit_date"))
  list(dx=dx,mmse=mmse)
}

match_PET_dates <- function(base,raw) {
  d <- copy(base)
  adni <- which(d$cohort=="ADNI"); habs <- which(d$cohort=="HABS-HD")
  tau_date <- amy_date <- as.Date(rep(NA_real_,nrow(d)),origin="1970-01-01")
  tau_date[adni] <- date_any(raw$adni_tau$SCANDATE)[as.integer(d$tau_source_row[adni])]
  amy_date[adni] <- date_any(raw$adni_amy$SCANDATE)[as.integer(d$amy_source_row[adni])]
  tau_date[habs] <- month_date(raw$htau$PI2620_PET_year,raw$htau$PI2620_PET_month)[as.integer(d$ptau_source_row[habs])]
  amy_date[habs] <- month_date(raw$hamy$FBB_PET_year,raw$hamy$FBB_PET_month)[as.integer(d$fbb_source_row[habs])]
  set(d,j="tau_date",value=tau_date);set(d,j="amyloid_date",value=amy_date)
  # Confirm source-row identities before using inherited date matches.
  aid <- paste0("ADNI_",col_id(raw$adni_tau,"RID"))
  hid <- paste0("HABS_",col_id(raw$htau,"Med_ID"))
  good_a <- which(is.finite(d$tau_raw[adni]))
  good_h <- which(is.finite(d$tau_raw[habs]))
  stopifnot(all(d$id[adni[good_a]]==aid[as.integer(d$tau_source_row[adni[good_a]])]),
            all(d$id[habs[good_h]]==hid[as.integer(d$ptau_source_row[habs[good_h]])]))
  aid <- paste0("ADNI_",col_id(raw$adni_amy,"RID"))
  hid <- paste0("HABS_",col_id(raw$hamy,"Med_ID"))
  good_a <- which(is.finite(d$amyloid_raw[adni]))
  good_h <- which(is.finite(d$amyloid_raw[habs]))
  stopifnot(all(d$id[adni[good_a]]==aid[as.integer(d$amy_source_row[adni[good_a]])]),
            all(d$id[habs[good_h]]==hid[as.integer(d$fbb_source_row[habs[good_h]])]))
  d
}

rebuild_baseline <- function(base,mri,clinical,raw) {
  d <- match_PET_dates(base,raw)
  keep <- c("id","cohort","baseline_date","age","sex_female","education","apoe4_carrier",
            "diagnosis","tau_raw","amyloid_pos","amyloid_raw","tau_date","amyloid_date",
            "hippocampus_raw","ashs_gap_days")
  for(v in setdiff(keep,names(d))) set(d,j=v,value=rep(NA_real_,nrow(d)))
  d <- d[,..keep]
  setnames(d,"hippocampus_raw","ASHS_legacy_hip_mean_mm3")
  d[cohort!="ADNI",ASHS_legacy_hip_mean_mm3:=NA_real_]
  d[, mri_date := baseline_date]
  d <- merge(d,mri[,.(id,mri_date,cortex_mm,hip_left_mm3,hip_right_mm3,hip_mean_mm3,icv_mm3,
    scanner,processing_version,qc_tier,eligible_MRI,pipeline,date_precision,image_id)],by=c("id","mri_date"),
    all.x=TRUE,sort=FALSE)
  d[, eligible_baseline := eligible_MRI %in% TRUE & is.finite(age) & age>=50 & age<=90 &
       is.finite(sex_female)]
  near <- merge(d[,.(id,mri_date)],clinical$mmse,by="id",allow.cartesian=TRUE)
  near[, gap_days := as.numeric(visit_date-mri_date)]
  near <- near[abs(gap_days)<=CONFIG$clinical_window_days]
  near[,abs_gap:=abs(gap_days)];setorderv(near,c("id","abs_gap","visit_date"))
  near <- near[,.SD[1L],by=id][,.(id,baseline_mmse=mmse,baseline_mmse_date=visit_date)]
  d <- merge(d,near,by="id",all.x=TRUE,sort=FALSE)
  dd <- merge(d[,.(id,mri_date)],clinical$dx,by="id",allow.cartesian=TRUE)
  dd[,gap:=abs(as.numeric(visit_date-mri_date))]
  dd <- dd[gap<=CONFIG$clinical_window_days];setorderv(dd,c("id","gap","visit_date"))
  dd <- dd[,.SD[1L],by=id][,.(id,diagnosis_current=diagnosis,diagnosis_date=visit_date)]
  d <- merge(d,dd,by="id",all.x=TRUE,sort=FALSE)
  d[,diagnosis:=diagnosis_current]
  d[, `:=`(log_cortex=positive_log(cortex_mm),log_hip_ratio=positive_log(hip_mean_mm3/icv_mm3),
     age_c=(age-70)/10,education_c=(education-12)/4)]
  d[, `:=`(tau_z=cohort_z(tau_raw),amyloid_z=cohort_z(amyloid_raw)),by=cohort]
  d[, initial_landmark:=date_max(mri_date,diagnosis_date,amyloid_date)]
  d[, tau_landmark:=date_max(initial_landmark,tau_date)]
  d[, `:=`(tau_MRI_gap_days=as.numeric(tau_date-mri_date),
            amyloid_MRI_gap_days=as.numeric(amyloid_date-mri_date))]
  stopifnot(all(abs(d$tau_MRI_gap_days[is.finite(d$tau_raw)])<=CONFIG$pet_window_days,na.rm=TRUE),
            all(abs(d$amyloid_MRI_gap_days[is.finite(d$amyloid_raw)])<=CONFIG$pet_window_days,na.rm=TRUE))
  d
}

fit_reference <- function(base) {
  models <- list(); reference <- list(); params <- list()
  for(co in c("ADNI","HABS-HD")) {
    d <- copy(base[cohort==co & eligible_baseline==TRUE & diagnosis=="CN" & amyloid_pos==0])
    setorder(d,id)
    set.seed(CONFIG$seed+match(co,c("ADNI","HABS-HD")))
    ntrain <- floor(nrow(d)*CONFIG$reference_fraction)
    if(ntrain<CONFIG$minimum_reference_n) stop("Insufficient normative reference: ",co,"; n=",ntrain)
    train_ids <- sample(d$id,ntrain)
    d[,reference_partition:=ifelse(id%in%train_ids,"training","held_out_control")]
    reference[[co]] <- d[,.(id,cohort,reference_partition,age,sex_female)]
  }
  ref <- rbindlist(reference)
  d <- merge(base,ref[,.(id,reference_partition)],by="id",all.x=TRUE,sort=FALSE)
  d[is.na(reference_partition),reference_partition:="study_participant"]
  d
}
make_MRI_long <- function(base,mri,norms) {
  b <- base[eligible_baseline==TRUE & diagnosis%in%c("CN","MCI") & is.finite(tau_z) &
              !is.na(tau_date) & !is.na(amyloid_date) & amyloid_pos%in%0:1]
  keys <- c("id","cohort","baseline_date","age_c","sex_female","education_c","apoe4_carrier",
            "diagnosis","tau_z","amyloid_z","amyloid_pos","tau_landmark","icv_mm3",
            "scanner","processing_version","qc_tier","reference_partition")
  b <- b[,..keys]
  setnames(b,c("icv_mm3","scanner","processing_version","qc_tier"),
    c("baseline_icv_mm3","baseline_scanner","baseline_version","baseline_qc"))
  x <- merge(b,mri[eligible_MRI==TRUE,.(id,mri_date,cortex_mm,hip_mean_mm3,
    scanner,processing_version,qc_tier,image_id)],by="id",allow.cartesian=TRUE)
  x[,time_years:=as.numeric(mri_date-baseline_date)/365.25]
  x <- x[time_years==0 | (time_years>=CONFIG$minimum_MRI_interval_years &
     time_years<=CONFIG$followup_years & mri_date>tau_landmark)]
  # Scanner and software stability is a paired, within-participant restriction.
  x[,stable_pipeline:=scanner==baseline_scanner &
      ((is.na(processing_version)&is.na(baseline_version))|processing_version==baseline_version)]
  x[is.na(stable_pipeline),stable_pipeline:=FALSE]
  x[, `:=`(log_cortex=positive_log(cortex_mm),log_hip_ratio=positive_log(hip_mean_mm3/baseline_icv_mm3))]
  x[,scan_key:=paste(id,as.character(mri_date),sep="__")]
  x
}

input <- load_inputs()
raw <- input$raw
rebuilt <- make_MRI(raw)
mri <- rebuilt$mri
clinical <- make_clinical(raw)
base <- rebuild_baseline(input$base,mri,clinical,raw)
base <- fit_reference(base)
long <- make_MRI_long(base,mri,NULL)
save_tsv(base,"Baseline_whole_hippocampus_scored.tsv","analysis_data")
save_tsv(mri,"All_MRI_whole_hippocampus.tsv","analysis_data")
save_tsv(long,"MRI_longitudinal_analysis.tsv","analysis_data")
samples <- list()
for (co in c("ADNI","HABS-HD")) for (scope in c("Apositive_primary","Anegative_secondary",
 "Apositive_stable_pipeline","Apositive_reported_QC")) {
 d <- copy(long[cohort==co & reference_partition!="training"])
 ab <- if(scope=="Anegative_secondary") 0L else 1L
 d <- d[amyloid_pos==ab]
 if(scope=="Apositive_stable_pipeline") d <- d[stable_pipeline==TRUE]
 if(scope=="Apositive_reported_QC") d <- d[baseline_qc%in%c("REPORTED_REGIONAL_PASS","INCLUSION_RELEASE") &
  qc_tier%in%c("REPORTED_REGIONAL_PASS","INCLUSION_RELEASE")]
 eligible <- d[,.(n_dates=uniqueN(mri_date),has_baseline=any(time_years==0),has_future=any(time_years>0)),by=id]
 d <- d[id%in%eligible[n_dates>=2 & has_baseline & has_future,id]]
 covars <- c("sex_female","education_c","apoe4_carrier","diagnosis","scanner","processing_version")
 covars <- covars[vapply(covars,function(v) uniqueN(d[[v]][!is.na(d[[v]])])>1L,logical(1))]
 fields <- c("id","time_years","tau_z","amyloid_z","age_c","log_cortex",covars)
 cc <- complete.cases(d[,..fields])
 for(v in fields) if(is.numeric(d[[v]])) cc <- cc & is.finite(d[[v]])
 d <- d[which(cc)]
 for(tissue in c("cortex","hippocampus")) samples[[paste(co,scope,tissue,"tau_time",sep="__")]] <- list(data=d)
}
saveRDS(list(MRI=samples),file.path(OUT,"models","Longitudinal_models.rds"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"))
})
