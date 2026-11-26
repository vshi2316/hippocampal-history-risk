# Extract the released acquisition metadata tables.
local({
dest <- file.path(.study$output,"oasis_tau_mri_gate_27","run_frozen","audit","input_csv_cache")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
members <- c("MUSE_volume_OASIS_04Nov2026.csv","OASIS3_MR_scan_summary.csv")
archives <- list.files(file.path(.study$data,"oasis3"),"\\.zip$",full.names=TRUE)
for(n in members) {
 found <- FALSE
 for(z in archives) {
  names <- unzip(z,list=TRUE)$Name
  hit <- names[basename(names)==n]
  if(length(hit)==1L) {unzip(z,files=hit,exdir=dest,junkpaths=TRUE);found <- TRUE;break}
 }
 if(!found) stop("Missing archive member: ",n)
}
})
