# Retrieve the public ENIGMA anatomical reference assets.
local({
commit <- "b08974b55243060cbc1fad12c87048037446e8f7"
dest <- file.path(.study$work,"inputs","ENIGMA_BigBrain")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
relative <- c("enigmatoolbox/histology/bb_gradient_aparc.csv",
 "enigmatoolbox/datasets/summary_statistics/allepi_case-controls_CortThick.csv")
for(rel in relative) {
 path <- file.path(dest,basename(rel))
 if(!file.exists(path)) download.file(paste0("https://raw.githubusercontent.com/MICA-MNI/ENIGMA/",commit,"/",rel),path,mode="wb")
}
})
