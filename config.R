# Run from the repository root. Change these paths to your authorized local inputs.
.study <- list(
 root = normalizePath(getwd(),winslash="/",mustWork=TRUE),
 data = Sys.getenv("HIPPO_HISTORY_DATA",unset=file.path(getwd(),"data")),
 output = Sys.getenv("HIPPO_HISTORY_RESULTS",unset=file.path(getwd(),"results"))
)
.study$work <- file.path(.study$output,"reference_work")
# Keep seeds, formulas and inclusion windows in the analysis scripts unchanged.
