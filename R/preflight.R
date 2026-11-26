# Check required packages and cohort input files.
check_environment <- function() {
 p <- read.delim(file.path(.study$root,"docs/package-versions.tsv"),stringsAsFactors=FALSE)$package
 absent <- p[!vapply(p,requireNamespace,logical(1),quietly=TRUE)]
 if(length(absent)) stop("Missing packages: ",paste(absent,collapse=", "),". See README.md.")
 invisible(TRUE)
}
check_inputs <- function(stages=c("biology","clinical")) {
 x <- read.delim(file.path(.study$root,"docs/input-files.tsv"),stringsAsFactors=FALSE)
 x <- x[x$used_by%in%stages,]
 missing <- x$relative_path[!file.exists(file.path(.study$data,x$relative_path))]
 if(length(missing)) stop("Missing authorized cohort inputs under HIPPO_HISTORY_DATA:\n",paste(missing,collapse="\n"))
 invisible(TRUE)
}
