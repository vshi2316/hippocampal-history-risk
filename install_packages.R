# Install the packages used by the statistical analyses.
packages <- read.delim("docs/package-versions.tsv", stringsAsFactors = FALSE)$package
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
