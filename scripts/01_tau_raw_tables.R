local({
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

CONFIG <- list(
  seed = 20260831L,
  mri_pet_window_days = 180L,
  mri_cognition_window_days = 90L,
  min_ct_fields = 8L,
  min_habs_cognition_tests = 4L,
  infer_amyloid_from_centiloid = FALSE,
  centiloid_positive_cutoff = 20,
  out_dir = file.path(.study$output,"tau_baseline_inputs","run_frozen"),
  paths = list(
    adni_clin = file.path(.study$data,"adni/Clinical/ADNIMERGE_03Jan2026.csv"),
    adni_fs = file.path(.study$data,"adni/MRI/MRI_FreeSurfer.csv"),
    adni_fs_dict = file.path(.study$data,"adni/MRI/MRI_FreeSurfer_Data_Dictionary.csv"),
    adni_region_map = file.path(.study$root,"resources","ST_to_Brain_Region_Mapping.csv"),
    adni_ashs = file.path(.study$data,"adni/MRI/MRI_Hippocampus_ASHS.csv"),
    adni_amyloid = file.path(.study$data,"adni/PET/PET_Amyloid_Berkeley.csv"),
    adni_tau = file.path(.study$data,"adni/PET/PET_Tau_Berkeley.csv"),

    habs_ct = file.path(.study$data,"habs_hd/HABS_HD_RP_8_CorticalThickness_Data_BraskieLab_20260501.xlsx"),
    habs_hippo = file.path(.study$data,"habs_hd/HABS_HD_RP_8_Hippocampal_Volume_HippoDeep_Data_BraskieLab_20260501.xlsx"),
    habs_fbb = file.path(.study$data,"habs_hd/HABS_HD_RP_8_FBB_PET_data_BraskieLab_20260501.xlsx"),
    habs_tau = file.path(.study$data,"habs_hd/HABS_HD_RP_8_PI2620_PET_Data_BraskieLab_20260501.xlsx"),
    habs_genomics = file.path(.study$data,"habs_hd/RP_HD_7_Genomics.xlsx"),
    habs_demographics = file.path(.study$data,"habs_hd/interviewdemographics_12Aug2026.csv"),
    habs_dx = file.path(.study$data,"habs_hd/dx_12Aug2026.csv"),
    habs_zscores = file.path(.study$data,"habs_hd/zscores_12Aug2026.csv"),
    habs_mmse = file.path(.study$data,"habs_hd/mmse_12Aug2026.csv"),
    habs_cdr = file.path(.study$data,"habs_hd/consensuscdr_12Aug2026.csv"),
    habs_clinical_wide = file.path(.study$data,"habs_hd/RP_HD_7_Clinical.csv")
  )
)
dir.create(CONFIG$out_dir,recursive=TRUE,showWarnings=FALSE)
set.seed(CONFIG$seed)
required_packages <- c("data.table", "readxl", "lme4")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install required R packages before running: ",
    paste(missing_packages, collapse = ", "),
    "\nExample: install.packages(c(\"data.table\",\"readxl\",\"lme4\"))",
    call. = FALSE
  )
}
if (requireNamespace("lmerTest", quietly = TRUE)) {
  suppressPackageStartupMessages(library(lmerTest))
} else {
  suppressPackageStartupMessages(library(lme4))
}
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))

log_msg <- function(...) {
  msg <- paste0(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", file = file.path(CONFIG$out_dir, "10_analysis_log.txt"), append = TRUE, sep = "")
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_write <- function(x, file, sep = "\t") {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::as.data.table(x), file = file, sep = sep,
    quote = TRUE, na = "NA", bom = TRUE
  )
}

write_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = file, compress = "xz")
}

resolve_first_existing <- function(paths, label, required = FALSE) {
  paths <- unique(paths[!is.na(paths) & nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit)) return(hit[1L])
  if (required) {
    stop("Missing required input for ", label, ":\n", paste(paths, collapse = "\n"), call. = FALSE)
  }
  NA_character_
}

read_csv_dt <- function(path, required = FALSE) {
  if (is.na(path) || !file.exists(path)) {
    if (required) stop("Missing CSV: ", path, call. = FALSE)
    return(data.table())
  }
  data.table::fread(
    path, check.names = FALSE,
    na.strings = c("", "NA", "N/A", "NULL", "NaN", ".", "-1", "-4", "-9", "-99", "-999", "-9999", "-8888", "-777777")
  )
}

# Read only data/inclusion sheets. Dictionary/protocol/exclusion/pending sheets
# are excluded unless no usable sheet is found.
read_xlsx_union <- function(path, id_hint = NULL, required = FALSE) {
  if (is.na(path) || !file.exists(path)) {
    if (required) stop("Missing XLSX: ", path, call. = FALSE)
    return(data.table())
  }
  sheets <- readxl::excel_sheets(path)
  sheet_low <- tolower(sheets)
  inclusion_sheets <- sheets[grepl("inclusion|included", sheet_low) &
                               !grepl("exclusion|excluded|pending|fail", sheet_low)]
  preferred <- if (length(inclusion_sheets)) inclusion_sheets else sheets[
    !grepl("dictionary|codebook|protocol|readme|exclusion|excluded|pending|fail|qc report", sheet_low)
  ]
  if (!length(preferred)) preferred <- sheets
  read_one <- function(s) {
    z <- tryCatch(
      data.table::as.data.table(readxl::read_excel(
        path, sheet = s, col_types = "text", na = c("", "NA", "N/A", "NULL", "NaN"),
        .name_repair = "unique_quiet"
      )),
      error = function(e) data.table()
    )
    if (!nrow(z)) return(NULL)
    z[, source_sheet := s]
    z
  }
  dat <- Filter(Negate(is.null), lapply(preferred, read_one))
  if (!length(dat)) return(data.table())
  if (!is.null(id_hint)) {
    with_id <- dat[vapply(dat, function(z) id_hint %in% names(z), logical(1))]
    if (length(with_id)) dat <- with_id
  }
  data.table::rbindlist(dat, fill = TRUE, use.names = TRUE)
}

norm_name <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))

column_mapping_audit <- new.env(parent = emptyenv())
column_mapping_audit$rows <- list()
resolve_col <- function(x, candidates = character(), regex = NULL,
                        required = FALSE, label = "column") {
  nms <- names(x)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  selected <- NA_character_
  method <- "NOT_FOUND"
  n_matches <- 0L
  input_table <- paste(deparse(substitute(x)), collapse = "")
  on.exit({
    column_mapping_audit$rows[[length(column_mapping_audit$rows) + 1L]] <- data.table(
      input_table = input_table, label = label,
      candidates = paste(candidates, collapse = "|"),
      regex = if (is.null(regex)) NA_character_ else regex,
      selected = selected, method = method, n_matches = n_matches
    )
  }, add = TRUE)
  for (candidate in candidates) {
    hit <- which(nms == candidate)
    if (length(hit) == 1L) {
      selected <- nms[hit]; method <- "EXACT"; n_matches <- 1L
      return(selected)
    }
  }
  normalized <- norm_name(nms)
  for (candidate in candidates) {
    hit <- which(normalized == norm_name(candidate))
    if (length(hit) == 1L) {
      selected <- nms[hit]; method <- "NORMALIZED_UNIQUE"; n_matches <- 1L
      return(selected)
    }
    if (length(hit) > 1L) {
      method <- "AMBIGUOUS_NORMALIZED"; n_matches <- length(hit)
      if (required) stop("Ambiguous field: ", label, ": ", paste(nms[hit], collapse = ", "))
      return(NA_character_)
    }
  }
  if (!is.null(regex)) {
    hit <- grep(regex, nms, ignore.case = TRUE, perl = TRUE)
    n_matches <- length(hit)
    if (length(hit) == 1L) {
      selected <- nms[hit]; method <- "REGEX_UNIQUE"
      return(selected)
    }
    if (length(hit) > 1L) {
      method <- "AMBIGUOUS_REGEX"
      if (required) stop("Ambiguous field: ", label, ": ", paste(nms[hit], collapse = ", "))
      return(NA_character_)
    }
  }
  if (required) stop("Could not resolve ", label, ". Candidates: ", paste(candidates, collapse = ", "))
  NA_character_
}

get_num <- function(x, col, n = nrow(x)) {
  if (is.na(col) || !col %in% names(x)) return(rep(NA_real_, n))
  num_clean(x[[col]])
}

get_chr <- function(x, col, n = nrow(x)) {
  if (is.na(col) || !col %in% names(x)) return(rep(NA_character_, n))
  char_clean(x[[col]])
}

char_clean <- function(x) {
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", ".") |
       grepl("^-(1|4|9|99|999|9999|8888|777777)(\\.0+)?$", z)] <- NA_character_
  z
}

num_clean <- function(x) {
  z <- char_clean(x)
  z <- gsub(",", "", z, fixed = TRUE)
  suppressWarnings(as.numeric(z))
}

parse_date_any <- function(x) {
  if (inherits(x, "Date")) return(as.Date(x))
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) {
    out <- as.Date(x, origin = "1899-12-30")
    ok <- is.finite(x) & x > 1000 & x < 100000
    out[!ok] <- as.Date(NA)
    return(out)
  }
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "N/A", "NULL", ".", "-9999", "-8888")] <- NA_character_
  out <- as.Date(rep(NA_character_, length(z)))
  formats <- c(
    "%Y-%m-%d", "%Y/%m/%d", "%m/%d/%Y", "%d/%m/%Y",
    "%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S",
    "%m/%d/%Y %H:%M:%S", "%Y-%m-%dT%H:%M:%S"
  )
  for (fmt in formats) {
    ii <- is.na(out) & !is.na(z)
    if (any(ii)) out[ii] <- as.Date(z[ii], format = fmt)
  }
  serial <- suppressWarnings(as.numeric(z))
  ii <- is.na(out) & is.finite(serial) & serial > 1000 & serial < 100000
  if (any(ii)) out[ii] <- as.Date(serial[ii], origin = "1899-12-30")
  out
}

month_date <- function(year, month) {
  y <- num_clean(year)
  mchar <- trimws(as.character(month))
  m <- suppressWarnings(as.integer(mchar))
  miss <- is.na(m)
  if (any(miss)) m[miss] <- match(toupper(substr(mchar[miss], 1L, 3L)), toupper(month.abb))
  out <- as.Date(rep(NA_character_, length(y)))
  ok <- is.finite(y) & is.finite(m) & y >= 1900 & y <= 2100 & m >= 1 & m <= 12
  out[ok] <- as.Date(sprintf("%04d-%02d-15", as.integer(y[ok]), as.integer(m[ok])))
  out
}

standard_sex <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(z))
  # HABS codebook convention: 0=Male, 1=Female. Textual ADNI M/F labels are
  # handled in the same direction.
  out[grepl("^(female|f|1|woman|female sex)$", z)] <- "Female"
  out[grepl("^(male|m|0|man|male sex)$", z)] <- "Male"
  out
}

standard_dx_adni <- function(x) {
  z <- toupper(trimws(as.character(x)))
  # ADNIMERGE may encode transitions (e.g., "MCI to Dementia"); use the
  # destination state for the MRI-anchored visit.
  z <- sub("^.*[[:space:]]TO[[:space:]]", "", z)
  out <- rep(NA_character_, length(z))
  out[grepl("CN|NORMAL|CONTROL", z)] <- "CN"
  out[grepl("MCI|EMCI|LMCI|SMC", z)] <- "MCI"
  out[grepl("AD|DEMENT|ALZ", z)] <- "AD"
  out
}

binary_status <- function(x) {
  z <- toupper(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[grepl("^(pos|positive|abnormal|yes|true|present|1(?:\\.0+)?)($|[ _-])", z)] <- 1L
  out[grepl("^(neg|negative|normal|no|false|absent|0(?:\\.0+)?)($|[ _-])", z)] <- 0L
  out
}

qc_pass_generic <- function(x, cols = character()) {
  cols <- cols[cols %in% names(x)]
  if (!length(cols)) return(rep(TRUE, nrow(x)))
  txt <- apply(as.data.frame(x[, cols, with = FALSE]), 1L, paste, collapse = " ")
  bad <- grepl("fail|failed|reject|excluded|exclude|poor|bad|partial|not.?pass|invalid", tolower(txt), perl = TRUE)
  !bad
}

row_mean_required <- function(x, fields, min_required = length(fields)) {
  fields <- fields[fields %in% names(x)]
  if (!length(fields)) return(rep(NA_real_, nrow(x)))
  mat <- as.matrix(x[, fields, with = FALSE])
  storage.mode(mat) <- "numeric"
  n_ok <- rowSums(is.finite(mat))
  ans <- rowMeans(mat, na.rm = TRUE)
  ans[n_ok < min_required] <- NA_real_
  ans[!is.finite(ans)] <- NA_real_
  ans
}

canonical_id <- function(x, prefix = "") {
  z <- char_clean(x)
  z <- sub("\\.0+$", "", z)
  out <- if (nzchar(prefix)) paste0(prefix, z) else z
  out[is.na(z)] <- NA_character_
  out
}

collapse_same_id_date <- function(x, id_col = "id", date_col = "date",
                                  qc_col = "qc_pass", value_cols = NULL) {
  if (!nrow(x)) return(x)
  x <- data.table::as.data.table(x)
  x <- data.table::copy(x)
  data.table::setalloccol(x, ncol(x) + 4L)
  if (is.null(value_cols)) value_cols <- setdiff(names(x), c(id_col, date_col, qc_col, "source_row"))
  present_values <- intersect(value_cols, names(x))
  x[, .n_complete__ := if (length(present_values)) rowSums(!is.na(.SD)) else 0L,
    .SDcols = present_values]
  x[, .qc_rank__ := if (qc_col %in% names(x)) as.integer(get(qc_col)) else 1L]
  x[, .date_order__ := as.numeric(get(date_col))]
  x[is.na(.date_order__), .date_order__ := Inf]
  data.table::setorderv(x, c(id_col, date_col, ".qc_rank__", ".n_complete__", ".date_order__"),
                        c(1L, 1L, -1L, -1L, 1L), na.last = TRUE)
  key_value <- paste(as.character(x[[id_col]]), as.character(x[[date_col]]), sep = "|")
  x <- x[which(!duplicated(key_value))]
  x[, c(".n_complete__", ".qc_rank__", ".date_order__") := NULL]
  x
}

# Match the nearest recorded date and retain date precision and match status.
nearest_attach <- function(base, add, id_col = "id", base_date_col = "mri_date",
                           add_date_col = "date", max_days = 180L,
                           prefix = "add", add_cols = NULL,
                           allow_id_only = FALSE) {
  base <- data.table::copy(data.table::as.data.table(base))
  add <- data.table::copy(data.table::as.data.table(add))
  if (is.null(add_cols)) add_cols <- setdiff(names(add), c(id_col, add_date_col, "source_row"))
  add_cols <- add_cols[add_cols %in% names(add)]
  # Reserve data.table capacity for the appended match fields.
  data.table::setalloccol(base, ncol(base) + length(add_cols) + 4L)
  n <- nrow(base)
  idx_by_id <- if (nrow(add) && id_col %in% names(add)) split(seq_len(nrow(add)), as.character(add[[id_col]])) else list()
  chosen <- rep(NA_integer_, n)
  gaps <- rep(NA_real_, n)
  n_candidates <- integer(n)
  status <- rep("NO_ID_OR_NO_DATE", n)
  source_row <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    sid <- as.character(base[[id_col]][i])
    cand <- idx_by_id[[sid]]
    if (is.null(cand) || !length(cand)) next
    bdate <- base[[base_date_col]][i]
    adate <- add[[add_date_col]][cand]
    finite_date <- !is.na(adate)
    if (is.na(bdate)) {
      if (allow_id_only && sum(finite_date) == 0L && length(cand) == 1L) {
        chosen[i] <- cand[1L]; n_candidates[i] <- 1L
        status[i] <- "ID_ONLY_NO_DATE"; source_row[i] <- add$source_row[cand[1L]] %||% cand[1L]
      }
      next
    }
    if (!any(finite_date)) {
      if (allow_id_only && length(cand) == 1L) {
        chosen[i] <- cand[1L]; n_candidates[i] <- 1L
        status[i] <- "ID_ONLY_NO_DATE"; source_row[i] <- add$source_row[cand[1L]] %||% cand[1L]
      }
      next
    }
    gap <- abs(as.numeric(adate - bdate))
    keep <- which(is.finite(gap) & gap <= max_days)
    n_candidates[i] <- length(keep)
    if (!length(keep)) {
      status[i] <- "OUTSIDE_WINDOW"
      next
    }
    cand2 <- cand[keep]
    gap2 <- gap[keep]
    qc2 <- if ("qc_pass" %in% names(add)) as.integer(add$qc_pass[cand2]) else rep(1L, length(cand2))
    d2 <- as.numeric(add[[add_date_col]][cand2])
    sr2 <- if ("source_row" %in% names(add)) add$source_row[cand2] else cand2
    ord <- order(gap2, -qc2, d2, sr2, na.last = TRUE)
    pick <- cand2[ord[1L]]
    chosen[i] <- pick
    gaps[i] <- gap[match(pick, cand)]
    status[i] <- "MATCHED_DATE"
    source_row[i] <- if ("source_row" %in% names(add)) add$source_row[pick] else pick
  }
  for (cc in add_cols) {
    out_name <- paste0(prefix, "_", cc)
    data.table::set(base, j = out_name, value = add[[cc]][chosen])
  }
  data.table::set(base, j = paste0(prefix, "_gap_days"), value = gaps)
  data.table::set(base, j = paste0(prefix, "_n_candidates"), value = n_candidates)
  data.table::set(base, j = paste0(prefix, "_match_status"), value = status)
  data.table::set(base, j = paste0(prefix, "_source_row"), value = source_row)
  base
}

z_within <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}
z_by_group <- function(x, g) {
  out <- rep(NA_real_, length(x))
  g <- as.character(g)
  for (gg in unique(g[!is.na(g)])) {
    ii <- !is.na(g) & g == gg
    out[ii] <- z_within(x[ii])
  }
  out[is.na(g)] <- z_within(x[is.na(g)])
  out
}

valid_positive <- function(x) is.finite(x) & x > 0
valid_nonnegative <- function(x) is.finite(x) & x >= 0

file_hash <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_character_)
  # tools::md5sum is available in base R installations.
  as.character(tools::md5sum(path))
}


save_schema <- function(dat, source_name, path, out_file) {
  if (!nrow(dat)) {
    z <- data.table(source = source_name, path = path, exists = file.exists(path),
                    n_rows = 0L, n_cols = 0L, columns = NA_character_)
  } else {
    z <- data.table(source = source_name, path = path, exists = TRUE,
                    n_rows = nrow(dat), n_cols = ncol(dat),
                    columns = paste(names(dat), collapse = "|"))
  }
  safe_write(z, out_file)
  invisible(z)
}

region_map <- data.table(
  region = c("entorhinal", "parahippocampal", "inferiortemporal", "middletemporal"),
  adni_left = c("ST24TA", "ST44TA", "ST32TA", "ST40TA"),
  adni_right = c("ST83TA", "ST103TA", "ST91TA", "ST99TA"),
  habs_left = c("L_entorhinal_thickavg", "L_parahippocampal_thickavg",
                "L_inferiortemporal_thickavg", "L_middletemporal_thickavg"),
  habs_right = c("R_entorhinal_thickavg", "R_parahippocampal_thickavg",
                 "R_inferiortemporal_thickavg", "R_middletemporal_thickavg")
)
region_map[, `:=`(
  adni_l_name = paste0("L_", region, "_thick"),
  adni_r_name = paste0("R_", region, "_thick"),
  habs_l_name = paste0("L_", region, "_thick"),
  habs_r_name = paste0("R_", region, "_thick")
)]
safe_write(region_map, file.path(CONFIG$out_dir, "frozen_region_map.tsv"))
log_msg("Loading ADNI tables")
adni_clin_raw <- read_csv_dt(CONFIG$paths$adni_clin, required = TRUE)
adni_fs_raw <- read_csv_dt(CONFIG$paths$adni_fs, required = TRUE)
adni_ashs_raw <- read_csv_dt(CONFIG$paths$adni_ashs, required = TRUE)
adni_amy_raw <- read_csv_dt(CONFIG$paths$adni_amyloid, required = TRUE)
adni_tau_raw <- read_csv_dt(CONFIG$paths$adni_tau, required = TRUE)

schema_objects <- list(
  ADNI_clinical = adni_clin_raw,
  ADNI_FreeSurfer = adni_fs_raw,
  ADNI_ASHS = adni_ashs_raw,
  ADNI_amyloid = adni_amy_raw,
  ADNI_tau = adni_tau_raw
)
schema_rows <- rbindlist(lapply(names(schema_objects), function(nm) {
  z <- schema_objects[[nm]]
  data.table(source = nm, n_rows = nrow(z), n_cols = ncol(z),
             columns = paste(names(z), collapse = "|"))
}), fill = TRUE)
safe_write(schema_rows, file.path(CONFIG$out_dir, "schema_audit_adni.tsv"))

# Retain *_bl fields in the raw import. Later matching anchors clinical
# records to each MRI date.
adni_clin <- data.table(
  id = canonical_id(get_chr(adni_clin_raw, "RID"), "ADNI_"),
  rid = get_chr(adni_clin_raw, "RID"),
  ptid = get_chr(adni_clin_raw, "PTID"),
  date = parse_date_any(adni_clin_raw[[resolve_col(adni_clin_raw, c("EXAMDATE"), required = TRUE, label = "ADNI clinical EXAMDATE")]]),
  visit = get_chr(adni_clin_raw, resolve_col(adni_clin_raw, c("VISCODE2", "VISCODE"))),
  age_baseline = get_num(adni_clin_raw, "AGE"),
  years_from_baseline = get_num(adni_clin_raw, resolve_col(adni_clin_raw, c("Years_bl", "YEARS_BL"))),
  sex = standard_sex(get_chr(adni_clin_raw, "PTGENDER")),
  education = get_num(adni_clin_raw, "PTEDUCAT"),
  apoe4 = get_num(adni_clin_raw, "APOE4"),
  diagnosis_raw = get_chr(adni_clin_raw, resolve_col(adni_clin_raw, c("DX", "DX_bl"))),
  diagnosis = standard_dx_adni(get_chr(adni_clin_raw, resolve_col(adni_clin_raw, c("DX", "DX_bl")))),
  adas13 = get_num(adni_clin_raw, "ADAS13"),
  cdrsb = get_num(adni_clin_raw, "CDRSB"),
  mmse = get_num(adni_clin_raw, "MMSE"),
  faq = get_num(adni_clin_raw, "FAQ"),
  icv_clin = get_num(adni_clin_raw, "ICV"),
  site = get_chr(adni_clin_raw, "SITE"),
  source_row = seq_len(nrow(adni_clin_raw))
)
adni_clin[, `:=`(
  age = ifelse(is.finite(age_baseline) & is.finite(years_from_baseline),
               age_baseline + years_from_baseline, age_baseline),
  apoe4_carrier = ifelse(is.finite(apoe4), as.integer(apoe4 > 0), NA_integer_),
  cognition_primary = adas13,
  cognition_primary_direction = "higher_worse"
)]
adni_clin <- adni_clin[!is.na(id)]
adni_clin <- collapse_same_id_date(
  adni_clin,
  value_cols = c("age", "sex", "education", "apoe4", "diagnosis", "adas13", "cdrsb", "mmse", "faq")
)

# Decode explicit FreeSurfer thickness fields through the ADNI dictionary
# and DK anatomical mapping.
fs_dict <- read_csv_dt(CONFIG$paths$adni_fs_dict, required = TRUE)
fs_map <- read_csv_dt(CONFIG$paths$adni_region_map, required = TRUE)
if (!all(c("ST_Code", "Region_Name", "dk_region", "hemi") %in% names(fs_map))) {
  stop("ADNI region mapping must contain ST_Code, Region_Name, dk_region and hemi.", call. = FALSE)
}
fs_ta_map <- fs_map[grepl("CV$", ST_Code), .(
  dk_region = tolower(gsub("[^a-z0-9]", "", dk_region)),
  hemi = tolower(hemi), field = sub("CV$", "TA", ST_Code),
  mapped_name = Region_Name
)]
safe_write(fs_ta_map, file.path(CONFIG$out_dir, "adni_freesurfer_TA_mapping_audit.tsv"))
dict_field_col <- resolve_col(fs_dict, c("FLDNAME", "FIELD", "FIELD_NAME"), required = FALSE)
dict_text_col <- resolve_col(fs_dict, c("TEXT", "DESCRIPTION", "LABEL"), required = FALSE)
dict_audit <- data.table(
  field = c(region_map$adni_left, region_map$adni_right),
  in_data = c(region_map$adni_left, region_map$adni_right) %in% names(adni_fs_raw),
  in_dictionary = if (!is.na(dict_field_col)) c(region_map$adni_left, region_map$adni_right) %in% fs_dict[[dict_field_col]] else NA,
  dictionary_text = if (!is.na(dict_field_col) && !is.na(dict_text_col)) {
    as.character(fs_dict[[dict_text_col]][match(c(region_map$adni_left, region_map$adni_right), fs_dict[[dict_field_col]])])
  } else NA_character_
)
safe_write(dict_audit, file.path(CONFIG$out_dir, "adni_freesurfer_dictionary_audit.tsv"))

adni_fs <- data.table(
  id = canonical_id(get_chr(adni_fs_raw, "RID"), "ADNI_"),
  rid = get_chr(adni_fs_raw, "RID"),
  ptid = get_chr(adni_fs_raw, "PTID"),
  date = parse_date_any(adni_fs_raw[[resolve_col(adni_fs_raw, c("EXAMDATE"), required = TRUE, label = "ADNI FreeSurfer EXAMDATE")]]),
  visit = get_chr(adni_fs_raw, resolve_col(adni_fs_raw, c("VISCODE2", "VISCODE"))),
  field_strength = num_clean(sub("[[:space:]]*[Tt]$", "",
                                get_chr(adni_fs_raw, "FIELD_STRENGTH"))),
  fsver = get_chr(adni_fs_raw, "FSVER"),
  fs_status = tolower(char_clean(get_chr(adni_fs_raw, "STATUS"))),
  fs_overall_qc = get_chr(adni_fs_raw, "OVERALLQC"),
  qc_pass = {
    status0 <- tolower(char_clean(get_chr(adni_fs_raw, "STATUS")))
    explicit_fail <- !qc_pass_generic(adni_fs_raw, c("OVERALLQC", "TEMPQC", "FRONTQC", "PARQC", "INSULAQC", "OCCQC", "BGQC", "CWMQC", "VENTQC", "HIPPOQC"))
    # Legacy "partial" rows have no supplied QC grade; the MRI anchor below
    # additionally requires all eight prespecified cortical measurements.
    status_ok <- is.na(status0) | status0 %in% c("complete", "partial")
    status_ok & !explicit_fail
  },
  source_row = seq_len(nrow(adni_fs_raw))
)
for (j in seq_len(nrow(region_map))) {
  reg <- region_map$region[j]
  lf <- region_map$adni_left[j]
  rf <- region_map$adni_right[j]
  if (!all(c(lf, rf) %in% names(adni_fs_raw))) {
    stop("Required ADNI FreeSurfer thickness field missing for ", reg,
         ": ", lf, ", ", rf, call. = FALSE)
  }
  data.table::set(adni_fs, j = paste0("L_", reg, "_thick"), value = num_clean(adni_fs_raw[[lf]]))
  data.table::set(adni_fs, j = paste0("R_", reg, "_thick"), value = num_clean(adni_fs_raw[[rf]]))
}
ct_fields_adni <- unlist(lapply(region_map$region, function(r) {
  c(paste0("L_", r, "_thick"), paste0("R_", r, "_thick"))
}))
ct_matrix_adni <- as.matrix(adni_fs[, ct_fields_adni, with = FALSE])
storage.mode(ct_matrix_adni) <- "numeric"
adni_fs[, n_ct_fields := rowSums(is.finite(ct_matrix_adni))]
adni_fs[, ct_mtl_raw := row_mean_required(.SD, ct_fields_adni, CONFIG$min_ct_fields),
        .SDcols = ct_fields_adni]
adni_fs <- adni_fs[!is.na(id) & !is.na(date)]
adni_fs <- collapse_same_id_date(adni_fs, value_cols = ct_fields_adni)

adni_ashs <- data.table(
  id = canonical_id(get_chr(adni_ashs_raw, "RID"), "ADNI_"),
  rid = get_chr(adni_ashs_raw, "RID"),
  ptid = get_chr(adni_ashs_raw, "PTID"),
  date = parse_date_any(adni_ashs_raw[[resolve_col(adni_ashs_raw, c("EXAMDATE"), required = TRUE, label = "ADNI ASHS EXAMDATE")]]),
  left_hippo = get_num(adni_ashs_raw, "LEFT_HIPP_VOL"),
  right_hippo = get_num(adni_ashs_raw, "RIGHT_HIPP_VOL"),
  icv = get_num(adni_ashs_raw, "ICV"),
  qc_pass = qc_pass_generic(adni_ashs_raw, c("T2_QC_L", "T2_QC_R")),
  source_row = seq_len(nrow(adni_ashs_raw))
)
adni_ashs[, hippo_raw := ifelse(is.finite(left_hippo) & is.finite(right_hippo),
                                (left_hippo + right_hippo) / 2, NA_real_)]
adni_ashs[, hippo_ratio_raw := NA_real_]
adni_ashs[valid_positive(hippo_raw) & valid_positive(icv),
          hippo_ratio_raw := log(hippo_raw / icv)]
adni_ashs <- adni_ashs[!is.na(id) & !is.na(date)]
adni_ashs <- collapse_same_id_date(adni_ashs, value_cols = c("left_hippo", "right_hippo", "icv", "hippo_raw"))


adni_amy <- data.table(
  id = canonical_id(get_chr(adni_amy_raw, "RID"), "ADNI_"),
  rid = get_chr(adni_amy_raw, "RID"),
  ptid = get_chr(adni_amy_raw, "PTID"),
  date = parse_date_any(adni_amy_raw[[resolve_col(adni_amy_raw, c("SCANDATE"), required = TRUE, label = "ADNI amyloid SCANDATE")]]),
  centiloid = get_num(adni_amy_raw, "CENTILOIDS"),
  amyloid_suvr = get_num(adni_amy_raw, resolve_col(adni_amy_raw, c("SUMMARY_SUVR", "COMPOSITE_REF_SUVR"))),
  amyloid_status_raw = get_chr(adni_amy_raw, resolve_col(adni_amy_raw, c("AMYLOID_STATUS", "AMYLOID_STATUS_COMPOSITE_REF"))),
  amyloid_pos = binary_status(get_chr(adni_amy_raw, resolve_col(adni_amy_raw, c("AMYLOID_STATUS", "AMYLOID_STATUS_COMPOSITE_REF")))),
  tracer = get_chr(adni_amy_raw, "TRACER"),
  qc_pass = qc_pass_generic(adni_amy_raw, c("qc_flag", "TRACER_SUVR_WARNING")),
  source_row = seq_len(nrow(adni_amy_raw))
)
if (CONFIG$infer_amyloid_from_centiloid) {
  ii <- is.na(adni_amy$amyloid_pos) & is.finite(adni_amy$centiloid)
  adni_amy[ii, amyloid_pos := as.integer(centiloid >= CONFIG$centiloid_positive_cutoff)]
}
adni_amy[, amyloid_raw := ifelse(is.finite(centiloid), centiloid, amyloid_suvr)]
adni_amy[, amyloid_scale := ifelse(is.finite(centiloid), "Centiloid", "SUVR_fallback")]
adni_amy <- adni_amy[!is.na(id) & !is.na(date)]
adni_amy <- collapse_same_id_date(adni_amy, value_cols = c("centiloid", "amyloid_suvr", "amyloid_pos", "amyloid_raw"))

adni_tau <- data.table(
  id = canonical_id(get_chr(adni_tau_raw, "RID"), "ADNI_"),
  rid = get_chr(adni_tau_raw, "RID"),
  ptid = get_chr(adni_tau_raw, "PTID"),
  date = parse_date_any(adni_tau_raw[[resolve_col(adni_tau_raw, c("SCANDATE"), required = TRUE, label = "ADNI tau SCANDATE")]]),
  tau_meta_temporal = get_num(adni_tau_raw, "META_TEMPORAL_SUVR"),
  tau_entorhinal = get_num(adni_tau_raw, "CTX_ENTORHINAL_SUVR"),
  tracer = get_chr(adni_tau_raw, "TRACER"),
  qc_pass = qc_pass_generic(adni_tau_raw, c("TRACER_SUVR_WARNING")),
  source_row = seq_len(nrow(adni_tau_raw))
)
adni_tau[, tau_raw := ifelse(is.finite(tau_meta_temporal), tau_meta_temporal, tau_entorhinal)]
adni_tau <- adni_tau[!is.na(id) & !is.na(date) & is.finite(tau_raw)]
adni_tau <- collapse_same_id_date(adni_tau, value_cols = c("tau_raw", "tracer"))

# Filter failed QC records and retain exclusion counts alongside the raw import.
adni_ashs <- adni_ashs[qc_pass == TRUE]
adni_amy <- adni_amy[qc_pass == TRUE]
adni_tau <- adni_tau[qc_pass == TRUE]


# ----------------------------------------------------------------------------
# 04. ADNI MRI ANCHOR AND DETERMINISTIC NEAREST-DATE MATCHING
# ----------------------------------------------------------------------------
log_msg("Freezing ADNI baseline MRI anchors")
adni_mri_valid <- adni_fs[
  qc_pass == TRUE & is.finite(ct_mtl_raw) & n_ct_fields >= CONFIG$min_ct_fields
]
safe_write(adni_mri_valid[, .(images = .N, participants = uniqueN(id)),
                          by = .(fs_status, fs_overall_qc)],
           file.path(CONFIG$out_dir, "adni_MRI_QC_status_audit.tsv"))
data.table::setorderv(adni_mri_valid, c("id", "date", "source_row"), c(1L, 1L, 1L))
adni_mri_base <- adni_mri_valid[, .SD[1L], by = id]
adni_mri_base[, mri_date := date]
adni_mri_base[, mri_date_precision := "day"]

# Match clinical variables and each imaging modality independently to the
# MRI date, retaining the match intervals and status.
adni_base <- nearest_attach(
  adni_mri_base, adni_clin, max_days = CONFIG$mri_cognition_window_days,
  prefix = "clin", add_cols = c("age", "sex", "education", "apoe4", "apoe4_carrier",
                                "diagnosis_raw", "diagnosis", "adas13", "cdrsb", "mmse", "faq", "site", "source_row"),
  allow_id_only = FALSE
)
adni_base <- nearest_attach(
  adni_base, adni_ashs, max_days = CONFIG$mri_pet_window_days,
  prefix = "ashs", add_cols = c("left_hippo", "right_hippo", "hippo_raw", "icv", "hippo_ratio_raw", "source_row"),
  allow_id_only = FALSE
)
adni_base <- nearest_attach(
  adni_base, adni_amy, max_days = CONFIG$mri_pet_window_days,
  prefix = "amy", add_cols = c("centiloid", "amyloid_suvr", "amyloid_raw", "amyloid_pos", "amyloid_status_raw", "amyloid_scale", "tracer", "source_row"),
  allow_id_only = FALSE
)
adni_base <- nearest_attach(
  adni_base, adni_tau, max_days = CONFIG$mri_pet_window_days,
  prefix = "tau", add_cols = c("tau_raw", "tracer", "source_row"),
  allow_id_only = FALSE
)
adni_base <- data.table::copy(adni_base)
data.table::setalloccol(adni_base, ncol(adni_base) + 40L)

adni_base[, `:=`(
  cohort = "ADNI",
  participant_id = id,
  age = clin_age,
  sex = factor(clin_sex),
  education = clin_education,
  apoe4 = clin_apoe4,
  apoe4_carrier = clin_apoe4_carrier,
  diagnosis = factor(clin_diagnosis, levels = c("CN", "MCI", "AD")),
  diagnosis_raw = clin_diagnosis_raw,
  site = factor(clin_site),
  scanner = factor(ifelse(is.finite(field_strength), paste0("", field_strength, "T"), "unknown")),
  cognition_adas13 = clin_adas13,
  cognition_cdrsb = clin_cdrsb,
  cognition_mmse = clin_mmse,
  cognition_faq = clin_faq,
  hippocampus_raw = as.numeric(ashs_hippo_raw),
  hippocampus_icv = as.numeric(ashs_icv),
  hippocampus_ratio_raw = as.numeric(ashs_hippo_ratio_raw),
  amyloid_raw = as.numeric(amy_amyloid_raw),
  amyloid_pos = as.integer(amy_amyloid_pos),
  amyloid_scale = amy_amyloid_scale,
  amyloid_tracer = get("amy_tracer"),
  tau_raw = as.numeric(tau_tau_raw),
  tau_tracer = get("tau_tracer"),
  baseline_date = mri_date
)]
adni_base[, amyloid_pos := fifelse(!is.na(amy_amyloid_pos), as.integer(amy_amyloid_pos),
                                   fifelse(is.finite(amyloid_pos), as.integer(amyloid_pos), NA_integer_))]
adni_base[, `:=`(
  sex_female = ifelse(is.na(sex), NA_real_, as.numeric(sex == "Female")),
  scanner_group = factor(as.character(scanner))
)]

# Remove internal helper columns only after all match fields have been copied.
adni_base[, c("date") := NULL]
safe_write(adni_base, file.path(CONFIG$out_dir, "adni_baseline_anchored_unscored.tsv"))
habs_ct_raw <- read_xlsx_union(CONFIG$paths$habs_ct, id_hint = "Med_ID", required = TRUE)
habs_hippo_raw <- read_xlsx_union(CONFIG$paths$habs_hippo, id_hint = "Med_ID", required = TRUE)
habs_fbb_raw <- read_xlsx_union(CONFIG$paths$habs_fbb, id_hint = "Med_ID", required = TRUE)
habs_tau_raw <- read_xlsx_union(CONFIG$paths$habs_tau, id_hint = "Med_ID", required = TRUE)
habs_gen_raw <- read_xlsx_union(CONFIG$paths$habs_genomics, id_hint = "Med_ID", required = TRUE)
habs_dem_raw <- read_csv_dt(CONFIG$paths$habs_demographics, required = TRUE)
habs_dx_raw <- read_csv_dt(CONFIG$paths$habs_dx, required = TRUE)
habs_z_raw <- read_csv_dt(CONFIG$paths$habs_zscores, required = TRUE)
habs_mmse_raw <- read_csv_dt(CONFIG$paths$habs_mmse, required = TRUE)
habs_cdr_raw <- read_csv_dt(CONFIG$paths$habs_cdr, required = FALSE)


# HABS month/year visits use day 15 as a proxy. Retain month-level date precision.
habs_dem <- data.table(
  id = canonical_id(get_chr(habs_dem_raw, "MED_ID"), "HABS_"),
  med_id = get_chr(habs_dem_raw, "MED_ID"),
  visit = get_chr(habs_dem_raw, "VISIT_ID"),
  date = month_date(get_chr(habs_dem_raw, "VISIT_YEAR"), get_chr(habs_dem_raw, "VISIT_MONTH")),
  date_precision = "month",
  age = get_num(habs_dem_raw, "AGE"),
  sex = standard_sex(get_chr(habs_dem_raw, "ID_GENDER")),
  ethnicity = get_chr(habs_dem_raw, "ETHNICITY"),
  source_row = seq_len(nrow(habs_dem_raw))
)
habs_dem[, sex_female := ifelse(is.na(sex), NA_real_, as.numeric(sex == "Female"))]
habs_dem <- habs_dem[!is.na(id) & !is.na(date)]
habs_dem <- collapse_same_id_date(habs_dem, value_cols = c("age", "sex", "ethnicity"))

# Extract date-anchored education from the wide HABS clinical release when available.
habs_education <- data.table()
if (!is.na(CONFIG$paths$habs_clinical_wide) && file.exists(CONFIG$paths$habs_clinical_wide)) {
  habs_wide_raw <- read_csv_dt(CONFIG$paths$habs_clinical_wide, required = FALSE)
  edu_col <- resolve_col(habs_wide_raw,
                         c("ID_Education", "ID_Education_Degree", "Education", "EDUCATION"),
                         regex = "(^|_)EDUCATION($|_)", required = FALSE)
  if (!is.na(edu_col) && "MED_ID" %in% names(habs_wide_raw)) {
    habs_education <- data.table(
      id = canonical_id(get_chr(habs_wide_raw, "MED_ID"), "HABS_"),
      date = month_date(get_chr(habs_wide_raw, resolve_col(habs_wide_raw, c("VISIT_YEAR", "MRI_year"))),
                        get_chr(habs_wide_raw, resolve_col(habs_wide_raw, c("VISIT_MONTH", "MRI_month")))),
      education = get_num(habs_wide_raw, edu_col),
      source_row = seq_len(nrow(habs_wide_raw))
    )
    habs_education <- habs_education[!is.na(id)]
    habs_education <- collapse_same_id_date(habs_education, value_cols = c("education"))
  }
}

habs_dx <- data.table(
  id = canonical_id(get_chr(habs_dx_raw, "MED_ID"), "HABS_"),
  med_id = get_chr(habs_dx_raw, "MED_ID"),
  visit = get_chr(habs_dx_raw, "VISIT_ID"),
  date = month_date(get_chr(habs_dx_raw, "VISIT_YEAR"), get_chr(habs_dx_raw, "VISIT_MONTH")),
  date_precision = "month",
  diagnosis_code = get_num(habs_dx_raw, "CDX_COG"),
  source_row = seq_len(nrow(habs_dx_raw))
)
habs_dx <- habs_dx[!is.na(id) & !is.na(date)]
habs_dx <- collapse_same_id_date(habs_dx, value_cols = c("diagnosis_code"))

# Release-8 clinical codebook: 0=cognitively unimpaired, 1=MCI, 2=dementia,
# 9=undetermined. The raw code is retained alongside the harmonized label.
habs_dx[, diagnosis := fifelse(diagnosis_code == 0, "CN",
                        fifelse(diagnosis_code == 1, "MCI",
                        fifelse(diagnosis_code == 2, "Dementia", NA_character_)))]

# Construct the HABS cognitive composite from the specified release-7 z-score fields.
habs_cognitive_fields <- c(
  "DS_ZSCORE", "LM1_AB_ZSCORE", "LM2_AB_ZSCORE", "FAS_ZSCORE",
  "ANIMAL_ZSCORE", "DIGIT_SYMBOL_SUBSTITUTION_ZSCORE", "TRAILS_A_ZSCORE",
  "TRAILS_B_ZSCORE", "SEVLT_T1235_ZSCORE", "SEVLT_DR_ZSCORE"
)
habs_cognitive_fields <- habs_cognitive_fields[habs_cognitive_fields %in% names(habs_z_raw)]
if (length(habs_cognitive_fields) < CONFIG$min_habs_cognition_tests) {
  log_msg("WARNING: fewer than four prespecified HABS cognitive fields found; HABS cognition will be unavailable.")
}
habs_z <- data.table(
  id = canonical_id(get_chr(habs_z_raw, "MED_ID"), "HABS_"),
  med_id = get_chr(habs_z_raw, "MED_ID"),
  visit = get_chr(habs_z_raw, "VISIT_ID"),
  date = month_date(get_chr(habs_z_raw, "VISIT_YEAR"), get_chr(habs_z_raw, "VISIT_MONTH")),
  date_precision = "month",
  source_row = seq_len(nrow(habs_z_raw))
)
for (cc in habs_cognitive_fields) {
  data.table::set(habs_z, j = cc, value = num_clean(habs_z_raw[[cc]]))
}
habs_z[, n_cognitive_tests := if (length(habs_cognitive_fields)) {
  rowSums(is.finite(as.matrix(.SD)))
} else integer(.N), .SDcols = habs_cognitive_fields]
habs_z[, cognition_z := if (length(habs_cognitive_fields)) {
  row_mean_required(.SD, habs_cognitive_fields, CONFIG$min_habs_cognition_tests)
} else rep(NA_real_, .N), .SDcols = habs_cognitive_fields]
habs_z[, cognition_impairment := -cognition_z]
habs_z <- habs_z[!is.na(id) & !is.na(date)]
habs_z <- collapse_same_id_date(habs_z, value_cols = habs_cognitive_fields)

habs_mmse <- data.table(
  id = canonical_id(get_chr(habs_mmse_raw, "MED_ID"), "HABS_"),
  med_id = get_chr(habs_mmse_raw, "MED_ID"),
  visit = get_chr(habs_mmse_raw, "VISIT_ID"),
  date = month_date(get_chr(habs_mmse_raw, "VISIT_YEAR"), get_chr(habs_mmse_raw, "VISIT_MONTH")),
  date_precision = "month",
  mmse = get_num(habs_mmse_raw, "MMSE_TOTAL"),
  source_row = seq_len(nrow(habs_mmse_raw))
)
habs_mmse[mmse < 0 | mmse > 30, mmse := NA_real_]
habs_mmse <- habs_mmse[!is.na(id) & !is.na(date) & is.finite(mmse)]
habs_mmse <- collapse_same_id_date(habs_mmse, value_cols = c("mmse"))

habs_cdr <- data.table()
if (nrow(habs_cdr_raw)) {
  habs_cdr <- data.table(
    id = canonical_id(get_chr(habs_cdr_raw, "MED_ID"), "HABS_"),
    med_id = get_chr(habs_cdr_raw, "MED_ID"),
    visit = get_chr(habs_cdr_raw, "VISIT_ID"),
    date = month_date(get_chr(habs_cdr_raw, "VISIT_YEAR"), get_chr(habs_cdr_raw, "VISIT_MONTH")),
    date_precision = "month",
    cdr_global = get_num(habs_cdr_raw, resolve_col(habs_cdr_raw, c("CDR_GLOBAL", "CDR_GLOBAL_SCORE"))),
    cdr_sum = get_num(habs_cdr_raw, resolve_col(habs_cdr_raw, c("CDR_SUM", "CDR_SUM_OF_BOXES"))),
    source_row = seq_len(nrow(habs_cdr_raw))
  )
  habs_cdr <- habs_cdr[!is.na(id) & !is.na(date)]
  habs_cdr <- collapse_same_id_date(habs_cdr, value_cols = c("cdr_global", "cdr_sum"))
}

# HABS cortical thickness: keep all inclusion scans for longitudinal work.
# The exact Release-8 field names are checked dynamically, while the four
# primary regions are fixed in region_map.
if (!"Med_ID" %in% names(habs_ct_raw)) {
  stop("HABS cortical-thickness workbook has no Med_ID column.", call. = FALSE)
}
habs_ct <- data.table(
  id = canonical_id(get_chr(habs_ct_raw, "Med_ID"), "HABS_"),
  med_id = get_chr(habs_ct_raw, "Med_ID"),
  source_sheet = get_chr(habs_ct_raw, "source_sheet"),
  visit = get_chr(habs_ct_raw, resolve_col(habs_ct_raw, c("Time_point_MRI", "VISIT_ID", "Timepoint"))),
  date = month_date(get_chr(habs_ct_raw, "MRI_year"), get_chr(habs_ct_raw, "MRI_month")),
  date_precision = "month",
  age = get_num(habs_ct_raw, "Age"),
  scanner_raw = get_chr(habs_ct_raw, resolve_col(habs_ct_raw, c("Scanner_MRI", "Scanner", "MRI_Scanner"))),
  qc_pass = qc_pass_generic(habs_ct_raw, c("QC", "MRI_QC", "QC_status", "Status")),
  source_row = seq_len(nrow(habs_ct_raw))
)
for (j in seq_len(nrow(region_map))) {
  reg <- region_map$region[j]
  lf <- resolve_col(habs_ct_raw, c(region_map$habs_left[j]), required = TRUE,
                    label = paste0("HABS ", reg, " left thickness"))
  rf <- resolve_col(habs_ct_raw, c(region_map$habs_right[j]), required = TRUE,
                    label = paste0("HABS ", reg, " right thickness"))
  data.table::set(habs_ct, j = paste0("L_", reg, "_thick"), value = num_clean(habs_ct_raw[[lf]]))
  data.table::set(habs_ct, j = paste0("R_", reg, "_thick"), value = num_clean(habs_ct_raw[[rf]]))
}
ct_fields_habs <- unlist(lapply(region_map$region, function(r) {
  c(paste0("L_", r, "_thick"), paste0("R_", r, "_thick"))
}))
ct_matrix_habs <- as.matrix(habs_ct[, ct_fields_habs, with = FALSE])
storage.mode(ct_matrix_habs) <- "numeric"
habs_ct[, n_ct_fields := rowSums(is.finite(ct_matrix_habs))]
habs_ct[, ct_mtl_raw := row_mean_required(.SD, ct_fields_habs, CONFIG$min_ct_fields),
        .SDcols = ct_fields_habs]
habs_ct <- habs_ct[!is.na(id) & !is.na(date)]
habs_ct <- collapse_same_id_date(habs_ct, value_cols = ct_fields_habs)

habs_hippo <- data.table(
  id = canonical_id(get_chr(habs_hippo_raw, "Med_ID"), "HABS_"),
  med_id = get_chr(habs_hippo_raw, "Med_ID"),
  source_sheet = get_chr(habs_hippo_raw, "source_sheet"),
  visit = get_chr(habs_hippo_raw, resolve_col(habs_hippo_raw, c("Time_point", "Time_point_MRI", "VISIT_ID"))),
  date = month_date(get_chr(habs_hippo_raw, "MRI_year"), get_chr(habs_hippo_raw, "MRI_month")),
  date_precision = "month",
  age = get_num(habs_hippo_raw, "Age"),
  scanner_raw = get_chr(habs_hippo_raw, resolve_col(habs_hippo_raw, c("Scanner_MRI", "Scanner"))),
  left_hippo = get_num(habs_hippo_raw, resolve_col(habs_hippo_raw, c("L_hippocampus", "LEFT_HIPPOCAMPUS"), required = TRUE, label = "HABS left hippocampus")),
  right_hippo = get_num(habs_hippo_raw, resolve_col(habs_hippo_raw, c("R_hippocampus", "RIGHT_HIPPOCAMPUS"), required = TRUE, label = "HABS right hippocampus")),
  icv = get_num(habs_hippo_raw, resolve_col(habs_hippo_raw, c("ICV_HD", "ICV"), required = TRUE, label = "HABS hippocampal ICV")),
  source_row = seq_len(nrow(habs_hippo_raw))
)
habs_hippo[, hippo_raw := ifelse(is.finite(left_hippo) & is.finite(right_hippo),
                                 (left_hippo + right_hippo) / 2, NA_real_)]
habs_hippo[, hippo_ratio_raw := NA_real_]
habs_hippo[valid_positive(hippo_raw) & valid_positive(icv),
           hippo_ratio_raw := log(hippo_raw / icv)]
habs_hippo <- habs_hippo[!is.na(id) & !is.na(date)]
habs_hippo <- collapse_same_id_date(habs_hippo, value_cols = c("left_hippo", "right_hippo", "icv", "hippo_raw"))

habs_fbb <- data.table(
  id = canonical_id(get_chr(habs_fbb_raw, "Med_ID"), "HABS_"),
  med_id = get_chr(habs_fbb_raw, "Med_ID"),
  visit = get_chr(habs_fbb_raw, resolve_col(habs_fbb_raw, c("Time_point_PET_FBB", "VISIT_ID", "Time_point"))),
  date = month_date(get_chr(habs_fbb_raw, "FBB_PET_year"), get_chr(habs_fbb_raw, "FBB_PET_month")),
  date_precision = "month",
  centiloid = get_num(habs_fbb_raw, resolve_col(habs_fbb_raw, c("Centiloid_Units_GlobalAmyloid", "CENTILOIDS"))),
  amyloid_suvr = get_num(habs_fbb_raw, resolve_col(habs_fbb_raw, c("Global_FBB_SUVR", "FBB_SUVR"))),
  amyloid_pos = binary_status(get_chr(habs_fbb_raw, resolve_col(habs_fbb_raw, c("AB_pos_FBB", "AMYLOID_STATUS")))),
  mri_match = get_chr(habs_fbb_raw, "MRI_match"),
  months_between_mri_pet = get_num(habs_fbb_raw, "Months_between_MRI_PET"),
  tracer = "FBB",
  source_row = seq_len(nrow(habs_fbb_raw))
)
habs_fbb[, amyloid_raw := ifelse(is.finite(centiloid), centiloid, amyloid_suvr)]
habs_fbb[, amyloid_scale := ifelse(is.finite(centiloid), "Centiloid", "FBB_SUVR")]
habs_fbb[, qc_pass := TRUE]
habs_fbb <- habs_fbb[!is.na(id) & !is.na(date)]
habs_fbb <- collapse_same_id_date(habs_fbb, value_cols = c("centiloid", "amyloid_suvr", "amyloid_pos", "amyloid_raw"))

habs_tau <- data.table(
  id = canonical_id(get_chr(habs_tau_raw, "Med_ID"), "HABS_"),
  med_id = get_chr(habs_tau_raw, "Med_ID"),
  visit = get_chr(habs_tau_raw, resolve_col(habs_tau_raw, c("Time_point_PET_PI2620", "Time_point_PET", "Time_point_PI2620", "VISIT_ID", "Time_point"))),
  date = month_date(get_chr(habs_tau_raw, resolve_col(habs_tau_raw, c("PI2620_PET_year", "PET_year", "Year"))),
                    get_chr(habs_tau_raw, resolve_col(habs_tau_raw, c("PI2620_PET_month", "PET_month", "Month")))),
  date_precision = "month",
  tau_raw = get_num(habs_tau_raw, resolve_col(habs_tau_raw,
                                              c("Median_Medial_Temporal_PI2620_SUVR",
                                                "Tennant_MsTL_taucutpoint_weighted_SUVR"),
                                              regex = "Medial.*Temporal.*SUVR|MsTL.*SUVR",
                                              required = TRUE, label = "HABS PI2620 medial temporal tau")),
  mri_match = get_chr(habs_tau_raw, "MRI_match"),
  tracer = "PI-2620",
  source_row = seq_len(nrow(habs_tau_raw))
)
habs_tau[, qc_pass := TRUE]
habs_tau <- habs_tau[!is.na(id) & !is.na(date) & is.finite(tau_raw)]
habs_tau <- collapse_same_id_date(habs_tau, value_cols = c("tau_raw"))

habs_gen <- data.table(
  id = canonical_id(get_chr(habs_gen_raw, "Med_ID"), "HABS_"),
  med_id = get_chr(habs_gen_raw, "Med_ID"),
  apoe4_carrier = binary_status(get_chr(habs_gen_raw, resolve_col(habs_gen_raw, c("APOE4_Positivity", "APOE4_Carrier")))),
  apoe4_raw = get_num(habs_gen_raw, resolve_col(habs_gen_raw, c("APOE4_Positivity", "APOE4_Carrier"))),
  education = get_num(habs_gen_raw, resolve_col(habs_gen_raw, c("ID_Education", "Education", "EDUCATION"))),
  age = get_num(habs_gen_raw, resolve_col(habs_gen_raw, c("Age", "AGE"))),
  sex = standard_sex(get_chr(habs_gen_raw, resolve_col(habs_gen_raw, c("ID_Gender", "Gender", "SEX")))),
  source_row = seq_len(nrow(habs_gen_raw))
)
habs_gen <- habs_gen[!is.na(id)]
habs_gen[is.na(apoe4_carrier) & is.finite(apoe4_raw), apoe4_carrier := as.integer(apoe4_raw > 0)]
habs_gen[, valid_gen := as.integer(is.finite(apoe4_raw) | !is.na(apoe4_carrier))]
data.table::setorderv(habs_gen, c("id", "valid_gen", "source_row"), c(1L, -1L, 1L))
habs_gen <- habs_gen[!duplicated(id)]
habs_gen[, valid_gen := NULL]


# ----------------------------------------------------------------------------
# 06. HABS MRI ANCHOR AND DATE MATCHING
# ----------------------------------------------------------------------------
log_msg("Freezing HABS-HD baseline MRI anchors")
habs_ct_valid <- habs_ct[
  qc_pass == TRUE & is.finite(ct_mtl_raw) & n_ct_fields >= CONFIG$min_ct_fields
]
scan1 <- habs_ct_valid[
  grepl("scan[ _-]*1|time[ _-]*point[ _-]*1|^1$", tolower(paste(source_sheet, visit)), perl = TRUE)
]
if (!nrow(scan1)) scan1 <- habs_ct_valid
data.table::setorderv(scan1, c("id", "date", "source_row"), c(1L, 1L, 1L))
habs_mri_base <- scan1[, .SD[1L], by = id]
habs_mri_base[, `:=`(mri_date = date, mri_date_precision = "month")]

habs_base <- nearest_attach(
  habs_mri_base, habs_dem, max_days = CONFIG$mri_cognition_window_days,
  prefix = "dem", add_cols = c("age", "sex", "sex_female", "ethnicity", "date_precision", "source_row"),
  allow_id_only = FALSE
)
habs_base <- nearest_attach(
  habs_base, habs_dx, max_days = CONFIG$mri_cognition_window_days,
  prefix = "dx", add_cols = c("diagnosis_code", "diagnosis", "date_precision", "source_row"),
  allow_id_only = FALSE
)
habs_base <- nearest_attach(
  habs_base, habs_z, max_days = CONFIG$mri_cognition_window_days,
  prefix = "cog", add_cols = c("cognition_z", "cognition_impairment", "n_cognitive_tests", "date_precision", "source_row"),
  allow_id_only = FALSE
)
habs_base <- nearest_attach(
  habs_base, habs_mmse, max_days = CONFIG$mri_cognition_window_days,
  prefix = "mmse", add_cols = c("mmse", "date_precision", "source_row"),
  allow_id_only = FALSE
)
if (nrow(habs_cdr)) {
  habs_base <- nearest_attach(
    habs_base, habs_cdr, max_days = CONFIG$mri_cognition_window_days,
    prefix = "cdr", add_cols = c("cdr_global", "cdr_sum", "date_precision", "source_row"),
    allow_id_only = FALSE
  )
} else {
  habs_base[, `:=`(cdr_cdr_global = NA_real_, cdr_cdr_sum = NA_real_,
                   cdr_gap_days = NA_real_, cdr_n_candidates = 0L,
                   cdr_match_status = "NOT_FOUND", cdr_source_row = NA_integer_)]
}
habs_base <- nearest_attach(
  habs_base, habs_hippo, max_days = CONFIG$mri_pet_window_days,
  prefix = "hip", add_cols = c("left_hippo", "right_hippo", "hippo_raw", "icv", "hippo_ratio_raw", "source_row"),
  allow_id_only = FALSE
)
habs_base <- nearest_attach(
  habs_base, habs_fbb, max_days = CONFIG$mri_pet_window_days,
  prefix = "fbb", add_cols = c("centiloid", "amyloid_suvr", "amyloid_raw", "amyloid_pos", "amyloid_scale", "mri_match", "months_between_mri_pet", "source_row"),
  allow_id_only = FALSE
)
habs_base <- nearest_attach(
  habs_base, habs_tau, max_days = CONFIG$mri_pet_window_days,
  prefix = "ptau", add_cols = c("tau_raw", "tracer", "mri_match", "source_row"),
  allow_id_only = FALSE
)
if (nrow(habs_education)) {
  habs_base <- nearest_attach(
    habs_base, habs_education, max_days = CONFIG$mri_cognition_window_days,
    prefix = "edu", add_cols = c("education", "source_row"), allow_id_only = FALSE
  )
} else {
  habs_base[, `:=`(edu_education = NA_real_, edu_gap_days = NA_real_,
                   edu_n_candidates = 0L, edu_match_status = "NOT_FOUND",
                   edu_source_row = NA_integer_)]
}
habs_base <- merge(
  habs_base,
  habs_gen[, .(id, gen_apoe4_carrier = apoe4_carrier,
               gen_education = education, gen_age = age, gen_sex = sex)],
  by = "id", all.x = TRUE, sort = FALSE
)
habs_base <- data.table::copy(habs_base)
data.table::setalloccol(habs_base, ncol(habs_base) + 40L)


# Coalesce direct CT values with matched demographic values, preserving the
# MRI-table value when available. HABS diagnosis uses the Release-7 codebook
# mapping (0=CN, 1=MCI, 2=dementia; 9=undetermined is kept missing).
habs_base[, `:=`(
  cohort = "HABS-HD",
  participant_id = id,
  age_final = ifelse(is.finite(dem_age), dem_age, ifelse(is.finite(gen_age), gen_age, age)),
  sex_final = ifelse(!is.na(dem_sex), dem_sex, gen_sex),
  sex_female = ifelse(!is.na(dem_sex_female), dem_sex_female, NA_real_),
  education = ifelse(is.finite(edu_education), edu_education, gen_education),
  apoe4_carrier = gen_apoe4_carrier,
  diagnosis_code = dx_diagnosis_code,
  diagnosis = factor(dx_diagnosis, levels = c("CN", "MCI", "Dementia")),
  site = factor(NA_character_),
  scanner_group = factor(as.character(scanner_raw)),
  cognition_z_final = cog_cognition_z,
  cognition_impairment_final = cog_cognition_impairment,
  cognition_mmse = mmse_mmse,
  cognition_cdr_sum = cdr_cdr_sum,
  hippocampus_raw = hip_hippo_raw,
  hippocampus_icv = hip_icv,
  hippocampus_ratio_raw = hip_hippo_ratio_raw,
  amyloid_raw = fbb_amyloid_raw,
  amyloid_pos = as.integer(fbb_amyloid_pos),
  amyloid_scale = fbb_amyloid_scale,
  tau_raw = ptau_tau_raw,
  tau_tracer = get("ptau_tracer"),
  baseline_date = mri_date
)]
habs_base[, `:=`(age=age_final,sex=factor(sex_final),cognition_impairment=cognition_impairment_final)]
habs_base[is.na(sex_female) & !is.na(gen_sex), sex_female := as.numeric(gen_sex == "Female")]
habs_base[, c("date") := NULL]
safe_write(habs_base, file.path(CONFIG$out_dir, "habs_baseline_anchored_unscored.tsv"))

})
