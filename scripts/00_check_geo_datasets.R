# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 00: check the availability of the public data
#
# This script does not download beta matrices. It does four things:
#   (1) list the supplementary files (name, URL) of each GEO series and confirm that they are public
#   (2) determine whether raw IDAT / RAW.tar files are provided
#   (3) read the header of the series matrix and count the actual samples and platforms
#   (4) list the sample characteristic fields to see whether age, sex, smoking, medication and disease duration exist
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/00_check_geo_datasets.R
# Dependency: GEOquery (Bioconductor)
# ==============================================================

# Project root of this study (data, scripts and output all live below it)
PROJ_ROOT <- "C:/blood_methylation_project"

# The root directory must already exist. It is deliberately not created automatically: with a wrong path,
# recursive creation would silently build a wrong directory tree that is only noticed several steps later.
if (!dir.exists(PROJ_ROOT)) {
  stop("Project root does not exist: ", PROJ_ROOT,
       "\nPlease correct PROJ_ROOT at the top of this script.")
}

DATA_DIR <- file.path(PROJ_ROOT, "data")
OUT_DIR  <- file.path(PROJ_ROOT, "output", "00_data_check")

for (d in c(DATA_DIR, OUT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GEOquery is missing. Please run first: ",
       'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); ',
       'BiocManager::install("GEOquery")')
}

options(timeout = 600)

# --------------------------------------------------------------
# Data sets to be checked
# --------------------------------------------------------------
datasets <- data.frame(
  gse  = c("GSE145361", "GSE111629", "GSE72774",
           "GSE87648",  "GSE112611",
           "GSE152027", "GSE84727",  "GSE80417"),
  role = c("PD discovery", "PD replication 1", "PD replication 2",
           "IBD discovery", "IBD replication",
           "negative control (SCZ) 1", "negative control (SCZ) 2",
           "negative control (SCZ) 3"),
  stringsAsFactors = FALSE
)

# --------------------------------------------------------------
# List of supplementary files (the files themselves are not downloaded)
# --------------------------------------------------------------
list_supp <- function(gse) {
  res <- tryCatch(
    GEOquery::getGEOSuppFiles(gse, makeDirectory = FALSE,
                              baseDir = OUT_DIR, fetch_files = FALSE),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(ok = FALSE, msg = conditionMessage(res), files = character(0)))
  }
  fnames <- if (is.data.frame(res)) {
    if ("fname" %in% names(res)) as.character(res$fname) else rownames(res)
  } else {
    as.character(res)
  }
  list(ok = TRUE, msg = "", files = fnames, raw = res)
}

# --------------------------------------------------------------
# Read only the header of the series matrix
# --------------------------------------------------------------
read_matrix_header <- function(gse, max_lines = 300L) {
  stub <- sub("[0-9]{3}$", "nnn", gse)
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", stub, gse)
  # Large series are often split by platform; try each candidate file name
  candidates <- paste0(base, gse,
                       c("_series_matrix.txt.gz",
                         "-GPL13534_series_matrix.txt.gz",  # 450K
                         "-GPL21145_series_matrix.txt.gz",  # EPIC
                         "-GPL8490_series_matrix.txt.gz"))  # 27K
  for (u in candidates) {
    hdr <- tryCatch({
      con <- gzcon(url(u, open = "rb"))
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      readLines(con, n = max_lines, warn = FALSE)
    }, error = function(e) NULL)
    if (!is.null(hdr) && length(hdr) > 0L) {
      return(list(url = u, lines = hdr))
    }
  }
  NULL
}

clean_fields <- function(line) {
  f <- strsplit(line, "\t", fixed = TRUE)[[1]]
  gsub('^"|"$', "", f)
}

parse_header <- function(lines) {
  get_line <- function(tag) {
    i <- grep(paste0("^", tag), lines, fixed = FALSE)
    if (length(i) == 0L) NULL else lines[i]
  }

  acc_line <- get_line("!Sample_geo_accession")
  n_samp <- if (is.null(acc_line)) NA_integer_ else length(clean_fields(acc_line[1])) - 1L

  plat_line <- get_line("!Sample_platform_id")
  platform <- if (is.null(plat_line)) NA_character_ else {
    paste(unique(clean_fields(plat_line[1])[-1]), collapse = ";")
  }

  title_line <- get_line("!Series_title")
  title <- if (is.null(title_line)) NA_character_ else clean_fields(title_line[1])[2]

  ch_lines <- get_line("!Sample_characteristics_ch1")
  ch_tags <- character(0)
  if (!is.null(ch_lines)) {
    ch_tags <- vapply(ch_lines, function(l) {
      v <- clean_fields(l)[-1]
      v <- v[nzchar(v) & v != "NA"]
      if (length(v) == 0L) return(NA_character_)
      first <- v[1]
      if (grepl(":", first, fixed = TRUE)) {
        trimws(sub(":.*$", "", first))
      } else {
        paste0("(no tag) ", substr(first, 1, 30))
      }
    }, character(1), USE.NAMES = FALSE)
    ch_tags <- unique(ch_tags[!is.na(ch_tags)])
  }

  list(n_samples = n_samp, platform = platform, title = title,
       characteristics = ch_tags)
}

# --------------------------------------------------------------
# Main loop
# --------------------------------------------------------------
summary_rows <- list()
detail_lines <- character(0)

for (i in seq_len(nrow(datasets))) {
  gse  <- datasets$gse[i]
  role <- datasets$role[i]
  message("\n==== ", gse, "  (", role, ") ====")

  supp <- list_supp(gse)
  if (!supp$ok) {
    message("  Failed to obtain the list of supplementary files: ", supp$msg)
  } else {
    message("  Number of supplementary files: ", length(supp$files), "")
    for (f in supp$files) message("    - ", f)
  }
  has_idat <- any(grepl("idat|RAW\\.tar", supp$files, ignore.case = TRUE))

  hdr <- read_matrix_header(gse)
  if (is.null(hdr)) {
    message("  Failed to read the series matrix header (it may be split by platform under a non-standard file name; inspect the directory manually)")
    meta <- list(n_samples = NA_integer_, platform = NA_character_,
                 title = NA_character_, characteristics = character(0))
    matrix_url <- NA_character_
  } else {
    meta <- parse_header(hdr$lines)
    matrix_url <- hdr$url
    message("  Samples: ", meta$n_samples, "   Platform: ", meta$platform)
    message("  Sample characteristic fields: ",
            if (length(meta$characteristics)) paste(meta$characteristics, collapse = " | ")
            else "(none)")
  }

  summary_rows[[gse]] <- data.frame(
    gse             = gse,
    role            = role,
    n_samples       = meta$n_samples,
    platform        = meta$platform,
    n_supp_files    = length(supp$files),
    raw_idat_avail  = has_idat,
    matrix_url      = matrix_url,
    characteristics = paste(meta$characteristics, collapse = " | "),
    stringsAsFactors = FALSE
  )

  detail_lines <- c(detail_lines,
                    paste0("### ", gse, " (", role, ")"),
                    paste0("title: ", meta$title),
                    paste0("supp files: ", paste(supp$files, collapse = "; ")),
                    "")
}

res <- do.call(rbind, summary_rows)
rownames(res) <- NULL
print(res[, c("gse", "role", "n_samples", "platform",
              "raw_idat_avail", "n_supp_files")])

csv_path <- file.path(OUT_DIR, "geo_dataset_check.csv")
write.csv(res, csv_path, row.names = FALSE)
writeLines(detail_lines, file.path(OUT_DIR, "geo_dataset_supp_files.txt"),
           useBytes = TRUE)

message("\nCheck table written: ", csv_path)
message("Supplementary file details: ", file.path(OUT_DIR, "geo_dataset_supp_files.txt"))
message("\nNext: review geo_dataset_check.csv before finalising the study design.")

