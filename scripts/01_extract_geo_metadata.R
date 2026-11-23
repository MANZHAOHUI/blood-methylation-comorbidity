# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 01: extract per-sample metadata of each data set
#
# Step 00 only gave total sample numbers; three things are still needed to finalise the design:
#   (1) the number of cases and controls in each data set (actual distribution of the group field)
#   (2) how many values the cell type field of GSE87648 really takes (whole blood or not)
#   (3) the age distribution of each set, to judge comparability of the PD and IBD cohorts
# The supplementary file names are also listed in full to confirm whether IDATs are really public.
#
# This script still does not download beta matrices; it only reads series matrix headers.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/01_extract_geo_metadata.R
# Dependency: GEOquery (Bioconductor)
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"

if (!dir.exists(PROJ_ROOT)) {
  stop("Project root does not exist: ", PROJ_ROOT)
}

OUT_DIR <- file.path(PROJ_ROOT, "output", "01_metadata")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GEOquery is missing. Please run first: ",
       'BiocManager::install("GEOquery")')
}

options(timeout = 900)

GSE_LIST <- c("GSE145361", "GSE111629", "GSE72774",
              "GSE87648",  "GSE112611",
              "GSE152027", "GSE84727",  "GSE80417")

MATRIX_SUFFIX <- c("_series_matrix.txt.gz",
                   "-GPL13534_series_matrix.txt.gz",
                   "-GPL21145_series_matrix.txt.gz",
                   "-GPL8490_series_matrix.txt.gz")

# --------------------------------------------------------------
# Read the header in chunks and stop at the data start marker, so that the matrix body is not read into memory
# --------------------------------------------------------------
fetch_header <- function(gse) {
  cache <- file.path(OUT_DIR, paste0(gse, "_header.txt"))
  if (file.exists(cache)) {
    return(readLines(cache, warn = FALSE, encoding = "UTF-8"))
  }
  stub <- sub("[0-9]{3}$", "nnn", gse)
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", stub, gse)
  for (u in paste0(base, gse, MATRIX_SUFFIX)) {
    lines <- tryCatch({
      con <- gzcon(url(u, open = "rb"))
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      acc <- character(0)
      repeat {
        chunk <- readLines(con, n = 20L, warn = FALSE)
        if (length(chunk) == 0L) break
        acc <- c(acc, chunk)
        if (any(grepl("!series_matrix_table_begin", chunk, fixed = TRUE))) break
        if (length(acc) >= 600L) break
      }
      acc
    }, error = function(e) NULL)
    if (!is.null(lines) && length(lines) > 0L) {
      hit <- which(grepl("!series_matrix_table_begin", lines, fixed = TRUE))
      if (length(hit)) lines <- lines[seq_len(hit[1])]
      writeLines(lines, cache, useBytes = TRUE)
      return(lines)
    }
  }
  stop("Cannot read the series matrix header of ", gse, "; please check the file names in its matrix directory")
}

split_fields <- function(line) {
  gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]])
}

# --------------------------------------------------------------
# Turn the !Sample_* header lines back into a per-sample table
# --------------------------------------------------------------
parse_samples <- function(lines) {
  acc_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
  if (length(acc_line) == 0L) stop("The header lacks the !Sample_geo_accession line")
  acc <- split_fields(acc_line[1])[-1]
  n <- length(acc)

  df <- data.frame(geo_accession = acc, stringsAsFactors = FALSE)

  title_line <- grep("^!Sample_title", lines, value = TRUE)
  if (length(title_line)) {
    v <- split_fields(title_line[1])[-1]
    df$sample_title <- if (length(v) == n) v else NA_character_
  }

  ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  for (k in seq_along(ch)) {
    v <- split_fields(ch[k])[-1]
    if (length(v) != n) v <- c(v, rep(NA_character_, n))[seq_len(n)]
    has_tag <- grepl(":", v, fixed = TRUE)
    tags <- ifelse(has_tag, trimws(sub(":.*$", "", v)), NA_character_)
    tab <- sort(table(tags[!is.na(tags)]), decreasing = TRUE)
    cname <- if (length(tab)) make.names(names(tab)[1]) else paste0("field_", k)
    while (cname %in% names(df)) cname <- paste0(cname, ".b")
    df[[cname]] <- ifelse(has_tag, trimws(sub("^[^:]*:", "", v)), v)
  }
  df
}

is_groupish <- function(x) {
  grepl("diagnos|disease|status|group|cell|type|tissue|case|control",
        x, ignore.case = TRUE)
}
is_ageish <- function(x) grepl("^age|age$|ageat", x, ignore.case = TRUE)

# --------------------------------------------------------------
# Main loop
# --------------------------------------------------------------
report <- character(0)
add <- function(...) report <<- c(report, paste0(...))

for (gse in GSE_LIST) {
  message("\n==== ", gse, " ====")
  lines <- fetch_header(gse)
  df <- parse_samples(lines)

  csv_path <- file.path(OUT_DIR, paste0(gse, "_sample_metadata.csv"))
  write.csv(df, csv_path, row.names = FALSE)

  add("### ", gse, "   samples = ", nrow(df))
  message("  Samples: ", nrow(df), ", fields: ", paste(names(df), collapse = ", "))

  for (nm in names(df)) {
    if (nm %in% c("geo_accession", "sample_title")) next
    v <- df[[nm]]
    if (is_ageish(nm)) {
      num <- suppressWarnings(as.numeric(v))
      if (sum(!is.na(num)) > 0) {
        add("  [age] ", nm,
            ": n=", sum(!is.na(num)),
            "  min=", round(min(num, na.rm = TRUE), 1),
            "  median=", round(median(num, na.rm = TRUE), 1),
            "  max=", round(max(num, na.rm = TRUE), 1))
        message("  [age] ", nm, " median=", round(median(num, na.rm = TRUE), 1))
        next
      }
    }
    if (is_groupish(nm) || length(unique(v)) <= 12L) {
      tb <- sort(table(v, useNA = "ifany"), decreasing = TRUE)
      tb <- tb[seq_len(min(12L, length(tb)))]
      add("  [", nm, "] ",
          paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; "))
      message("  [", nm, "] ",
              paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; "))
    }
  }

  supp <- tryCatch(
    GEOquery::getGEOSuppFiles(gse, makeDirectory = FALSE,
                              baseDir = OUT_DIR, fetch_files = FALSE),
    error = function(e) NULL)
  if (!is.null(supp)) {
    fn <- if (is.data.frame(supp)) {
      if ("fname" %in% names(supp)) as.character(supp$fname) else rownames(supp)
    } else as.character(supp)
    add("  [supp files] ", paste(fn, collapse = "; "))
    message("  Supplementary files: ", paste(fn, collapse = "; "))
  } else {
    add("  [supp files] could not be obtained")
  }
  add("")
}

rep_path <- file.path(OUT_DIR, "01_metadata_report.txt")
writeLines(report, rep_path, useBytes = TRUE)

message("\nSummary report: ", rep_path)
message("Per-sample tables: ", OUT_DIR, "  (one *_sample_metadata.csv per data set)")
message("\nNext: review 01_metadata_report.txt.")

