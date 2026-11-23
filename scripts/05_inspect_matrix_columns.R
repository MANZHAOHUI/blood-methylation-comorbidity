# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 05: check the column structure of the matrices
#
# Step 04 reported three sample counts that disagree with the metadata:
#   GSE111629  573 columns, 572 samples in the metadata
#   GSE112611  804 columns, 402 samples in the metadata (exactly twice)
#   GSE152027 1600 columns, 800 samples in the metadata (exactly twice)
# In addition, 43 columns of GSE72774 were classified as non-numeric.
#
# "Values between 0 and 1" does not prove that a column holds beta values - detection P values fall in the
# same range. This step therefore ignores the values and looks at the column names, to establish what the
# columns of each file really are before deciding which ones to take.
#
# The script only reads the header and the first 2000 rows; it is fast and uses little memory.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/05_inspect_matrix_columns.R
# ==============================================================

DATA_ROOT <- "D:/blood_methylation_data"
PROJ_ROOT <- "C:/blood_methylation_project"

if (!dir.exists(DATA_ROOT)) stop("Data directory does not exist: ", DATA_ROOT)
OUT_DIR <- file.path(PROJ_ROOT, "output", "05_columns")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

library(data.table)

FILES <- list(
  GSE87648  = "GSE87648/GSE87648_processedMethWb.txt.gz",
  GSE72774  = "GSE72774/GSE72774_datBetaNormalized.csv.gz",
  GSE111629 = "GSE111629/GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz",
  GSE112611 = "GSE112611/GSE112611_beta_values.txt.gz",
  GSE152027 = "GSE152027/GSE152027_IOP_processed_signals.csv.gz"
)

report <- character(0)
add <- function(...) report <<- c(report, paste0(...))

for (gse in names(FILES)) {
  path <- file.path(DATA_ROOT, FILES[[gse]])
  message("\n==== ", gse, " ====")
  if (!file.exists(path)) {
    message("  Source file does not exist, skipped"); add("### ", gse, "  Source file does not exist"); add(""); next
  }

  # ---- header ----
  con <- gzfile(path, "rt")
  hdr <- readLines(con, n = 1L, warn = FALSE)
  close(con)

  sep <- if (lengths(regmatches(hdr, gregexpr("\t", hdr))) >
             lengths(regmatches(hdr, gregexpr(",", hdr)))) "\t" else ","
  nm <- gsub('^"|"$', "", strsplit(hdr, sep, fixed = TRUE)[[1]])

  message("  Delimiter: ", if (sep == "\t") "tab" else "comma",
          "   header fields: ", length(nm))
  message("  First 8: ", paste(head(nm, 8), collapse = " | "))
  message("  Last 4: ", paste(tail(nm, 4), collapse = " | "))

  add("### ", gse)
  add("  Delimiter: ", if (sep == "\t") "tab" else "comma", "   header fields: ", length(nm))
  add("  First 8 fields: ", paste(head(nm, 8), collapse = " | "))
  add("  Last 4 fields: ", paste(tail(nm, 4), collapse = " | "))

  # ---- are there paired columns? ----
  pat <- c("[Dd]etection", "[Pp]val", "[Pp]_?value", "Signal_?[AB]",
           "[Uu]nmethylated", "[Mm]ethylated", "\\.1$", "_1$")
  for (p in pat) {
    k <- sum(grepl(p, nm))
    if (k > 0) {
      message("  Fields containing \"", p, "\": ", k, "")
      add("  Fields containing the pattern ", p, ": ", k, "")
    }
  }

  # Count unique stems after removing common suffixes, to decide whether each sample has two columns
  stem <- sub("[._-]?([Dd]etection[._ ]?[Pp]val.*|[Pp]val.*|[Dd]etection.*)$", "", nm)
  nstem <- length(unique(stem))
  message("  Unique stems after removing suffixes: ", nstem, " (fields in total ", length(nm), ")")
  add("  Unique stems after removing suffixes: ", nstem, " / fields in total ", length(nm))
  if (nstem > 1 && abs(nstem - length(nm) / 2) <= 2) {
    message("  >>> each sample has two columns; only one kind must be taken")
    add("  >>> Decision: two columns per sample, selection needed")
  }

  # ---- non-numeric columns in the first 2000 rows ----
  dt <- fread(path, nrows = 2000L, showProgress = FALSE)
  bad <- names(dt)[!vapply(dt, is.numeric, logical(1))]
  bad <- setdiff(bad, names(dt)[1])   # The first column is the probe ID and is character by nature
  message("  Non-numeric sample columns in the first 2000 rows: ", length(bad))
  add("  Non-numeric sample columns in the first 2000 rows: ", length(bad))
  if (length(bad)) {
    add("    column names (at most 10): ", paste(head(bad, 10), collapse = ", "))
    for (b in head(bad, 3)) {
      v <- as.character(dt[[b]])
      nonnum <- v[is.na(suppressWarnings(as.numeric(v))) & !is.na(v) & nzchar(v)]
      add("    ", b, ": examples of values that cannot be converted to numeric: ",
          paste(head(unique(nonnum), 5), collapse = " | "),
          "  ( ", length(unique(nonnum)), " distinct)")
      message("    ", b, " examples of abnormal values: ",
              paste(head(unique(nonnum), 5), collapse = " | "))
    }
  }
  rm(dt); invisible(gc())

  writeLines(nm, file.path(OUT_DIR, paste0(gse, "_colnames.txt")), useBytes = TRUE)
  add("")
}

rep_path <- file.path(OUT_DIR, "05_columns_report.txt")
writeLines(report, rep_path, useBytes = TRUE)
message("\nReport: ", rep_path)
message("The complete column names of each file were written separately for checking.")
message("\nNext: review 05_columns_report.txt.")
