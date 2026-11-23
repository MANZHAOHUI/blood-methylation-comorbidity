# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 33a: look at what the inputs really are before changing 33
#
# 33 stopped at three places, each written on an assumption that did not hold:
#
#   (1) The probe universe came out as 0. Probe names were taken with rownames(anno), but this
#       annotation object is an S4 structure; colnames works (all three MAF columns were recognised),
#       but rownames returns no cg IDs - the probe names are most likely in a column, not in the row names.
#
#   (2) The column names of GSE84727 and GSE80417 match neither GSMs nor barcodes; both
#       counts are 0. These are the submitters' own normalisedBetas files with a third naming
#       scheme, and the actual column names were not printed, so the log does not show what they are.
#
#   (3) GSE198904 reported "adjacent column names differ". The repetition rate is exactly 0.50, so
#       columns are indeed paired, but the layout is not beta and P value alternating; more likely
#       the first half is all beta and the second half all P values. Only the alternating layout had been implemented.
#
# This step processes nothing and only prints three things as they are: the structure of the annotation object,
# the first column names and probe names of each data file, and how column names correspond to GSMs / barcodes.
# Changing 33 after seeing this is more reliable than guessing again.
#
# Only the first lines of each file are read; everything finishes within a minute or two.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/33a_diagnose_inputs.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "33a_diagnose_inputs"

RAW_DIR   <- file.path(DATA_ROOT, "raw")
PROBE_DIR <- file.path(DATA_ROOT, "probe_lists")

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
suppressPackageStartupMessages(library(data.table))

GSES <- c("GSE84727", "GSE80417", "GSE198904", "GSE273970", "GSE338167", "GSE144858")

# ==============================================================
say("=== 1. Structure of the annotation object ===")
af <- file.path(PROBE_DIR, "anno450k.rds")
if (!file.exists(af)) {
  say("  [FAIL] not found: ", af)
} else {
  anno <- readRDS(af)
  say("  Class: ", paste(class(anno), collapse = ", "))
  say("  Dimensions: ", paste(dim(anno), collapse = " x "))
  cn <- colnames(anno)
  say("  Column names (first 20): ", paste(utils::head(cn, 20), collapse = ", "))
  rn <- rownames(anno)
  say("  Row names present: ", !is.null(rn), "; first 5: ",
      if (is.null(rn)) "(none)" else paste(utils::head(rn, 5), collapse = ", "))
  # Find the column holding the cg IDs
  hit <- character(0)
  for (k in cn) {
    v <- tryCatch(as.character(anno[[k]]), error = function(e) character(0))
    if (!length(v)) next
    if (mean(grepl("^cg[0-9]{8}$", utils::head(v, 1000))) > 0.8) hit <- c(hit, k)
  }
  say("  Columns containing cg IDs: ", if (length(hit)) paste(hit, collapse = ", ") else "(none)")
  if (length(hit)) {
    v <- as.character(anno[[hit[1]]])
    say("    e.g. ", paste(utils::head(v, 5), collapse = ", "), "; total ", length(v), " rows")
  }
  # Chromosome column
  chrcol <- cn[grepl("^chr$|chromosome", cn, ignore.case = TRUE)]
  say("  Chromosome column: ", if (length(chrcol)) paste(chrcol, collapse = ", ") else "(none)")
  if (length(chrcol)) {
    tb <- table(as.character(anno[[chrcol[1]]]))
    say("    values (first 8): ", paste(sprintf("%s=%d", names(tb)[1:min(8, length(tb))],
                                            as.integer(tb)[1:min(8, length(tb))]),
                                    collapse = "; "))
  }
}
say("")

# ==============================================================
say("=== 2. Column names and probe names of each data file ===")
sl <- function(x) trimws(gsub('^"|"$', "", strsplit(x, "\t", fixed = TRUE)[[1]])[-1])
for (g in GSES) {
  say("---------- ", g, " ----------")
  hf <- file.path(RAW_DIR, g, paste0(g, "_series_matrix.txt.gz"))
  gsm <- ttl <- character(0)
  if (file.exists(hf)) {
    con <- gzfile(hf, "rt"); hd <- readLines(con, n = 400L, warn = FALSE); close(con)
    gsm <- sl(grep("^!Sample_geo_accession", hd, value = TRUE)[1])
    ttl <- sl(grep("^!Sample_title", hd, value = TRUE)[1])
    say("  Header: GSM ", length(gsm), " (", paste(utils::head(gsm, 2), collapse = ", "),
        "); title examples ", paste(utils::head(ttl, 2), collapse = ", "))
  } else say("  [NOTE] header file not found")

  fs <- list.files(file.path(RAW_DIR, g), pattern = "\\.(txt|csv)\\.gz$", full.names = TRUE)
  fs <- fs[file.info(fs)$size > 5 * 1024^2]
  if (!length(fs)) { say("  [NOTE] data file not found"); next }
  f <- fs[which.max(file.info(fs)$size)]
  dt <- data.table::fread(f, header = TRUE, check.names = FALSE,
                          nrows = 5L, showProgress = FALSE)
  cn <- names(dt)
  say("  Data file ", basename(f), ": ", ncol(dt), " columns")
  say("  Name of first column: ", cn[1], "; probe name examples: ",
      paste(utils::head(as.character(dt[[1]]), 3), collapse = ", "))
  say("  First 6 column names: ", paste(cn[2:min(7, length(cn))], collapse = " | "))
  say("  Last 3 column names: ", paste(utils::tail(cn, 3), collapse = " | "))

  d <- cn[-1]
  say("  Column names matching GSMs: ", sum(d %in% gsm), "; matching titles: ", sum(d %in% ttl), "")
  # If neither matches, try again after removing common prefixes and suffixes
  strip <- gsub("^X|\\.(AVG_Beta|Detection_Pval|beta|pval)$|_(Grn|Red)$", "", d,
                ignore.case = TRUE)
  say("  After removing common prefixes and suffixes: matching GSMs ", sum(strip %in% gsm),
      "; matching titles: ", sum(strip %in% ttl), "")

  ndup <- sum(duplicated(d))
  say("  Repeated column names: ", ndup, " ( ", sprintf("%.2f", mean(duplicated(d))), ")")
  if (ndup > 0) {
    n2 <- length(d) / 2
    interleaved <- all(d[seq(1, length(d), 2)] == d[seq(2, length(d), 2)])
    blocked <- length(d) %% 2 == 0 && all(d[1:n2] == d[(n2 + 1):length(d)])
    say("  Layout decision: alternating ", interleaved, "; in blocks ", blocked)
    if (!interleaved && !blocked)
      say("  [NOTE] neither layout fits; the complete column names need manual inspection")
  }
  rm(dt); invisible(gc())
  say("")
}
say("Diagnosis finished. Review this log before revising step 33.")
log_stop()
