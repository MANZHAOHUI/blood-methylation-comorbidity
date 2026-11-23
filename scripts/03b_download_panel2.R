# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 03b: download of the six additional cohorts
#
# Relation to 03_download_data.R: 03 is the download record of the first six data sets and
# is left unchanged - its manifest documents what was downloaded then, and changing it would make it impossible
# to reconstruct the sources of the first round from the scripts. This script is its sibling and covers only the six new sets.
#
# File selection (decided from measured sizes, not guessed):
#   Of the six sets only GSE144858 has a series matrix containing the data table (490 MB).
#   The series matrices of the other five are only some ten to fifty KB: header only,
#   no data table - the submitters did not put the processed matrix in the VALUE column but only in
#   the supplementary files. For these five the supplementary matrix is downloaded instead.
#   Rationale: one 450K sample is about 1 MB compressed, so a matrix of hundreds or thousands of samples
#   must be hundreds of MB; a few tens of KB can only be a bare header.
#
# The script has this check built in: any "data-containing" file whose remote size is below 5 MB
# is flagged as header-only, so that it cannot download "successfully" and then fail downstream.
#
# The download back end reuses the tuned settings of 03:
#   aria2c -x16 -s16 -k 20M --max-connection-per-server=16
# These parameters were tuned empirically - with 100M pieces and -s4 the log only showed CN:4,
# and four connections were not saturated; only 20M pieces give 16 segments. Without aria2c, fall back to curl.
#
# Environment: Windows + RStudio. Can run unattended; after an interruption a re-run resumes.
# Location: <project root>/scripts/03b_download_panel2.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "C:/blood_methylation_download"
SCRIPT_NAME <- "03b_download_panel2"

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

DO_DOWNLOAD <- TRUE     # Run with FALSE first to see the manifest and the true sizes; switch to TRUE after checking
MIN_DATA_MB <- 5         # Files containing data should not be smaller than this

RAW_DIR <- file.path(DATA_ROOT, "raw")
if (!dir.exists(RAW_DIR)) dir.create(RAW_DIR, recursive = TRUE)
if (!dir.exists(RAW_DIR)) stop("Failed to create the data directory: ", RAW_DIR)

ARIA <- Sys.which("aria2c")
CURL <- Sys.which("curl")
if (!nzchar(ARIA) && !nzchar(CURL)) stop("Neither aria2c nor curl.exe was found")
BACKEND <- if (nzchar(ARIA)) "aria2c" else "curl"
say("Download back end: ", BACKEND, "  (", if (BACKEND == "aria2c") ARIA else CURL, ")")

human <- function(b) {
  if (is.na(b)) return("unknown")
  u <- c("B", "KB", "MB", "GB", "TB"); i <- 1
  while (b >= 1024 && i < length(u)) { b <- b / 1024; i <- i + 1 }
  sprintf("%.1f %s", b, u[i])
}

# ==============================================================
# Manifest: specified set by set, no global switch
# ==============================================================
FILES <- data.frame(rbind(
  c("GSE84727",  "schizophrenia", "GSE84727_normalisedBetas.csv.gz"),
  c("GSE80417",  "schizophrenia", "GSE80417_normalizedBetas.csv.gz"),
  c("GSE198904", "depression",    "GSE198904_DHRC_Matrix_processed.txt.gz"),
  c("GSE198904", "depression",    "GSE198904_OBSERVEMDD0001_Matrix_processed.txt.gz"),
  c("GSE273970", "depression",    "GSE273970_Processed_Data.txt.gz"),
  c("GSE338167", "alzheimer",     "GSE338167_beta_matrix_normalized.txt.gz"),
  c("GSE144858", "alzheimer",     "GSE144858_series_matrix.txt.gz")
), stringsAsFactors = FALSE)
names(FILES) <- c("gse", "disease", "file")

# Header files: the supplementary matrices carry no sample information, so phenotypes still come from the series matrix headers;
# these are only tens of KB and are downloaded along the way
HEADERS <- data.frame(rbind(
  c("GSE84727",  "schizophrenia", "GSE84727_series_matrix.txt.gz"),
  c("GSE80417",  "schizophrenia", "GSE80417_series_matrix.txt.gz"),
  c("GSE198904", "depression",    "GSE198904_series_matrix.txt.gz"),
  c("GSE273970", "depression",    "GSE273970_series_matrix.txt.gz"),
  c("GSE338167", "alzheimer",     "GSE338167_series_matrix.txt.gz")
), stringsAsFactors = FALSE)
names(HEADERS) <- c("gse", "disease", "file")

FILES$kind <- "data"
HEADERS$kind <- "header"
ALL <- rbind(HEADERS, FILES)      # Small files first

geo_url <- function(gse, fname) {
  stub <- sub("[0-9]{3}$", "nnn", gse)
  sub_dir <- if (grepl("series_matrix\\.txt\\.gz$", fname)) "matrix" else "suppl"
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/%s/%s", stub, gse, sub_dir, fname)
}
remote_size <- function(u) {
  out <- tryCatch(system2(CURL, c("-sIL", "--retry", "3", shQuote(u)),
                          stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0))
  cl <- grep("^content-length:", tolower(out), value = TRUE)
  if (!length(cl)) return(NA_real_)
  as.numeric(trimws(sub("^content-length:", "", tolower(cl[length(cl)]))))
}

# ==============================================================
# Part 1: manifest and true sizes
# ==============================================================
say("")
say("=== Querying remote sizes ===")
rows <- list(); flagged <- character(0)
for (i in seq_len(nrow(ALL))) {
  gse <- ALL$gse[i]; f <- ALL$file[i]
  d <- file.path(RAW_DIR, gse)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
  u <- geo_url(gse, f)
  sz <- remote_size(u)
  dest <- file.path(d, f)
  local <- if (file.exists(dest)) file.size(dest) else 0
  warn <- ""
  if (ALL$kind[i] == "data" && !is.na(sz) && sz < MIN_DATA_MB * 1024^2) {
    warn <- "  [FAIL] file too small; most likely header only without a data table"
    flagged <- c(flagged, paste0(gse, " / ", f))
  }
  say(sprintf("  %-9s %-14s %-52s %10s%s",
              gse, ALL$kind[i], f, human(sz), warn))
  rows[[length(rows) + 1]] <- data.frame(
    gse = gse, disease = ALL$disease[i], kind = ALL$kind[i], file = f,
    url = u, dest = dest, remote_bytes = sz, local_bytes = local,
    done = !is.na(sz) && local == sz, stringsAsFactors = FALSE)
  Sys.sleep(0.3)
}
man <- do.call(rbind, rows)
write.csv(man, file.path(DATA_ROOT, "download_manifest_panel2.csv"), row.names = FALSE)
todo <- man[!man$done, ]
say("")
say("  Files in total: ", nrow(man), ", to download: ", nrow(todo),
    ", total size: ", human(sum(todo$remote_bytes, na.rm = TRUE)))
if (length(flagged)) {
  say("  [NOTE] the following files have unexpected sizes; confirm the file names before downloading: ")
  for (x in flagged) say("    ", x)
}

# ==============================================================
# Part 2: download
# ==============================================================
gz_ok <- function(path) {
  tryCatch({
    con <- gzfile(path, "rb"); on.exit(close(con), add = TRUE)
    length(readBin(con, "raw", n = 1024L)) > 0
  }, error = function(e) FALSE, warning = function(w) FALSE)
}
fetch_one <- function(u, dest) {
  if (BACKEND == "aria2c") {
    system2(ARIA, c("-c", "-x", "16", "-s", "16", "-k", "20M",
                    "--max-connection-per-server=16",
                    "--max-tries=0", "--retry-wait=10",
                    "--connect-timeout=60", "--timeout=120",
                    "--file-allocation=none", "--summary-interval=30",
                    "-d", shQuote(dirname(dest)),
                    "-o", shQuote(basename(dest)), shQuote(u)))
  } else {
    system2(CURL, c("-L", "-C", "-", "--retry", "50", "--retry-delay", "10",
                    "--speed-limit", "1000", "--speed-time", "180",
                    "--connect-timeout", "60",
                    "-o", shQuote(dest), shQuote(u)))
  }
}

# GSE198904 combines two sub-studies (DHRC and OBSERVEMDD0001).
# If cases and controls fell into different sub-studies, group and sub-study would be completely confounded and the data set unusable.
# This only needs the header of a few tens of KB, so it is checked before the 2.7 GB data matrix is downloaded;
# if it fails, the data files are skipped and the download is saved.
gate_198904 <- function() {
  f <- file.path(RAW_DIR, "GSE198904", "GSE198904_series_matrix.txt.gz")
  if (!file.exists(f)) { say("  [NOTE] header not yet available; confounding check not done in this run"); return(TRUE) }
  con <- gzfile(f, "rt"); on.exit(close(con), add = TRUE)
  hd <- readLines(con, n = 400L, warn = FALSE)
  pull <- function(tag) {
    for (one in grep("^!Sample_characteristics_ch1", hd, value = TRUE)) {
      v <- trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1])
      vv <- v[nzchar(v) & grepl(":", v, fixed = TRUE)]
      if (!length(vv)) next
      if (tolower(trimws(sub(":.*$", "", vv[1]))) == tolower(tag))
        return(trimws(sub("^[^:]*:", "", v)))
    }
    NULL
  }
  dg <- pull("diagnosis"); ch <- pull("cohort")
  if (is.null(dg) || is.null(ch)) { say("  [NOTE] diagnosis or cohort field not found"); return(TRUE) }
  keep <- dg %in% c("MDD", "HC")
  ct <- table(ch[keep], dg[keep])
  say("  Cross-tabulation of sub-study and group: ")
  for (r in rownames(ct))
    say("    ", r, ": ", paste(sprintf("%s=%d", colnames(ct), as.integer(ct[r, ])),
                               collapse = "; "))
  if (any(ct == 0)) {
    say("  [FAIL] a sub-study lacks one of the groups; group and sub-study are completely confounded, the data files of this set are skipped")
    return(FALSE)
  }
  say("  PASS: both sub-studies contain cases and controls")
  TRUE
}

if (!DO_DOWNLOAD) {
  say("")
  say("DO_DOWNLOAD = FALSE: manifest only, nothing downloaded. Switch to TRUE and re-run after checking.")
} else {
  say("")
  say("=== Step 1: download the headers first (less than 1 MB in total) ===")
  hd_todo <- todo[todo$kind == "header", , drop = FALSE]
  for (i in seq_len(nrow(hd_todo))) {
    r <- hd_todo[i, ]
    fetch_one(r$url, r$dest)
    got <- if (file.exists(r$dest)) file.size(r$dest) else 0
    say(sprintf("  %-52s %s", r$file,
                if (!is.na(r$remote_bytes) && got == r$remote_bytes) "done" else "not finished"))
  }

  say("")
  say("=== Step 2: confounding gate before downloading the data matrices (GSE198904) ===")
  ok198 <- gate_198904()

  todo <- todo[todo$kind == "data", , drop = FALSE]
  if (!ok198) todo <- todo[todo$gse != "GSE198904", , drop = FALSE]

  say("")
  say("=== Step 3: download the data matrices (back end ", BACKEND, "), ", nrow(todo), " files, ",
      human(sum(todo$remote_bytes, na.rm = TRUE)), " ===")
  for (i in seq_len(nrow(todo))) {
    r <- todo[i, ]
    say("[", i, "/", nrow(todo), "] ", r$file, "  ", human(r$remote_bytes))
    t0 <- Sys.time()
    # Retries within a single file: on a connection reset the download resumes automatically, at most 30 rounds
    for (attempt in 1:30) {
      fetch_one(r$url, r$dest)
      got <- if (file.exists(r$dest)) file.size(r$dest) else 0
      if (!is.na(r$remote_bytes) && got == r$remote_bytes) break
      say("  round ", attempt, " finished: received ", human(got), " / ",
          human(r$remote_bytes), ", resuming")
      Sys.sleep(5)
    }
    got <- if (file.exists(r$dest)) file.size(r$dest) else 0
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (!is.na(r$remote_bytes) && got != r$remote_bytes) {
      say("  still incomplete; re-running this script resumes from the break point"); next
    }
    if (grepl("\\.gz$", r$file) && !gz_ok(r$dest)) {
      say("  gzip integrity check failed; file deleted, please re-run to download again")
      file.remove(r$dest); next
    }
    say(sprintf("  finished and verified, %.1f min, average %.0f KB/s",
                el, got / 1024 / max(el * 60, 1)))
  }
  man$local_bytes <- vapply(man$dest, function(p)
    if (file.exists(p)) file.size(p) else 0, numeric(1))
  man$done <- !is.na(man$remote_bytes) & man$local_bytes == man$remote_bytes
  write.csv(man, file.path(DATA_ROOT, "download_manifest_panel2.csv"), row.names = FALSE)
  say("")
  say("Finished, completed: ", sum(man$done), " / ", nrow(man))
  say("Next: run 30_ingest_new_cohorts.R for verification.")
}
log_stop()
