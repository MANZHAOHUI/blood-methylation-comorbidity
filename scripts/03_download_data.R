# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 03: download of the public data
#
# Part 0: compare Sentrix barcodes between GSE145361 and GSE111629 (zero cost, reads the local cache)
# Part 1: list the files to download and their remote sizes by tier, and build a manifest
# Part 2: download with resume, checking size and integrity file by file
#
# Download back end: aria2c is preferred (multi-connection segmented download, usually an order of magnitude faster than a single connection on international links);
#           if it is not found, fall back to the curl.exe shipped with Windows. Both support resuming,
#           so after an interruption re-running this script continues the download.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/03_download_data.R
# ==============================================================

DATA_ROOT <- "C:/blood_methylation_download"
HDR_DIR   <- "C:/blood_methylation_project/output/01_metadata"

# Download up to which tier:
#   1 = enough to run the whole pipeline (IBD primary analysis set + PD confounding adjudication set)
#   2 = the remaining cohorts needed for the primary analysis
#   3 = plus the large GSE145361 matrix (large, not essential)
#   4 = plus all RAW.tar files (IDAT, about 32 GB; no longer needed by this design)
MAX_TIER    <- 3L
DO_DOWNLOAD <- TRUE   # Run with FALSE first to see the manifest; switch to TRUE after checking it

# Data sets already downloaded on the server are listed here and are not downloaded again locally
SKIP_GSE <- c("GSE87648")

if (!dir.exists(DATA_ROOT)) dir.create(DATA_ROOT, recursive = TRUE)
if (!dir.exists(DATA_ROOT)) stop("Failed to create the data directory; please confirm that drive D is writable: ", DATA_ROOT)

ARIA <- Sys.which("aria2c")
CURL <- Sys.which("curl")
if (!nzchar(ARIA) && !nzchar(CURL)) stop("Neither aria2c nor curl.exe was found")
BACKEND <- if (nzchar(ARIA)) "aria2c" else "curl"
message("Download back end: ", BACKEND, "  (", if (BACKEND == "aria2c") ARIA else CURL, ")")

human <- function(b) {
  if (is.na(b)) return("unknown")
  u <- c("B", "KB", "MB", "GB", "TB"); i <- 1
  while (b >= 1024 && i < length(u)) { b <- b / 1024; i <- i + 1 }
  sprintf("%.1f %s", b, u[i])
}

# ==============================================================
# Part 0: Sentrix barcode comparison
# ==============================================================
sentrix_check <- function() {
  f1 <- file.path(HDR_DIR, "GSE145361_header.txt")
  f2 <- file.path(HDR_DIR, "GSE111629_header.txt")
  if (!file.exists(f1) || !file.exists(f2)) {
    message("\n[Part 0] cached headers not found; barcode comparison skipped"); return(invisible(NULL))
  }
  get_titles <- function(f) {
    l <- readLines(f, warn = FALSE, encoding = "UTF-8")
    tl <- grep("^!Sample_title", l, value = TRUE)
    if (!length(tl)) return(character(0))
    trimws(gsub('^"|"$', "", strsplit(tl[1], "\t", fixed = TRUE)[[1]])[-1])
  }
  # GSE111629 titles carry a fixed prefix that must be stripped before comparing barcodes
  norm <- function(x) {
    x <- sub("^genomic DNA from\\s+", "", x)
    x <- toupper(trimws(x))
    x[grepl("^[0-9]{8,}_R[0-9]{2}C[0-9]{2}$", x)]
  }
  a <- norm(get_titles(f1)); b <- norm(get_titles(f2))
  inter <- intersect(a, b)
  message("\n[Part 0] Sentrix barcode comparison")
  message("  GSE145361 barcodes parsed: ", length(a), "")
  message("  GSE111629 barcodes parsed: ", length(b), "")
  message("  barcodes shared by both: ", length(inter), "")
  if (length(inter)) message("  e.g. ", paste(head(inter, 8), collapse = ", "))
  writeLines(c("GSE145361 vs GSE111629 Sentrix barcode overlap",
               paste0("GSE145361 parsed: ", length(a)),
               paste0("GSE111629 parsed: ", length(b)),
               paste0("shared: ", length(inter)), inter),
             file.path(DATA_ROOT, "00_sentrix_overlap.txt"), useBytes = TRUE)
  invisible(inter)
}
sentrix_check()

# ==============================================================
# Files to download (by tier)
# ==============================================================
# Source other than GEO: public summary statistics and sample QC information of Vallerga 2020
# (PCTG site, University of Queensland). The ancestry file is small but important: of the
# 1889 individuals in GSE145361 only 1638 are unrelated, and this file identifies those 1638.
EXTRA <- data.frame(rbind(
  c(1, "VALLERGA", "Vallerga2020_NCOMMS_ancestry.txt"),
  c(1, "VALLERGA", "Vallerga2020_NCOMMS_MWAS-meta-analysis-MOA.mlma"),
  c(2, "VALLERGA", "Vallerga2020_NCOMMS_MWAS-meta-analysis-MOMENT.mlma")
), stringsAsFactors = FALSE)
names(EXTRA) <- c("tier", "gse", "file")
EXTRA$tier <- as.integer(EXTRA$tier)
EXTRA$url <- paste0("https://cnsgenomics.com/data/vallerga_et_al_2020_nc/", EXTRA$file)

FILES <- data.frame(rbind(
  c(1, "GSE87648",  "GSE87648_processedMethWb.txt.gz"),
  c(1, "GSE72774",  "GSE72774_datBetaNormalized.csv.gz"),
  c(1, "GSE72774",  "GSE72774_Explanation_of_characteristic_variables.pdf"),
  c(2, "GSE111629", "GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz"),
  c(2, "GSE112611", "GSE112611_beta_values.txt.gz"),
  c(2, "GSE152027", "GSE152027_IOP_processed_signals.csv.gz"),
  c(3, "GSE145361", "GSE145361_Vallerga2020_NCOMMS_AvgBeta_Matrix-file.txt.gz"),
  c(4, "GSE87648",  "GSE87648_RAW.tar"),
  c(4, "GSE72774",  "GSE72774_RAW.tar"),
  c(4, "GSE111629", "GSE111629_RAW.tar"),
  c(4, "GSE112611", "GSE112611_RAW.tar"),
  c(4, "GSE152027", "GSE152027_RAW.tar"),
  c(4, "GSE145361", "GSE145361_RAW.tar")
), stringsAsFactors = FALSE)
names(FILES) <- c("tier", "gse", "file")
FILES$tier <- as.integer(FILES$tier)
FILES$url <- NA_character_          # NA means the URL is built according to the GEO path rules
FILES <- rbind(FILES, EXTRA)
FILES <- FILES[FILES$tier <= MAX_TIER, , drop = FALSE]
if (length(SKIP_GSE)) FILES <- FILES[!(FILES$gse %in% SKIP_GSE), , drop = FALSE]

supp_url <- function(gse, fname) {
  stub <- sub("[0-9]{3}$", "nnn", gse)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s", stub, gse, fname)
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
# Part 1: manifest
# ==============================================================
message("\n[Part 1] querying remote sizes (tier <= ", MAX_TIER, ")")
rows <- list()
for (i in seq_len(nrow(FILES))) {
  gse <- FILES$gse[i]; f <- FILES$file[i]
  d <- file.path(DATA_ROOT, gse)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
  u <- if (is.na(FILES$url[i])) supp_url(gse, f) else FILES$url[i]
  sz <- remote_size(u)
  dest <- file.path(d, f)
  local <- if (file.exists(dest)) file.size(dest) else 0
  message(sprintf("  T%d  %-12s %-62s %10s", FILES$tier[i], gse, f, human(sz)))
  rows[[length(rows) + 1]] <- data.frame(
    tier = FILES$tier[i], gse = gse, file = f, url = u, dest = dest,
    remote_bytes = sz, local_bytes = local,
    done = !is.na(sz) && local == sz, stringsAsFactors = FALSE)
}
man <- do.call(rbind, rows)
write.csv(man, file.path(DATA_ROOT, "download_manifest.csv"), row.names = FALSE)
todo <- man[!man$done, ]
message("\n  Files in total: ", nrow(man), ", to download: ", nrow(todo),
        ", total size: ", human(sum(todo$remote_bytes, na.rm = TRUE)))

# ==============================================================
# Part 2: download
# ==============================================================
# The whole stream must be decompressed before a file passes. Reading only the first bytes cannot detect holes in the middle:
# with parallel segments aria2c first extends the file to full length, unwritten segments are holes,
# the size matches and the beginning decompresses, and it only fails when that segment is read.
gz_ok <- function(path) {
  tryCatch({
    con <- gzfile(path, "rb"); on.exit(close(con), add = TRUE)
    repeat {
      b <- readBin(con, "raw", n = 8L * 1024L^2)
      if (length(b) == 0L) break
    }
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
}

# aria2c removes the .aria2 control file only on a clean finish; if it is still there the download is incomplete
aria_incomplete <- function(dest) file.exists(paste0(dest, ".aria2"))

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
                    "--speed-limit", "1000", "--speed-time", "60",
                    "--connect-timeout", "60",
                    "-o", shQuote(dest), shQuote(u)))
  }
}

if (!DO_DOWNLOAD) {
  message("\nDO_DOWNLOAD = FALSE: manifest only, nothing downloaded. Switch to TRUE and re-run after checking.")
} else {
  log_path <- file.path(DATA_ROOT, "download_log.txt")
  logline <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
    message(msg); cat(msg, "\n", file = log_path, append = TRUE)
  }
  logline("Starting download (back end ", BACKEND, "), to download: ", nrow(todo), "")

  for (i in seq_len(nrow(todo))) {
    r <- todo[i, ]
    logline("[", i, "/", nrow(todo), "] ", r$file, "  ", human(r$remote_bytes))

    # Retries within a single file. No fixed limit on the number of rounds: connection resets are normal, and as long as
    # each round brings in new bytes the download continues. Only after STALL_LIMIT consecutive rounds without a single byte
    # is the file considered stuck rather than slow; the script then stops and leaves the partial file, and a re-run resumes.
    STALL_LIMIT <- 12L
    stall <- 0L; attempt <- 0L
    repeat {
      attempt <- attempt + 1L
      before <- if (file.exists(r$dest)) file.size(r$dest) else 0
      st <- fetch_one(r$url, r$dest)
      got <- if (file.exists(r$dest)) file.size(r$dest) else 0
      done_now <- !is.na(r$remote_bytes) && got == r$remote_bytes &&
        !(BACKEND == "aria2c" && aria_incomplete(r$dest))
      if (done_now) break
      stall <- if (got > before) 0L else stall + 1L
      logline("  round ", attempt, ": received ", human(got), " / ",
              human(r$remote_bytes),
              " (new in this round ", human(got - before), ", idle rounds ", stall, ")")
      if (stall >= STALL_LIMIT) {
        logline("  no progress in ", STALL_LIMIT, " consecutive rounds; this file is paused and a re-run continues it")
        break
      }
      Sys.sleep(5)
    }

    got <- if (file.exists(r$dest)) file.size(r$dest) else 0
    if (!is.na(r$remote_bytes) && got != r$remote_bytes) {
      logline("  still incomplete; re-running this script resumes from the break point"); next
    }
    if (BACKEND == "aria2c" && aria_incomplete(r$dest)) {
      logline("  size reached but the .aria2 control file is still present, so segments are incomplete; re-run this script to resume")
      next
    }
    if (grepl("\\.gz$", r$file) && !gz_ok(r$dest)) {
      logline("  gzip integrity check failed; file deleted, please re-run to download again")
      file.remove(r$dest); next
    }
    logline("  finished and verified")
  }

  man$local_bytes <- vapply(man$dest, function(p)
    if (file.exists(p)) file.size(p) else 0, numeric(1))
  man$done <- !is.na(man$remote_bytes) & man$local_bytes == man$remote_bytes
  write.csv(man, file.path(DATA_ROOT, "download_manifest.csv"), row.names = FALSE)
  logline("Finished, completed: ", sum(man$done), " / ", nrow(man))
}

