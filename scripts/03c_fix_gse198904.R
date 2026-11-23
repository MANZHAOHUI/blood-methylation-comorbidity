# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 03c: splitting GSE198904 and downloading the usable part
#
# The gate in 03b stopped this data set, rightly, but the action it prescribed was too strict:
#
#   Cohort 1: HC 32, MDD 0      - controls only; another 191 samples have the diagnosis recorded as OBSERVEMDD0001
#   Cohort 2: HC 68, MDD 361    - both cases and controls
#
# So the data set is not unusable as a whole: Cohort 2 alone is a valid case-control set
# (361 vs 68, both sides above our minimum). The rule "skip the whole set if any cell is zero"
# was hard-coded; the correct action is to drop the problematic sub-study and keep the usable one.
# The action attached to the criterion was wrong, not the criterion itself - it did catch a real problem.
#
# This step: determine from the header which supplementary matrix corresponds to Cohort 2 and download only that one.
# The two candidates are DHRC (1.4 GB) and OBSERVEMDD0001 (1.3 GB); a wrong guess costs
# one or two hours of download, so the correspondence between sample titles and sub-studies is checked first.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/03c_fix_gse198904.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "C:/blood_methylation_download"
SCRIPT_NAME <- "03c_fix_gse198904"

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

DO_DOWNLOAD <- FALSE      # Run with FALSE first to see the decision; switch to TRUE after checking
GSE <- "GSE198904"
RAW_DIR <- file.path(DATA_ROOT, "raw")
DIR <- file.path(RAW_DIR, GSE)
HDR <- file.path(DIR, paste0(GSE, "_series_matrix.txt.gz"))
if (!file.exists(HDR)) stop("Header file not found; please run 03b first: ", HDR)

CAND <- c(DHRC = "GSE198904_DHRC_Matrix_processed.txt.gz",
          OBSERVEMDD0001 = "GSE198904_OBSERVEMDD0001_Matrix_processed.txt.gz")

ARIA <- Sys.which("aria2c"); CURL <- Sys.which("curl")
BACKEND <- if (nzchar(ARIA)) "aria2c" else "curl"

con <- gzfile(HDR, "rt"); hd <- readLines(con, n = 400L, warn = FALSE); close(con)
sl <- function(x) trimws(gsub('^"|"$', "", strsplit(x, "\t", fixed = TRUE)[[1]])[-1])
gsm <- sl(grep("^!Sample_geo_accession", hd, value = TRUE)[1])
ttl <- sl(grep("^!Sample_title", hd, value = TRUE)[1])
pull <- function(tag) {
  for (one in grep("^!Sample_characteristics_ch1", hd, value = TRUE)) {
    v <- sl(one); vv <- v[nzchar(v) & grepl(":", v, fixed = TRUE)]
    if (!length(vv)) next
    if (tolower(trimws(sub(":.*$", "", vv[1]))) == tolower(tag))
      return(trimws(sub("^[^:]*:", "", v)))
  }
  NULL
}
dg <- pull("diagnosis"); ch <- pull("cohort")
say("Samples: ", length(gsm), "; diagnosis and sub-study fields both found: ", !is.null(dg) && !is.null(ch))

say("")
say("=== Sub-study x diagnosis ===")
ct <- table(ch, dg)
for (r in rownames(ct))
  say("  ", r, ": ", paste(sprintf("%s=%d", colnames(ct), as.integer(ct[r, ])),
                           collapse = "; "))

say("")
say("=== Sample title examples of each sub-study (to decide which matrix file it corresponds to) ===")
for (cc in unique(ch)) {
  ii <- which(ch == cc)
  say("  ", cc, " (", length(ii), " samples): ",
      paste(utils::head(ttl[ii], 4), collapse = " | "))
}

# Decision: usable sub-study = the one in which neither cases nor controls are zero
usable <- rownames(ct)[apply(ct[, intersect(colnames(ct), c("HC", "MDD")), drop = FALSE],
                             1, function(z) all(z > 0))]
say("")
if (!length(usable)) {
  say("[FAIL] no sub-study contains both cases and controls; this data set is unusable")
  log_stop()
} else {
say("Usable sub-study: ", paste(usable, collapse = ", "))
keep <- ch %in% usable & dg %in% c("HC", "MDD")
say("  Included: ", sum(keep), " samples: controls ", sum(keep & dg == "HC"),
    ", cases ", sum(keep & dg == "MDD"))

# The sample titles are array barcodes and the file names carry no study name, so matching by title does not work.
# A more direct way is used instead: fetch only the first 2 MB of each candidate file
# with an HTTP range request. gzip is a stream format, so the first 2 MB are enough to decompress the header line,
# and the header holds the sample column names. The decision then costs a few MB rather than 1.4 GB.
peek_header <- function(url) {
  tmp <- tempfile(fileext = ".gz")
  on.exit(unlink(tmp), add = TRUE)
  st <- system2(CURL, c("-sL", "-r", "0-2097151", "--retry", "3",
                        "-o", shQuote(tmp), shQuote(url)))
  if (st != 0L || !file.exists(tmp)) return(character(0))
  ln <- tryCatch(suppressWarnings({
    cn <- gzcon(file(tmp, "rb")); on.exit(try(close(cn), silent = TRUE), add = TRUE)
    readLines(cn, n = 1L, warn = FALSE)
  }), error = function(e) character(0))
  if (!length(ln)) return(character(0))
  x <- trimws(gsub('^"|"$', "", strsplit(ln[1], "[\t,]")[[1]]))
  x[nzchar(x)]
}

say("")
say("=== Fetching the headers of the two candidate files (2 MB each) ===")
bc2 <- ttl[keep]                       # Barcodes of the usable sub-study
bc1 <- ttl[!keep]                      # Barcodes of the other sub-study
gs2 <- gsm[keep]; gs1 <- gsm[!keep]
score <- setNames(numeric(length(CAND)), names(CAND))
for (k in names(CAND)) {
  u <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/%s/suppl/%s", GSE, CAND[[k]])
  cols <- peek_header(u)
  if (!length(cols)) { say("  ", k, ": failed to read the header"); next }
  hit2 <- sum(cols %in% bc2) + sum(cols %in% gs2)
  hit1 <- sum(cols %in% bc1) + sum(cols %in% gs1)
  say(sprintf("  %-16s columns %5d; %4d match Cohort 2, %4d match Cohort 1",
              k, length(cols), hit2, hit1))
  say("    column name examples: ", paste(utils::head(cols, 4), collapse = " | "))
  score[[k]] <- hit2 - hit1
}
pick <- if (max(score) > 0) names(score)[which.max(score)] else NA_character_
say("")
if (is.na(pick)) {
  say("  [FAIL] the column names of neither file match the sample IDs or barcodes; the files need to be inspected manually")
} else {
  say("  Decision: Cohort 2 corresponds to ", CAND[[pick]],
      " (lead in matches: ", max(score), ")")
}

write.csv(data.frame(gsm = gsm, title = ttl, cohort = ch, diagnosis = dg,
                     keep = keep, stringsAsFactors = FALSE),
          file.path(DATA_ROOT, "GSE198904_cohort2_samples.csv"), row.names = FALSE)
say("")
say("Sample list written, to be used by step 30 via the keep column")

if (DO_DOWNLOAD && !is.na(pick)) {
  f <- CAND[[pick]]
  u <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/%s/suppl/%s", GSE, f)
  dest <- file.path(DIR, f)
  say("")
  say("=== Downloading ", f, " ===")
  for (attempt in 1:30) {
    if (BACKEND == "aria2c") {
      system2(ARIA, c("-c", "-x", "16", "-s", "16", "-k", "20M",
                      "--max-connection-per-server=16", "--max-tries=0",
                      "--retry-wait=10", "--connect-timeout=60", "--timeout=120",
                      "--file-allocation=none", "--summary-interval=30",
                      "-d", shQuote(DIR), "-o", shQuote(f), shQuote(u)))
    } else {
      system2(CURL, c("-L", "-C", "-", "--retry", "50", "--retry-delay", "10",
                      "-o", shQuote(dest), shQuote(u)))
    }
    if (file.exists(dest) && file.size(dest) > 1e9) break
    say("  round ", attempt, ": received ",
        sprintf("%.2f GB", (if (file.exists(dest)) file.size(dest) else 0) / 1024^3),
        ", resuming")
    Sys.sleep(5)
  }
  ok <- tryCatch({ cn <- gzfile(dest, "rb"); on.exit(close(cn), add = TRUE)
                   length(readBin(cn, "raw", n = 1024L)) > 0 },
                 error = function(e) FALSE)
  say("  ", if (ok) "finished and passed the gzip check" else "[FAIL] gzip check failed; please re-run")
} else if (!DO_DOWNLOAD) {
  say("")
  say("DO_DOWNLOAD = FALSE: nothing downloaded. Switch to TRUE after checking the decision above.")
}
log_stop()
}   # end of the branch "a usable sub-study exists"
