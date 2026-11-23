# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - cross-disease atlas
# Step 19d: retry failed items + tissue composition + cohort independence
#
# The benchmark of 19c now really holds (GSE87648 returns simplified_diagnosis instead of
# smoking status); tissue decision and group field are credible. Three things remain:
#
#   1   for 201 series the header is still missing, but the error changed from 404 to 503 Service
#       Unavailable - requests were too frequent and were throttled; the files do exist. Retrying with back-off
#       recovers them, and this batch probably adds several tens of whole-blood cohorts.
#
#   2   the tissue flag is currently Boolean: "whole blood if any sample looks like whole blood". Series with mixed
#       tissues slip through (e.g. studies depositing brain tissue with paired peripheral blood).
#       Tissue is now decided per sample and the proportion of whole-blood samples is reported; series with a low proportion are flagged separately.
#
#   3   cohort independence has not been checked. 19c counted 4 Parkinson's disease cohorts, but
#       GSE72774 is known to be a subset of GSE111629; 7 were counted for schizophrenia, of which
#       GSE152026/GSE152027, GSE84727/GSE80417 and GSE116379/GSE116380
#       each very likely share a source. Without excluding double counting, "11 diseases pass the threshold" is inflated.
#       This step uses three lines of evidence for a common source: overlap of array barcodes (the hardest), submitting contact,
#       and similarity of series titles.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/19d_independence.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "19d_independence"

IN_DIR  <- file.path(PROJ_ROOT, "output", "19_inventory")
OUT_DIR <- IN_DIR
HDR_DIR <- file.path(DATA_ROOT, "geo_headers")

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
options(timeout = 600)

RETRY_N     <- 3L
RETRY_SLEEP <- 8       # 503 means throttling; the back-off must be long enough
BARCODE     <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
WB_MIN_FRAC <- 0.6     # Series with a proportion of whole-blood samples below this value are flagged separately

TIS_WHOLE  <- c("whole blood", "peripheral blood", "buffy coat",
                "blood leukocyte", "leukocyte", "^blood$", "\\bblood\\b")
TIS_PBMC   <- c("pbmc", "mononuclear")
TIS_SORTED <- c("cd4\\+", "cd8\\+", "cd14\\+", "cd19\\+", "monocyte",
                "granulocyte", "neutrophil", "t cell", "b cell", "nk cell", "sorted")
TIS_OTHER  <- c("brain", "cortex", "cerebell", "buccal", "saliva", "placenta",
                "cord blood", "adipose", "muscle", "liver", "colon", "skin",
                "biopsy", "tumor", "tumour", "kidney", "lung", "heart", "cell line")

cls1 <- function(v, pats) any(vapply(pats, function(p)
  grepl(p, v, ignore.case = TRUE), logical(1)))

list_matrix_files <- function(gse) {
  stub <- sub("[0-9]{3}$", "nnn", gse)
  u <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", stub, gse)
  html <- tryCatch(paste(readLines(url(u), warn = FALSE), collapse = " "),
                   error = function(e) NULL)
  if (is.null(html)) return(character(0))
  m <- regmatches(html, gregexpr('href="[^"]*series_matrix\\.txt\\.gz"', html))[[1]]
  f <- unique(sub('^href="', "", sub('"$', "", m)))
  if (!length(f)) return(character(0))
  paste0(u, basename(f))
}

fetch_header <- function(gse, retries = RETRY_N) {
  cache <- file.path(HDR_DIR, paste0(gse, "_header.txt"))
  if (file.exists(cache)) return(readLines(cache, warn = FALSE, encoding = "UTF-8"))
  for (k in seq_len(retries)) {
    fl <- list_matrix_files(gse)
    for (u in fl) {
      lines <- tryCatch({
        con <- gzcon(url(u, open = "rb"))
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        acc <- character(0)
        repeat {
          ch <- readLines(con, n = 20L, warn = FALSE)
          if (!length(ch)) break
          acc <- c(acc, ch)
          if (any(grepl("!series_matrix_table_begin", ch, fixed = TRUE))) break
          if (length(acc) >= 400L) break
        }
        acc
      }, error = function(e) NULL)
      if (!is.null(lines) && length(lines)) {
        h <- which(grepl("!series_matrix_table_begin", lines, fixed = TRUE))
        if (length(h)) lines <- lines[seq_len(h[1])]
        writeLines(lines, cache, useBytes = TRUE)
        return(lines)
      }
    }
    Sys.sleep(RETRY_SLEEP * k)      # Back-off
  }
  NULL
}

fld1 <- function(lines, tag) {
  x <- grep(paste0("^", tag), lines, value = TRUE)
  if (!length(x)) return(character(0))
  trimws(gsub('^"|"$', "", strsplit(x[1], "\t", fixed = TRUE)[[1]])[-1])
}
fldall <- function(lines, tag) {
  x <- grep(paste0("^", tag), lines, value = TRUE)
  if (!length(x)) return(list())
  lapply(x, function(one)
    trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1]))
}

# --------------------------------------------------------------
v2 <- read.csv(file.path(IN_DIR, "series_verified_v2.csv"), stringsAsFactors = FALSE)
todo <- v2$accession[v2$tissue == "header not available"]
say("Headers not obtained in the previous round: ", length(todo), "; retrying with back-off in this round")
got <- 0L
for (i in seq_along(todo)) {
  if (!is.null(fetch_header(todo[i]))) got <- got + 1L
  if (i %% 25 == 0) say("  Retried ", i, " / ", length(todo), ", succeeded ", got)
}
say("  Recovered ", got, "")
say("")

# --------------------------------------------------------------
say("=== Re-deciding tissue composition per sample ===")
rows <- list()
for (i in seq_len(nrow(v2))) {
  g <- v2$accession[i]
  f <- file.path(HDR_DIR, paste0(g, "_header.txt"))
  if (!file.exists(f)) {
    rows[[length(rows) + 1]] <- data.frame(
      accession = g, n_wb = NA_integer_, n_tot = NA_integer_,
      wb_frac = NA_real_, tissue2 = "header not available",
      contact = "", barcodes = "", stringsAsFactors = FALSE)
    next
  }
  h <- readLines(f, warn = FALSE, encoding = "UTF-8")
  src <- fld1(h, "!Sample_source_name_ch1")
  ch <- fldall(h, "!Sample_characteristics_ch1")
  tis_line <- NULL
  for (v in ch) if (any(grepl("^tissue\\s*:", v, ignore.case = TRUE))) tis_line <- v
  per <- if (length(src)) src else if (!is.null(tis_line)) tis_line else character(0)
  if (!is.null(tis_line) && length(src) == length(tis_line))
    per <- paste(src, tis_line)
  n_tot <- length(per)
  n_wb <- if (n_tot) sum(vapply(per, function(x)
    cls1(x, TIS_WHOLE) && !cls1(x, TIS_OTHER), logical(1))) else NA_integer_
  frac <- if (is.na(n_wb) || n_tot == 0) NA_real_ else n_wb / n_tot
  t2 <- if (is.na(frac)) "unclear"
  else if (frac >= WB_MIN_FRAC) "whole blood"
  else if (frac > 0) "mixed (contains whole blood)"
  else if (any(vapply(per, cls1, logical(1), TIS_PBMC))) "PBMC"
  else if (any(vapply(per, cls1, logical(1), TIS_SORTED))) "sorted cells"
  else "non-blood"

  ct <- paste(c(fld1(h, "!Series_contact_name"),
                fld1(h, "!Series_contact_institute")), collapse = " / ")
  ttl <- fld1(h, "!Sample_title")
  bc <- unique(unlist(regmatches(ttl, gregexpr(BARCODE, ttl))))
  if (!length(bc)) {
    allch <- unlist(ch)
    bc <- unique(unlist(regmatches(allch, gregexpr(BARCODE, allch))))
  }
  rows[[length(rows) + 1]] <- data.frame(
    accession = g, n_wb = n_wb, n_tot = n_tot, wb_frac = frac, tissue2 = t2,
    contact = ct, barcodes = paste(bc, collapse = ","), stringsAsFactors = FALSE)
}
ex <- do.call(rbind, rows)
m <- merge(v2, ex, by = "accession", all.x = TRUE)
write.csv(m, file.path(OUT_DIR, "series_verified_v3.csv"), row.names = FALSE)

tb <- sort(table(m$tissue2), decreasing = TRUE)
for (i in seq_along(tb)) say(sprintf("  %-18s %d", names(tb)[i], as.integer(tb[i])))
say("  of which mixed (whole blood present but proportion < ", WB_MIN_FRAC, ") can only be used after selecting samples")
say("")

# --------------------------------------------------------------
say("=== Cohort independence ===")
wb <- m[!is.na(m$tissue2) & m$tissue2 %in% c("whole blood", "mixed (contains whole blood)") &
          nzchar(m$group_field) & nzchar(m$disease) &
          !is.na(m$score) & m$score >= 2, , drop = FALSE]
say("Series entering the independence decision: ", nrow(wb))

norm_title <- function(x) {
  x <- tolower(gsub("[^a-z0-9 ]", " ", x))
  gsub("\\s+", " ", trimws(x))
}
jac <- function(a, b) {
  A <- unique(strsplit(a, " ")[[1]]); B <- unique(strsplit(b, " ")[[1]])
  if (!length(A) || !length(B)) return(0)
  length(intersect(A, B)) / length(union(A, B))
}

res <- list()
for (d in unique(unlist(strsplit(wb$disease, "|", fixed = TRUE)))) {
  s <- wb[grepl(d, wb$disease, fixed = TRUE), , drop = FALSE]
  if (nrow(s) < 1) next
  s <- s[order(-s$n_sample), , drop = FALSE]
  n <- nrow(s)
  grp <- seq_len(n)                      # Union-find style merging
  link <- character(0)
  if (n >= 2) for (i in 1:(n - 1)) for (j in (i + 1):n) {
    bi <- strsplit(s$barcodes[i], ",")[[1]]; bj <- strsplit(s$barcodes[j], ",")[[1]]
    bi <- bi[nzchar(bi)]; bj <- bj[nzchar(bj)]
    ov <- if (length(bi) && length(bj)) length(intersect(bi, bj)) else 0L
    same_ct <- nzchar(s$contact[i]) && s$contact[i] == s$contact[j]
    tsim <- jac(norm_title(s$title[i]), norm_title(s$title[j]))
    why <- character(0)
    if (ov > 0) why <- c(why, paste0("barcode overlap ", ov))
    if (same_ct) why <- c(why, "same submitter")
    if (tsim >= 0.6) why <- c(why, sprintf("title similarity %.2f", tsim))
    if (ov > 0 || (same_ct && tsim >= 0.4)) {
      old <- grp[j]; grp[grp == old] <- grp[i]
      link <- c(link, sprintf("%s~%s (%s)", s$accession[i], s$accession[j],
                              paste(why, collapse = ", ")))
    }
  }
  k <- length(unique(grp))
  res[[length(res) + 1]] <- data.frame(
    disease = d, n_series = n, n_independent = k,
    total_samples = sum(s$n_sample), max_series = max(s$n_sample),
    accessions = paste(s$accession, collapse = ", "),
    merged = paste(link, collapse = " | "), stringsAsFactors = FALSE)
}
rr <- do.call(rbind, res)
rr <- rr[order(-rr$n_independent, -rr$total_samples), , drop = FALSE]
for (i in seq_len(nrow(rr))) {
  say(sprintf("  %-28s series %2d -> independent cohorts %2d, samples in total %5d",
              rr$disease[i], rr$n_series[i], rr$n_independent[i], rr$total_samples[i]))
  say("        ", rr$accessions[i])
  if (nzchar(rr$merged[i])) say("        merged: ", rr$merged[i])
}
write.csv(rr, file.path(OUT_DIR, "disease_independence.csv"), row.names = FALSE)

say("")
say("Diseases with >= 2 independent cohorts: ", sum(rr$n_independent >= 2), "")
say("Diseases with >= 3 independent cohorts: ", sum(rr$n_independent >= 3), "")
say("Reminder: barcode overlap is hard evidence; same submitter plus similar titles is only a suspicion, ")
say("      and each pair in the final list still needs a manual look at the original descriptions.")
log_stop()
