# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - epigenetic test of the Brainstorm structure
# Step 19e: scan again with the disease list of the Brainstorm study
#
# Why the expansion is necessary: so far only two psychiatric and two neurological disorders were found. Two disorders
# give a single within-group comparison, and a structural statement such as "psychiatric disorders correlate within the group,
# neurological disorders do not" cannot rest on one pair - each group needs at least three disorders and more than three
# within-group comparisons, otherwise the so-called structure is just a chance difference between two numbers.
#
# The expansion is cheap: the disease dictionary decides what can be seen, and the previous dictionary
# lacked most of the Brainstorm psychiatric group - autism, ADHD, post-traumatic stress disorder,
# anorexia nervosa, obsessive-compulsive disorder and Tourette syndrome were all absent; the neurological group lacked epilepsy,
# stroke and migraine. The metadata of 3180 series and 947 headers are cached locally,
# so a new scan with another dictionary requires no new download.
#
# At the same time the lower limit of the sample number is relaxed from 80 to 60; headers are fetched only for newly matched diseases.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/19e_expand_brainstorm.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "19e_expand_brainstorm"

IN_DIR  <- file.path(PROJ_ROOT, "output", "19_inventory")
OUT_DIR <- IN_DIR
HDR_DIR <- file.path(DATA_ROOT, "geo_headers")

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
options(timeout = 600)

MIN_N <- 60L
RETRY_SLEEP <- 6
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
WB_MIN_FRAC <- 0.6

# Disease groups of Brainstorm. The core of the test is the contrast within and between groups,
# so the group of each disorder must be stated.
GROUP_OF <- c(
  schizophrenia = "psychiatric", "bipolar disorder" = "psychiatric",
  "major depression" = "psychiatric", "anxiety" = "psychiatric",
  PTSD = "psychiatric", ADHD = "psychiatric", autism = "psychiatric",
  "anorexia nervosa" = "psychiatric", OCD = "psychiatric",
  "Tourette" = "psychiatric",
  Alzheimer = "neurological", Parkinson = "neurological",
  "multiple sclerosis" = "neurological", epilepsy = "neurological",
  migraine = "neurological", stroke = "neurological", ALS = "neurological",
  "inflammatory bowel disease" = "comparator",
  "rheumatoid arthritis" = "comparator", "systemic lupus" = "comparator",
  "type 1 diabetes" = "comparator", asthma = "comparator",
  HIV = "comparator", "coronary artery disease" = "comparator")

DISEASE <- list(
  # psychiatric
  schizophrenia        = c("schizophrenia", "schizoaffective", "first.episode psychosis"),
  "bipolar disorder"   = c("bipolar disorder", "bipolar affective"),
  "major depression"   = c("major depress", "depressive disorder", "\\bMDD\\b"),
  "anxiety"            = c("anxiety disorder", "panic disorder"),
  PTSD                 = c("post.?traumatic stress", "\\bPTSD\\b"),
  ADHD                 = c("attention.deficit", "\\bADHD\\b"),
  autism               = c("autism", "autistic", "\\bASD\\b"),
  "anorexia nervosa"   = c("anorexia nervosa"),
  OCD                  = c("obsessive.compulsive"),
  "Tourette"           = c("tourette", "tic disorder"),
  # neurological
  Alzheimer            = c("alzheimer"),
  Parkinson            = c("parkinson"),
  "multiple sclerosis" = c("multiple sclerosis"),
  epilepsy             = c("epilep", "seizure"),
  migraine             = c("migraine"),
  stroke               = c("ischemic stroke", "ischaemic stroke", "cerebral infarct",
                           "intracerebral haemorrhage", "intracerebral hemorrhage"),
  ALS                  = c("amyotrophic lateral"),
  # non-brain control diseases
  "inflammatory bowel disease" = c("inflammatory bowel", "crohn", "ulcerative colitis"),
  "rheumatoid arthritis"       = c("rheumatoid arthritis"),
  "systemic lupus"             = c("lupus erythematosus"),
  "type 1 diabetes"            = c("type 1 diabetes", "type I diabetes"),
  asthma                       = c("asthma"),
  HIV                          = c("\\bHIV\\b"),
  "coronary artery disease"    = c("coronary artery disease", "myocardial infarction")
)

TIS_WHOLE  <- c("whole blood", "peripheral blood", "buffy coat",
                "blood leukocyte", "leukocyte", "^blood$", "\\bblood\\b")
TIS_OTHER  <- c("brain", "cortex", "cerebell", "buccal", "saliva", "placenta",
                "cord blood", "adipose", "muscle", "liver", "colon", "skin",
                "biopsy", "tumor", "tumour", "kidney", "lung", "heart", "cell line")
CTRL_WORDS <- "control|healthy|normal|\\bHC\\b|\\bNC\\b|unaffected|non-|negative"
GOOD_TAG   <- "diagnos|disease|status|group|condition|phenotype|affection|subject|case"
BAD_TAG    <- "smok|sex|gender|age|batch|slide|sentrix|array|position|ethnic|race|bmi|tissue|cell.?type|plate|well|chip|scan|storage|time|visit|treatment|medication|drug"

hit_any <- function(t, p) any(vapply(p, function(x)
  grepl(x, t, ignore.case = TRUE), logical(1)))
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
fetch_header <- function(gse, retries = 3L) {
  cache <- file.path(HDR_DIR, paste0(gse, "_header.txt"))
  if (file.exists(cache)) return(readLines(cache, warn = FALSE, encoding = "UTF-8"))
  for (k in seq_len(retries)) {
    for (u in list_matrix_files(gse)) {
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
    Sys.sleep(RETRY_SLEEP * k)
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
pick_group <- function(lines) {
  ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  best <- list(score = -99, tag = "", levels = "")
  for (one in ch) {
    v <- trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1])
    v <- v[nzchar(v) & grepl(":", v, fixed = TRUE)]
    if (!length(v)) next
    tag <- trimws(sub(":.*$", "", v[1]))
    val <- trimws(sub("^[^:]*:", "", v))
    tb <- sort(table(val), decreasing = TRUE)
    if (length(tb) < 2L || length(tb) > 8L) next
    sc <- 0
    if (grepl(GOOD_TAG, tag, ignore.case = TRUE)) sc <- sc + 1
    if (grepl(BAD_TAG,  tag, ignore.case = TRUE)) sc <- sc - 5
    if (any(grepl(CTRL_WORDS, names(tb), ignore.case = TRUE))) sc <- sc + 2
    if (length(tb) <= 3L) sc <- sc + 1
    if (min(tb) >= 20L) sc <- sc + 1
    if (sc > best$score) best <- list(score = sc, tag = tag,
      levels = paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; "))
  }
  best
}

# --------------------------------------------------------------
all <- read.csv(file.path(IN_DIR, "geo_all_series.csv"), stringsAsFactors = FALSE)
say("Cached series metadata: ", nrow(all), " records")
pool <- all[!is.na(all$n_sample) & all$n_sample >= MIN_N, , drop = FALSE]
say("sample number >= ", MIN_N, ": ", nrow(pool), " (the previous round used 80)")

txt0 <- paste(pool$title, pool$summary)
pool$disease0 <- vapply(txt0, function(t) {
  h <- names(DISEASE)[vapply(DISEASE, function(p) hit_any(t, p), logical(1))]
  if (!length(h)) "" else paste(h, collapse = "|")
}, character(1))
need <- pool[nzchar(pool$disease0), , drop = FALSE]
say("Matching a disease term under the new dictionary: ", nrow(need), "")
have <- file.exists(file.path(HDR_DIR, paste0(need$accession, "_header.txt")))
say("of which with cached header: ", sum(have), ", to fetch: ", sum(!have))
say("")

if (any(!have)) {
  tgt <- need$accession[!have]
  ok <- 0L
  for (i in seq_along(tgt)) {
    if (!is.null(fetch_header(tgt[i]))) ok <- ok + 1L
    if (i %% 25 == 0) say("  Fetched ", i, " / ", length(tgt), ", succeeded ", ok)
  }
  say("  Recovered ", ok, "")
  say("")
}

say("=== Deciding one by one ===")
rows <- list()
for (i in seq_len(nrow(need))) {
  g <- need$accession[i]
  f <- file.path(HDR_DIR, paste0(g, "_header.txt"))
  if (!file.exists(f)) next
  h <- readLines(f, warn = FALSE, encoding = "UTF-8")
  src <- fld1(h, "!Sample_source_name_ch1")
  ch <- fldall(h, "!Sample_characteristics_ch1")
  tis_line <- NULL
  for (v in ch) if (any(grepl("^tissue\\s*:", v, ignore.case = TRUE))) tis_line <- v
  per <- if (length(src)) src else if (!is.null(tis_line)) tis_line else character(0)
  if (!is.null(tis_line) && length(src) == length(tis_line)) per <- paste(src, tis_line)
  n_tot <- length(per)
  n_wb <- if (n_tot) sum(vapply(per, function(x)
    cls1(x, TIS_WHOLE) && !cls1(x, TIS_OTHER), logical(1))) else 0L
  frac <- if (n_tot == 0) NA_real_ else n_wb / n_tot
  gp <- pick_group(h)
  ttl <- fld1(h, "!Sample_title")
  bc <- unique(unlist(regmatches(ttl, gregexpr(BARCODE, ttl))))
  ct <- paste(c(fld1(h, "!Series_contact_name"),
                fld1(h, "!Series_contact_institute")), collapse = " / ")
  txt <- paste(need$title[i], need$summary[i], paste(utils::head(per, 40), collapse = " "))
  dz <- names(DISEASE)[vapply(DISEASE, function(p) hit_any(txt, p), logical(1))]
  rows[[length(rows) + 1]] <- data.frame(
    accession = g, n_sample = need$n_sample[i], n_wb = n_wb, wb_frac = frac,
    group_field = gp$tag, score = gp$score, levels = gp$levels,
    disease = paste(dz, collapse = "|"), contact = ct,
    barcodes = paste(bc, collapse = ","), title = need$title[i],
    stringsAsFactors = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, "brainstorm_candidates.csv"), row.names = FALSE)
say("Decided: ", nrow(res), "")

good <- res[!is.na(res$wb_frac) & res$wb_frac >= WB_MIN_FRAC &
              nzchar(res$group_field) & res$score >= 2 & res$n_wb >= MIN_N,
            , drop = FALSE]
say("Passing the tissue and group filters: ", nrow(good), "")
say("")

say("=== Summary by Brainstorm group ===")
out <- list()
for (d in names(DISEASE)) {
  s <- good[grepl(d, good$disease, fixed = TRUE), , drop = FALSE]
  if (!nrow(s)) next
  s <- s[order(-s$n_wb), , drop = FALSE]
  out[[length(out) + 1]] <- data.frame(
    disease = d, brainstorm_group = unname(GROUP_OF[d]),
    n_series = nrow(s), total_wb = sum(s$n_wb), max_wb = max(s$n_wb),
    accessions = paste(s$accession, collapse = ", "),
    stringsAsFactors = FALSE)
}
tt <- do.call(rbind, out)
tt <- tt[order(tt$brainstorm_group, -tt$n_series, -tt$total_wb), , drop = FALSE]
for (grp in c("psychiatric", "neurological", "comparator")) {
  say("---- ", grp, " ----")
  s <- tt[!is.na(tt$brainstorm_group) & tt$brainstorm_group == grp, , drop = FALSE]
  if (!nrow(s)) { say("   (none)"); next }
  for (i in seq_len(nrow(s))) {
    say(sprintf("   %-26s series %2d, whole-blood samples in total %5d, largest set %4d",
                s$disease[i], s$n_series[i], s$total_wb[i], s$max_wb[i]))
    say("        ", s$accessions[i])
  }
  say("   Diseases of this group with >= 2 series: ", sum(s$n_series >= 2))
}
write.csv(tt, file.path(OUT_DIR, "brainstorm_summary.csv"), row.names = FALSE)

say("")
say("Reminder: these are still candidates. The values of the group fields must be inspected one by one and independence verified pair by pair ")
say("      before a series enters the final list - 19d already misjudged the PD pair, for which the answer is known.")
log_stop()
