# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - cross-disease atlas
# Step 19b: stage 2 of the inventory, now reading the real metadata
#
# Why redone: the filter of 19 did not pass the self-check. Six data sets are known to qualify,
# and any credible filter must recover all of them; 19 did not -
# the largest set in the Parkinson's disease row has only 508 samples although GSE145361 has 1889;
# the largest in the schizophrenia row has 95 although GSE152027 has 800 and GSE84727
# has 847. Things known to be true were filtered out, so the output cannot be used.
#
# Both causes lie in the keyword logic of 19: the exclusion list is too broad (longitudinal,
# brain, biopsy and twin wrongly remove many proper case-control studies; 1103 of 3180
# series were rejected by it), and tissue was judged from title and summary only, whereas
# many series do not write "whole blood" in the summary at all but in the sample fields.
#
# This step reads the series matrix header of every series and takes the sample-level
# tissue and group fields directly - as was done for the first six data sets. Keywords
# are only used to guess the disease, no longer for inclusion decisions.
#
# The script carries a benchmark test: the six known data sets must all be recovered and classified correctly;
# a failure is reported explicitly in the log rather than silently delivering a list.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/19b_inventory_verify.R
# Note: headers are cached on disk; after an interruption a re-run only fetches the missing ones and repeats no request
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "19b_inventory_verify"

IN_DIR  <- file.path(PROJ_ROOT, "output", "19_inventory")
OUT_DIR <- IN_DIR
HDR_DIR <- file.path(DATA_ROOT, "geo_headers")
for (d in c(OUT_DIR, HDR_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

options(timeout = 600)

MIN_N     <- 80L     # Lower limit of the total sample number, the only pre-filter
MAX_FETCH <- 2000L   # Maximum number of headers fetched in this run; several runs complete the job

# Benchmark: six data sets known to qualify, all of which must be recovered
BENCH <- c("GSE145361", "GSE111629", "GSE72774",
           "GSE87648", "GSE112611", "GSE152027")

# Tissue is decided from sample fields, not from the summary
TIS_WHOLE  <- c("whole blood", "peripheral blood", "buffy coat",
                "blood leukocyte", "leukocytes", "peripheral whole blood")
TIS_PBMC   <- c("pbmc", "mononuclear")
TIS_SORTED <- c("cd4", "cd8", "cd14", "cd19", "monocyte", "granulocyte",
                "neutrophil", "t cell", "b cell", "nk cell", "sorted")
TIS_OTHER  <- c("brain", "cortex", "buccal", "saliva", "placenta", "cord",
                "adipose", "muscle", "liver", "colon", "skin", "biopsy",
                "tumor", "tumour", "kidney", "lung", "heart", "cell line")

GROUPISH <- "diagnos|disease|status|group|condition|phenotype|case|control|subject type|sample type|affection"

# Disease dictionary: several dangerous abbreviations of 19 were removed (MS often means mass spectrometry, RA often something else,
# TB is too short); only self-explanatory long terms are used
DISEASE <- list(
  "inflammatory bowel disease" = c("inflammatory bowel", "crohn", "ulcerative colitis"),
  "rheumatoid arthritis"       = c("rheumatoid arthritis"),
  "systemic lupus"             = c("lupus erythematosus", "\\bSLE\\b"),
  "psoriasis"                  = c("psoriasis", "psoriatic"),
  "multiple sclerosis"         = c("multiple sclerosis"),
  "type 1 diabetes"            = c("type 1 diabetes", "type I diabetes"),
  "type 2 diabetes"            = c("type 2 diabetes", "type II diabetes"),
  "obesity"                    = c("obesity", "body mass index"),
  "asthma"                     = c("asthma"),
  "atopic dermatitis"          = c("atopic dermatitis", "eczema"),
  "celiac disease"             = c("celiac", "coeliac"),
  "Sjogren"                    = c("sjogren", "sj\u00f6gren"),
  "systemic sclerosis"         = c("systemic sclerosis", "scleroderma"),
  "ankylosing spondylitis"     = c("ankylosing spondylitis"),
  "Graves disease"             = c("graves' disease", "graves disease"),
  "schizophrenia"              = c("schizophrenia", "psychosis", "psychotic"),
  "bipolar disorder"           = c("bipolar"),
  "major depression"           = c("depressive", "depression"),
  "Alzheimer"                  = c("alzheimer", "dementia"),
  "Parkinson"                  = c("parkinson"),
  "ALS"                        = c("amyotrophic lateral"),
  "coronary artery disease"    = c("coronary", "myocardial infarction", "atheroscler"),
  "heart failure"              = c("heart failure", "cardiomyopathy"),
  "chronic kidney disease"     = c("chronic kidney", "nephropathy"),
  "COPD"                       = c("chronic obstructive"),
  "NAFLD"                      = c("fatty liver", "NAFLD", "NASH", "MASLD"),
  "preeclampsia"               = c("preeclampsia", "pre-eclampsia"),
  "HIV"                        = c("\\bHIV\\b"),
  "tuberculosis"               = c("tuberculosis"),
  "COVID"                      = c("COVID", "SARS-CoV-2"),
  "rheumatic/other autoimmune" = c("autoimmun")
)
hit_any <- function(t, p) any(vapply(p, function(x)
  grepl(x, t, ignore.case = TRUE), logical(1)))
cls <- function(v, pats) any(vapply(pats, function(p)
  any(grepl(p, v, ignore.case = TRUE)), logical(1)))

# --------------------------------------------------------------
fetch_header <- function(gse) {
  cache <- file.path(HDR_DIR, paste0(gse, "_header.txt"))
  if (file.exists(cache)) return(readLines(cache, warn = FALSE, encoding = "UTF-8"))
  stub <- sub("[0-9]{3}$", "nnn", gse)
  base <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/", stub, gse)
  sufs <- c("_series_matrix.txt.gz", "-GPL13534_series_matrix.txt.gz",
            "-GPL21145_series_matrix.txt.gz", "-GPL33022_series_matrix.txt.gz")
  for (u in paste0(base, gse, sufs)) {
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
  NULL
}

fld <- function(lines, tag) {
  x <- grep(paste0("^", tag), lines, value = TRUE)
  if (!length(x)) return(character(0))
  unlist(lapply(x, function(one)
    trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1])))
}

# --------------------------------------------------------------
f_all <- file.path(IN_DIR, "geo_all_series.csv")
if (!file.exists(f_all)) stop("Output of 19 missing: ", f_all)
all <- read.csv(f_all, stringsAsFactors = FALSE)
say("Series retrieved by 19: ", nrow(all))

pool <- all[!is.na(all$n_sample) & all$n_sample >= MIN_N, , drop = FALSE]
pool <- pool[order(-pool$n_sample), , drop = FALSE]
say("sample number >= ", MIN_N, ": ", nrow(pool), " (the only pre-filter of this step)")

miss_bench <- setdiff(BENCH, all$accession)
if (length(miss_bench)) {
  say("[BENCHMARK WARNING] the following known data sets are not even in the complete table of 19: ",
      paste(miss_bench, collapse = ", "))
}
pool <- rbind(pool[pool$accession %in% BENCH, , drop = FALSE],
              pool[!pool$accession %in% BENCH, , drop = FALSE])

n_do <- min(nrow(pool), MAX_FETCH)
say("Processing the first ", n_do, "")
say("")

rows <- list()
t0 <- Sys.time()
for (i in seq_len(n_do)) {
  g <- pool$accession[i]
  h <- fetch_header(g)
  if (is.null(h)) {
    rows[[length(rows) + 1]] <- data.frame(
      accession = g, n_sample = pool$n_sample[i], ok = FALSE,
      tissue = "header not available", group_field = "", levels = "",
      disease = "", title = pool$title[i], stringsAsFactors = FALSE)
    next
  }
  src <- c(fld(h, "!Sample_source_name_ch1"), fld(h, "!Sample_characteristics_ch1"))
  src <- src[nzchar(src)]
  tis <- if (cls(src, TIS_WHOLE)) "whole blood"
  else if (cls(src, TIS_PBMC)) "PBMC"
  else if (cls(src, TIS_SORTED)) "sorted cells"
  else if (cls(src, TIS_OTHER)) "non-blood"
  else "unclear"

  ch <- grep("^!Sample_characteristics_ch1", h, value = TRUE)
  gf <- ""; lv <- ""
  for (one in ch) {
    v <- trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1])
    v <- v[nzchar(v) & grepl(":", v, fixed = TRUE)]
    if (!length(v)) next
    tag <- trimws(sub(":.*$", "", v[1]))
    if (!grepl(GROUPISH, tag, ignore.case = TRUE)) next
    val <- trimws(sub("^[^:]*:", "", v))
    tb <- sort(table(val), decreasing = TRUE)
    if (length(tb) < 2L || length(tb) > 8L) next
    gf <- tag
    lv <- paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
    break
  }
  txt <- paste(pool$title[i], pool$summary[i], paste(utils::head(src, 30), collapse = " "))
  dz <- names(DISEASE)[vapply(DISEASE, function(p) hit_any(txt, p), logical(1))]

  rows[[length(rows) + 1]] <- data.frame(
    accession = g, n_sample = pool$n_sample[i], ok = TRUE,
    tissue = tis, group_field = gf, levels = lv,
    disease = paste(dz, collapse = "|"),
    title = pool$title[i], stringsAsFactors = FALSE)

  if (i %% 50 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("  Processed ", i, " / ", n_do, ", time ", sprintf("%.1f", el),
        " min, estimated remaining ", sprintf("%.0f", el / i * (n_do - i)), " min")
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, "series_verified.csv"), row.names = FALSE)

# --------------------------------------------------------------
say("")
say("=== Benchmark test: six known data sets ===")
pass <- TRUE
for (g in BENCH) {
  r <- res[res$accession == g, , drop = FALSE]
  if (!nrow(r)) { say("  ", g, "  not in the results  [FAIL]"); pass <- FALSE; next }
  good <- r$tissue[1] == "whole blood" && nzchar(r$group_field[1])
  if (!good) pass <- FALSE
  say(sprintf("  %-10s n=%-5d tissue=%-12s group field=%-22s %s",
              g, r$n_sample[1], r$tissue[1],
              if (nzchar(r$group_field[1])) r$group_field[1] else "(not identified)",
              if (good) "pass" else "[FAIL]"))
  if (nzchar(r$levels[1])) say("             values: ", r$levels[1])
}
say(if (pass) "  Benchmark passed completely; the list below can be trusted"
    else "  [BENCHMARK NOT FULLY PASSED] the list cannot be trusted yet; the decision logic must be fixed first")
say("")

say("=== Tissue distribution ===")
tb <- sort(table(res$tissue), decreasing = TRUE)
for (i in seq_along(tb)) say(sprintf("  %-14s %d", names(tb)[i], as.integer(tb[i])))

wb <- res[res$ok & res$tissue == "whole blood" & nzchar(res$group_field) &
            nzchar(res$disease), , drop = FALSE]
say("")
say("=== Series with whole blood, a group field and a disease term: ", nrow(wb), " ===")
tab <- list()
for (d in names(DISEASE)) {
  s <- wb[grepl(d, wb$disease, fixed = TRUE), , drop = FALSE]
  if (!nrow(s)) next
  tab[[length(tab) + 1]] <- data.frame(
    disease = d, n_series = nrow(s), total_samples = sum(s$n_sample),
    max_series = max(s$n_sample),
    accessions = paste(s$accession[order(-s$n_sample)], collapse = ", "),
    stringsAsFactors = FALSE)
}
tt <- do.call(rbind, tab)
tt <- tt[order(-tt$n_series, -tt$total_samples), , drop = FALSE]
for (i in seq_len(nrow(tt))) {
  say(sprintf("  %-28s cohorts %2d, samples in total %5d, largest set %4d",
              tt$disease[i], tt$n_series[i], tt$total_samples[i], tt$max_series[i]))
}
write.csv(tt, file.path(OUT_DIR, "disease_summary_verified.csv"), row.names = FALSE)
write.csv(wb, file.path(OUT_DIR, "wholeblood_candidates.csv"), row.names = FALSE)

say("")
say("Diseases with at least two cohorts: ", sum(tt$n_series >= 2), "")
say("Reminder: the cohort numbers of this step do not yet exclude double counting of the same study population; ")
say("      true independence must be verified with array barcodes (the lesson from PEG).")

log_stop()
