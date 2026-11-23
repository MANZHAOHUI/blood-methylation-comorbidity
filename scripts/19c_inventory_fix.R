# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - cross-disease atlas
# Step 19c: close the three gaps exposed by 19b
#
# The tissue decision of 19b passed the benchmark (all six recovered and classified as whole blood), but two other
# points failed, the second of which the benchmark itself did not catch:
#
#   Gap 1   for 300 of 947 series no header was obtained (32%). The reason: file names were guessed
#           with only four suffixes tried, whereas file names in the GEO matrix directories are not uniform.
#           The missed ones include real candidates such as GSE87650 (the IBD
#           super-series of Ventham et al.) and GSE56047 (MESA monocytes). The directory contents are now listed
#           directly instead of guessing.
#
#   Gap 2   the wrong group field was taken. GSE87648 was assigned "smoking status" -
#           because the rule was "take the first field whose tag contains status", and
#           smoking status happens to contain status. Its true group field is
#           simplified_diagnosis. The benchmark let it pass because it only required
#           "a group field exists", not "the right one". A criterion that loose decides nothing.
#           Each candidate field is now scored: values containing terms such as control/healthy
#           add points, tags such as smoking/sex/age/batch subtract points,
#           the highest score wins, and all candidate fields are written out for manual review.
#
#   Gap 3   150 series were classified as unclear. The bare word blood is added and the tissue: tag is read.
#
# The benchmark is tightened at the same time: the six sets must not only be recovered and classified as whole blood, the group field
# must also be the correct one. GSE87648 must return the diagnosis field, not the smoking field, or it does not pass.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/19c_inventory_fix.R
# Note: cached headers are not fetched again; only the 300 that 19b did not obtain are added
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "19c_inventory_fix"

IN_DIR  <- file.path(PROJ_ROOT, "output", "19_inventory")
OUT_DIR <- IN_DIR
HDR_DIR <- file.path(DATA_ROOT, "geo_headers")
if (!dir.exists(HDR_DIR)) dir.create(HDR_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
options(timeout = 600)

MIN_N <- 80L

# The benchmark now also requires the correct group field: expected tag keywords are given
BENCH <- list(
  GSE145361 = "disease",  GSE111629 = "disease", GSE72774  = "disease",
  GSE87648  = "diagnos",  GSE112611 = "diagnos", GSE152027 = "status")

TIS_WHOLE  <- c("whole blood", "peripheral blood", "buffy coat",
                "blood leukocyte", "leukocyte", "^blood$", "\\bblood\\b")
TIS_PBMC   <- c("pbmc", "mononuclear")
TIS_SORTED <- c("cd4", "cd8", "cd14", "cd19", "monocyte", "granulocyte",
                "neutrophil", "t cell", "b cell", "nk cell", "sorted")
TIS_OTHER  <- c("brain", "cortex", "buccal", "saliva", "placenta", "cord blood",
                "adipose", "muscle", "liver", "colon", "skin", "biopsy",
                "tumor", "tumour", "kidney", "lung", "heart", "cell line")

CTRL_WORDS <- "control|healthy|normal|\\bHC\\b|\\bNC\\b|unaffected|non-|negative|baseline"
GOOD_TAG   <- "diagnos|disease|status|group|condition|phenotype|affection|subject|case"
BAD_TAG    <- "smok|sex|gender|age|batch|slide|sentrix|array|position|ethnic|race|bmi|tissue|cell.?type|plate|well|chip|scan|storage|time|visit|treatment|medication|drug"

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
  "COVID"                      = c("COVID", "SARS-CoV-2")
)
hit_any <- function(t, p) any(vapply(p, function(x)
  grepl(x, t, ignore.case = TRUE), logical(1)))
cls <- function(v, pats) any(vapply(pats, function(p)
  any(grepl(p, v, ignore.case = TRUE)), logical(1)))

# --------------------------------------------------------------
# List the matrix directory instead of guessing file names
# --------------------------------------------------------------
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

fetch_header <- function(gse) {
  cache <- file.path(HDR_DIR, paste0(gse, "_header.txt"))
  if (file.exists(cache)) return(readLines(cache, warn = FALSE, encoding = "UTF-8"))
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
  NULL
}

fld <- function(lines, tag) {
  x <- grep(paste0("^", tag), lines, value = TRUE)
  if (!length(x)) return(character(0))
  unlist(lapply(x, function(one)
    trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1])))
}

# --------------------------------------------------------------
# Group field: scored rather than taking the first
# --------------------------------------------------------------
pick_group <- function(lines) {
  ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  best <- list(score = -99, tag = "", levels = "")
  allc <- character(0)
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
    # Two or three levels with a smallest group that is not too small look more like case-control
    if (length(tb) <= 3L) sc <- sc + 1
    if (min(tb) >= 20L) sc <- sc + 1
    lv <- paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")
    allc <- c(allc, sprintf("[%s|score %d] %s", tag, sc, lv))
    if (sc > best$score) best <- list(score = sc, tag = tag, levels = lv)
  }
  list(tag = best$tag, levels = best$levels, score = best$score,
       all = paste(allc, collapse = " || "))
}

# --------------------------------------------------------------
all <- read.csv(file.path(IN_DIR, "geo_all_series.csv"), stringsAsFactors = FALSE)
pool <- all[!is.na(all$n_sample) & all$n_sample >= MIN_N, , drop = FALSE]
pool <- pool[order(-pool$n_sample), , drop = FALSE]
pool <- rbind(pool[pool$accession %in% names(BENCH), , drop = FALSE],
              pool[!pool$accession %in% names(BENCH), , drop = FALSE])
cached <- length(list.files(HDR_DIR, pattern = "_header\\.txt$"))
say("To process: ", nrow(pool), ", cached headers ", cached, "")
say("")

rows <- list(); t0 <- Sys.time()
for (i in seq_len(nrow(pool))) {
  g <- pool$accession[i]
  h <- fetch_header(g)
  if (is.null(h)) {
    rows[[length(rows) + 1]] <- data.frame(
      accession = g, n_sample = pool$n_sample[i], tissue = "header not available",
      group_field = "", levels = "", score = NA_integer_, all_fields = "",
      disease = "", title = pool$title[i], stringsAsFactors = FALSE)
    next
  }
  src <- c(fld(h, "!Sample_source_name_ch1"), fld(h, "!Sample_characteristics_ch1"))
  src <- src[nzchar(src)]
  tis <- if (cls(src, TIS_OTHER) && !cls(src, TIS_WHOLE)) "non-blood"
  else if (cls(src, TIS_WHOLE)) "whole blood"
  else if (cls(src, TIS_PBMC)) "PBMC"
  else if (cls(src, TIS_SORTED)) "sorted cells"
  else "unclear"

  gp <- pick_group(h)
  txt <- paste(pool$title[i], pool$summary[i], paste(utils::head(src, 40), collapse = " "))
  dz <- names(DISEASE)[vapply(DISEASE, function(p) hit_any(txt, p), logical(1))]

  rows[[length(rows) + 1]] <- data.frame(
    accession = g, n_sample = pool$n_sample[i], tissue = tis,
    group_field = gp$tag, levels = gp$levels, score = gp$score,
    all_fields = gp$all, disease = paste(dz, collapse = "|"),
    title = pool$title[i], stringsAsFactors = FALSE)

  if (i %% 100 == 0) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("  ", i, " / ", nrow(pool), "  Time: ", sprintf("%.1f", el),
        " min, estimated remaining ", sprintf("%.0f", el / i * (nrow(pool) - i)), " min")
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(OUT_DIR, "series_verified_v2.csv"), row.names = FALSE)

# --------------------------------------------------------------
say("")
say("=== Benchmark test (now also requiring the correct group field) ===")
pass <- TRUE
for (g in names(BENCH)) {
  r <- res[res$accession == g, , drop = FALSE]
  if (!nrow(r)) { say("  ", g, "  not present  [FAIL]"); pass <- FALSE; next }
  ok_t <- r$tissue[1] == "whole blood"
  ok_g <- grepl(BENCH[[g]], r$group_field[1], ignore.case = TRUE)
  if (!(ok_t && ok_g)) pass <- FALSE
  say(sprintf("  %-10s tissue=%-12s group field=%-24s %s", g, r$tissue[1],
              if (nzchar(r$group_field[1])) r$group_field[1] else "(not identified)",
              if (ok_t && ok_g) "pass" else
                paste0("[FAIL, expected tag to contain ", BENCH[[g]], "]")))
  say("             values: ", r$levels[1])
}
say(if (pass) "  Benchmark passed completely" else "  [BENCHMARK FAILED] the list below cannot be trusted")
say("")

say("=== Tissue distribution ===")
tb <- sort(table(res$tissue), decreasing = TRUE)
for (i in seq_along(tb)) say(sprintf("  %-14s %d", names(tb)[i], as.integer(tb[i])))

wb <- res[res$tissue == "whole blood" & nzchar(res$group_field) &
            nzchar(res$disease) & !is.na(res$score) & res$score >= 2,
          , drop = FALSE]
say("")
say("=== Whole blood, credible group field (score >= 2), disease term: ", nrow(wb), " ===")
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
  say("        ", tt$accessions[i])
}
write.csv(tt, file.path(OUT_DIR, "disease_summary_v2.csv"), row.names = FALSE)
write.csv(wb, file.path(OUT_DIR, "wholeblood_candidates_v2.csv"), row.names = FALSE)

say("")
say("Diseases with at least two cohorts: ", sum(tt$n_series >= 2), "")
say("This is still an upper bound: double counting of the same study population is not yet excluded and must be verified pair by pair with barcodes.")
log_stop()
