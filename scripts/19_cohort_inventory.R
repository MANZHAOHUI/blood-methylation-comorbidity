# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - cross-disease atlas
# Step 19: inventory of candidate cohorts (stage 1: discovery)
#
# Purpose: establish for how many diseases at least two mutually independent public case-control cohorts with
# peripheral blood methylation can be assembled. This number decides whether the atlas can carry a paper - without a second
# cohort there is no positive control, and a negative cell without a positive control cannot be interpreted; this is
# the hardest lesson of the PD x IBD round.
#
# Inclusion criteria (fixed before seeing results, not adjusted afterwards):
#   platform  450K (GPL13534) or EPIC v1 (GPL21145) / v2 (GPL33022)
#   tissue    peripheral blood. Whole blood, PBMC and sorted cells are counted separately, not merged
#   design    case-control, with a group field that can be parsed from the GEO metadata
#   size      all candidates are recorded, but only those with cases >= 50 and controls >= 30 pass the threshold
#   independence  subsets of the same study population do not count as two cohorts. At this stage suspects are flagged by study name and
#          authors; stage 2 verifies with array barcodes (the lesson from PEG)
#
# This stage only retrieves series-level metadata from GEO (accession, title, summary, sample number, platform)
# and downloads no matrix. Stage 2 reads the series matrix header of each candidate,
# parses the group field and counts cases and controls.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/19_cohort_inventory.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
SCRIPT_NAME <- "19_cohort_inventory"
OUT_DIR <- file.path(PROJ_ROOT, "output", "19_inventory")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

if (!requireNamespace("rentrez", quietly = TRUE)) {
  message("Installing rentrez")
  install.packages("rentrez")
}
if (!requireNamespace("rentrez", quietly = TRUE)) stop("rentrez could not be installed")
library(rentrez)

options(timeout = 1200)
PAUSE <- 0.4          # Without an API key NCBI allows 3 requests per second; leave a margin
BATCH <- 200L

PLATFORMS <- c(GPL13534 = "450K", GPL21145 = "EPIC v1", GPL33022 = "EPIC v2")

# Tissue keywords: a hit is required to become a candidate
BLOOD <- c("whole blood", "peripheral blood", "blood leukocyte", "leukocyte",
           "buffy coat", "PBMC", "peripheral blood mononuclear",
           "CD4", "CD8", "CD14", "monocyte", "granulocyte", "lymphocyte",
           "neutrophil", "T cell", "B cell")
WHOLE_BLOOD <- c("whole blood", "peripheral blood", "buffy coat",
                 "blood leukocyte")

# Exclusions: these are almost certainly not the case-control design we want
EXCLUDE <- c("cord blood", "newborn", "neonat", "placenta", "cell line",
             "in vitro", "treatment response", "time course", "longitudinal",
             "twin", "cancer cell", "tumor tissue", "biopsy", "brain",
             "buccal", "saliva", "sperm", "muscle", "adipose tissue")

# Dictionary of disease keywords. Better too broad: manual review is preferable to missing something.
DISEASE <- list(
  "inflammatory bowel disease" = c("inflammatory bowel", "crohn", "ulcerative colitis", "\\bIBD\\b"),
  "rheumatoid arthritis"       = c("rheumatoid arthritis", "\\bRA\\b"),
  "systemic lupus"             = c("lupus", "\\bSLE\\b"),
  "psoriasis"                  = c("psoriasis", "psoriatic"),
  "multiple sclerosis"         = c("multiple sclerosis", "\\bMS\\b"),
  "type 1 diabetes"            = c("type 1 diabetes", "\\bT1D\\b"),
  "type 2 diabetes"            = c("type 2 diabetes", "\\bT2D\\b"),
  "obesity"                    = c("obesity", "body mass index", "\\bBMI\\b"),
  "asthma"                     = c("asthma"),
  "atopic dermatitis"          = c("atopic dermatitis", "eczema"),
  "celiac disease"             = c("celiac", "coeliac"),
  "Sjogren"                    = c("sjogren", "sjögren"),
  "systemic sclerosis"         = c("systemic sclerosis", "scleroderma"),
  "ankylosing spondylitis"     = c("ankylosing spondylitis"),
  "Graves disease"             = c("graves"),
  "schizophrenia"              = c("schizophrenia", "psychosis"),
  "bipolar disorder"           = c("bipolar"),
  "major depression"           = c("depress"),
  "Alzheimer"                  = c("alzheimer", "dementia"),
  "Parkinson"                  = c("parkinson"),
  "ALS"                        = c("amyotrophic lateral"),
  "coronary artery disease"    = c("coronary", "myocardial infarction", "atheroscler"),
  "heart failure"              = c("heart failure", "cardiomyopathy"),
  "chronic kidney disease"     = c("chronic kidney", "renal failure", "nephropathy"),
  "COPD"                       = c("chronic obstructive", "\\bCOPD\\b"),
  "NAFLD"                      = c("fatty liver", "NAFLD", "NASH", "MASLD"),
  "osteoporosis"               = c("osteoporosis"),
  "preeclampsia"               = c("preeclampsia", "pre-eclampsia"),
  "HIV"                        = c("\\bHIV\\b"),
  "tuberculosis"               = c("tuberculosis", "\\bTB\\b"),
  "COVID"                      = c("COVID", "SARS-CoV-2")
)

hit_any <- function(txt, pats) any(vapply(pats, function(p)
  grepl(p, txt, ignore.case = TRUE), logical(1)))

# --------------------------------------------------------------
say("=== Searching candidate series ===")
all_ids <- character(0)
for (pl in names(PLATFORMS)) {
  q <- sprintf('%s[Accession] AND "Homo sapiens"[Organism] AND gse[Entry Type]', pl)
  r <- tryCatch(rentrez::entrez_search(db = "gds", term = q, retmax = 10000),
                error = function(e) NULL)
  if (is.null(r)) { say("  ", pl, " search failed"); next }
  say("  ", pl, " (", PLATFORMS[[pl]], ") series found: ", r$count)
  all_ids <- union(all_ids, r$ids)
  Sys.sleep(PAUSE)
}
say("  Candidates after de-duplication: ", length(all_ids), "")
say("")

# --------------------------------------------------------------
say("=== Retrieving series summaries ===")
recs <- list()
grp <- split(all_ids, ceiling(seq_along(all_ids) / BATCH))
for (i in seq_along(grp)) {
  s <- tryCatch(rentrez::entrez_summary(db = "gds", id = grp[[i]]),
                error = function(e) NULL)
  if (is.null(s)) { say("  round ", i, ": batch failed, skipped"); Sys.sleep(2); next }
  if (!is.null(s$accession)) s <- list(s)     # The structure differs for a single record
  for (x in s) {
    recs[[length(recs) + 1]] <- data.frame(
      accession = if (is.null(x$accession)) NA_character_ else x$accession,
      title     = if (is.null(x$title)) "" else x$title,
      summary   = if (is.null(x$summary)) "" else x$summary,
      n_sample  = if (is.null(x$n_samples)) NA_integer_ else as.integer(x$n_samples),
      platform  = if (is.null(x$gpl)) "" else paste(x$gpl, collapse = ";"),
      pdat      = if (is.null(x$pdat)) "" else x$pdat,
      stringsAsFactors = FALSE)
  }
  if (i %% 5 == 0) say("  Retrieved ", length(recs), " / ", length(all_ids))
  Sys.sleep(PAUSE)
}
df <- do.call(rbind, recs)
df <- df[!is.na(df$accession) & grepl("^GSE", df$accession), , drop = FALSE]
df <- df[!duplicated(df$accession), , drop = FALSE]
say("  Series retrieved: ", nrow(df))
say("")

# --------------------------------------------------------------
say("=== Tagging and filtering ===")
txt <- paste(df$title, df$summary)
df$blood       <- vapply(txt, hit_any, logical(1), BLOOD)
df$whole_blood <- vapply(txt, hit_any, logical(1), WHOLE_BLOOD)
df$excluded    <- vapply(txt, hit_any, logical(1), EXCLUDE)
df$disease <- vapply(txt, function(t) {
  h <- names(DISEASE)[vapply(DISEASE, function(p) hit_any(t, p), logical(1))]
  if (!length(h)) "" else paste(h, collapse = "|")
}, character(1))

say("  Mentioning blood: ", sum(df$blood), "; of which explicitly whole blood / peripheral blood: ", sum(df$whole_blood))
say("  Matching an exclusion term: ", sum(df$excluded))
say("  Matching a disease term: ", sum(nzchar(df$disease)))

cand <- df[df$blood & !df$excluded & nzchar(df$disease) &
             !is.na(df$n_sample) & df$n_sample >= 80, , drop = FALSE]
cand <- cand[order(-cand$n_sample), , drop = FALSE]
say("  Candidates (mention blood, no exclusion term, a disease term, samples >= 80): ", nrow(cand))
say("")

say("=== Number of candidate cohorts per disease ===")
tab <- list()
for (d in names(DISEASE)) {
  sub <- cand[grepl(d, cand$disease, fixed = TRUE), , drop = FALSE]
  if (!nrow(sub)) next
  tab[[length(tab) + 1]] <- data.frame(
    disease = d, n_series = nrow(sub),
    n_series_wholeblood = sum(sub$whole_blood),
    total_samples = sum(sub$n_sample),
    max_series = max(sub$n_sample),
    accessions = paste(utils::head(sub$accession[order(-sub$n_sample)], 8),
                       collapse = ", "),
    stringsAsFactors = FALSE)
}
tb <- do.call(rbind, tab)
tb <- tb[order(-tb$n_series_wholeblood, -tb$total_samples), , drop = FALSE]
for (i in seq_len(nrow(tb))) {
  say(sprintf("  %-28s candidates %2d (whole blood %2d), samples in total %5d, largest set %4d",
              tb$disease[i], tb$n_series[i], tb$n_series_wholeblood[i],
              tb$total_samples[i], tb$max_series[i]))
}

write.csv(df,   file.path(OUT_DIR, "geo_all_series.csv"), row.names = FALSE)
write.csv(cand, file.path(OUT_DIR, "geo_candidates.csv"), row.names = FALSE)
write.csv(tb,   file.path(OUT_DIR, "disease_summary.csv"), row.names = FALSE)

say("")
say("Three tables written: all series / candidate series / summary by disease")
say("Note: this stage relies entirely on keyword matching in titles and summaries and yields candidates only.")
say("True case and control numbers require stage 2, which reads the series matrix headers.")

log_stop()
