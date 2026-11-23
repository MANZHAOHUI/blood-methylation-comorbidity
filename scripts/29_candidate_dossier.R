# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - expansion to five or six disorders
# Step 29: complete dossier of the candidate cohorts, for manual adjudication
#
# The expansion requires looking at every candidate closely. The lesson of the three previous rounds of automatic decisions is clear:
# 19b let a smoking field pass as the group field, and 19d did not recognise that GSE72774 and GSE111629
# are the same people - automation failed twice on cases with a known correct answer, so
# the final list must be adjudicated manually.
#
# This step downloads no matrix. It lays out the complete metadata of each candidate from the locally cached headers:
# title, summary, overall design, submitter, linked publications, platform, sample size,
# the value distribution of every characteristic field, tissue composition, and names and sizes of the supplementary files.
# Only with these can it be judged whether a series is case-control, whether it is whole blood, and whether it shares
# individuals with other data sets of the same disease.
#
# Criteria (kept and tightened):
#   tissue        proportion of whole-blood samples >= 0.6
#   design        healthy controls; not before/after treatment, not remission vs active disease, not a subtype comparison
#   size          cases >= 50 and controls >= 30
#   independence  a population different from the other included cohorts of the same disease (submitter, publication, barcodes)
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/29_candidate_dossier.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "29_candidate_dossier"

HDR_DIR <- file.path(DATA_ROOT, "geo_headers")
OUT_DIR <- file.path(PROJ_ROOT, "output", "29_dossier")
if (!dir.exists(HDR_DIR)) stop("Directory of cached headers not found: ", HDR_DIR)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
options(timeout = 600)

GET_SIZES <- TRUE      # Whether to query the sizes of supplementary files online (a few HEAD requests per candidate)
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

CAND <- list(
  "included: Parkinson's disease" = c("GSE145361", "GSE111629", "GSE72774"),
  "included: inflammatory bowel disease" = c("GSE87648", "GSE112611", "GSE87650"),
  "schizophrenia" = c("GSE152027", "GSE84727", "GSE80417", "GSE152026",
                   "GSE147221", "GSE116378", "GSE116379", "GSE116380", "GSE41169"),
  "depression" = c("GSE198904", "GSE273970", "GSE113725", "GSE200954",
               "GSE200956", "GSE201287"),
  "Alzheimer's disease" = c("GSE338167", "GSE144858", "GSE53740", "GSE156994"),
  "HIV infection" = c("GSE107082", "GSE217633", "GSE67705", "GSE67751")
)

TIS_WHOLE <- c("whole blood", "peripheral blood", "buffy coat",
               "blood leukocyte", "leukocyte", "^blood$", "\\bblood\\b")
TIS_OTHER <- c("brain", "cortex", "cerebell", "buccal", "saliva", "placenta",
               "cord blood", "adipose", "muscle", "liver", "colon", "skin",
               "biopsy", "tumor", "tumour", "kidney", "lung", "heart", "cell line")
cls1 <- function(v, p) any(vapply(p, function(x) grepl(x, v, ignore.case = TRUE), logical(1)))

rd <- function(g) {
  f <- file.path(HDR_DIR, paste0(g, "_header.txt"))
  if (!file.exists(f)) return(NULL)
  readLines(f, warn = FALSE, encoding = "UTF-8")
}
f1 <- function(h, tag) {
  x <- grep(paste0("^", tag), h, value = TRUE)
  if (!length(x)) return(character(0))
  trimws(gsub('^"|"$', "", strsplit(x[1], "\t", fixed = TRUE)[[1]])[-1])
}
fall <- function(h, tag) {
  x <- grep(paste0("^", tag), h, value = TRUE)
  if (!length(x)) return(list())
  lapply(x, function(o) trimws(gsub('^"|"$', "", strsplit(o, "\t", fixed = TRUE)[[1]])[-1]))
}
supp_sizes <- function(g) {
  if (!GET_SIZES) return(character(0))
  stub <- sub("[0-9]{3}$", "nnn", g)
  u <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", stub, g)
  html <- tryCatch(paste(readLines(url(u), warn = FALSE), collapse = " "),
                   error = function(e) NULL)
  if (is.null(html)) return("(failed to read the supplementary file directory)")
  m <- regmatches(html, gregexpr('href="[^"]+"', html))[[1]]
  fn <- unique(basename(sub('^href="', "", sub('"$', "", m))))
  fn <- fn[grepl("\\.(gz|tar|csv|txt|xlsx?)$", fn, ignore.case = TRUE)]
  if (!length(fn)) return("(no supplementary files)")
  out <- character(0)
  cu <- Sys.which("curl")
  for (x in utils::head(fn, 8)) {
    sz <- NA_real_
    if (nzchar(cu)) {
      o <- tryCatch(system2(cu, c("-sIL", "--retry", "2", shQuote(paste0(u, x))),
                            stdout = TRUE, stderr = FALSE), error = function(e) character(0))
      cl <- grep("^content-length:", tolower(o), value = TRUE)
      if (length(cl)) sz <- as.numeric(trimws(sub("^content-length:", "",
                                                  tolower(cl[length(cl)]))))
    }
    out <- c(out, sprintf("%s (%s)", x,
                          if (is.na(sz)) "size unknown" else
                            sprintf("%.2f GB", sz / 1024^3)))
    Sys.sleep(0.3)
  }
  out
}

rows <- list()
for (grp in names(CAND)) {
  say("")
  say("##############################################################")
  say("###  ", grp)
  say("##############################################################")
  for (g in CAND[[grp]]) {
    h <- rd(g)
    say("")
    say("---------- ", g, " ----------")
    if (is.null(h)) { say("  No locally cached header; it must be fetched first"); next }
    ttl <- f1(h, "!Series_title")
    smy <- f1(h, "!Series_summary")
    dsg <- f1(h, "!Series_overall_design")
    pmid <- f1(h, "!Series_pubmed_id")
    ct <- c(f1(h, "!Series_contact_name"), f1(h, "!Series_contact_institute"))
    acc <- f1(h, "!Sample_geo_accession")
    plat <- unique(f1(h, "!Sample_platform_id"))
    say("  Title: ", paste(ttl, collapse = " "))
    if (length(smy)) say("  Summary: ", substr(paste(smy, collapse = " "), 1, 420))
    if (length(dsg)) say("  Design: ", substr(paste(dsg, collapse = " "), 1, 300))
    say("  Submitter: ", paste(ct, collapse = " / "),
        if (length(pmid)) paste0("; PubMed ", paste(pmid, collapse = ",")) else "")
    say("  Platform: ", paste(plat, collapse = ";"), "; samples ", length(acc))

    src <- f1(h, "!Sample_source_name_ch1")
    chl <- fall(h, "!Sample_characteristics_ch1")
    tl <- NULL
    for (v in chl) if (any(grepl("^tissue\\s*:", v, ignore.case = TRUE))) tl <- v
    per <- if (length(src)) src else if (!is.null(tl)) tl else character(0)
    if (!is.null(tl) && length(src) == length(tl)) per <- paste(src, tl)
    nwb <- if (length(per)) sum(vapply(per, function(x)
      cls1(x, TIS_WHOLE) && !cls1(x, TIS_OTHER), logical(1))) else NA_integer_
    frac <- if (!length(per)) NA_real_ else nwb / length(per)
    say(sprintf("  Tissue: whole-blood samples %s / %d (proportion %s)",
                ifelse(is.na(nwb), "unknown", nwb), length(per),
                ifelse(is.na(frac), "unknown", sprintf("%.2f", frac))))
    if (length(per)) {
      tb <- sort(table(per), decreasing = TRUE)
      say("    top four source values: ", paste(sprintf("%s=%d", names(tb)[1:min(4, length(tb))],
                                              as.integer(tb)[1:min(4, length(tb))]),
                                      collapse = "; "))
    }
    say("  All characteristic fields: ")
    for (v in chl) {
      vv <- v[nzchar(v) & grepl(":", v, fixed = TRUE)]
      if (!length(vv)) next
      tag <- trimws(sub(":.*$", "", vv[1]))
      val <- trimws(sub("^[^:]*:", "", vv))
      tb <- sort(table(val), decreasing = TRUE)
      if (length(tb) > 10L) {
        say(sprintf("    [%s] %d distinct values (numeric or identifiers), e.g. %s",
                    tag, length(tb), paste(utils::head(names(tb), 3), collapse = ", ")))
      } else {
        say(sprintf("    [%s] %s", tag,
                    paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; ")))
      }
    }
    smp_ttl <- f1(h, "!Sample_title")
    nbc <- length(unique(unlist(regmatches(smp_ttl, gregexpr(BARCODE, smp_ttl)))))
    say("  Sample title examples: ", paste(utils::head(smp_ttl, 3), collapse = " | "))
    say("  Array barcodes parsed: ", nbc, "", if (nbc == 0) " (overlap cannot be checked with barcodes)" else "")
    ss <- supp_sizes(g)
    say("  Supplementary files: ")
    for (x in ss) say("    ", x)

    rows[[length(rows) + 1]] <- data.frame(
      group = grp, accession = g, n_sample = length(acc),
      platform = paste(plat, collapse = ";"), wb_frac = frac,
      contact = paste(ct, collapse = " / "),
      pubmed = paste(pmid, collapse = ","), n_barcode = nbc,
      title = paste(ttl, collapse = " "), stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "candidate_dossier.csv"), row.names = FALSE)
say("")
say("Dossier written: ", file.path(OUT_DIR, "candidate_dossier.csv"))
say("Next: review this log in full and adjudicate each candidate to obtain the final list and the download manifest.")
log_stop()

