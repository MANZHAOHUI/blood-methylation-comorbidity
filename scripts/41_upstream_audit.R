# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 41: upstream audit - do samples and phenotypes really correspond?
#
# This step looks for a class of silent errors: no error, no warning, only a series of plausible-looking
# numbers. This has happened four times in this round (an outcome among the covariates, an unverified theoretical value,
# intermediate results not saved, a new script not reproducing the model of an old one), so similar upstream risks
# are checked one by one.
#
# The most serious is the correspondence between samples and phenotypes. A mapping success rate of 1.000 only proves that identifiers match,
# not that they refer to the same person. If the matrix columns of a set were shifted as a whole relative to the GSMs, the inflation factor would be normal,
# the FDR count zero and everything would look right, while all effect sizes of the cohort would be noise - which would also
# explain a failed within-disease replication; the two cannot be told apart from the results.
#
# Inferring sex from the methylation data themselves and comparing it with the recorded sex verifies this independently.
# Criterion: the proportion of X-chromosome probes with intermediate methylation (0.3 to 0.7) is clearly higher in women than in
# men - women have one inactivated X whose CpG islands are intermediately methylated; men have a single X and most
# probes are near 0 or 1.
#
# Sex-chromosome probes were filtered out of the analysis sets by D2, so the raw matrices must be read again.
# The raw matrices of the six first-round sets are under rds/, those of the six second-round sets in the archives under raw/.
#
# Two accompanying things are done:
#   (1) Random shuffling as a control. A verification that cannot fail is no verification: if the concordance after shuffling
#       is similar to that without shuffling, the test has no discriminating power in that data set and its conclusion is void.
#   (2) Summarise the validation evidence of each cohort that is independent of this pipeline and flag cohorts without any.
#       GSE273970 has age only in ten-year groups and is the only one of the ten validated on neither side.
#
# This step changes no data; it only audits.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/41_upstream_audit.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "41_upstream_audit"

RAW_DIR   <- file.path(DATA_ROOT, "raw")
RDS_DIR   <- file.path(DATA_ROOT, "rds")
RES_DIR   <- file.path(DATA_ROOT, "ewas_final")
COV_DIR   <- file.path(DATA_ROOT, "covariates")
PROBE_DIR <- file.path(DATA_ROOT, "probe_lists")
MAP_DIR   <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR   <- file.path(PROJ_ROOT, "output", "41_audit")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(data.table))

set.seed(20260918)
CONC_MIN  <- 0.95     # Lower limit of sex concordance
SEP_MIN   <- 0.10     # Minimum distance between the two cluster centres; below it the test has no discriminating power
EPICV2    <- "GSE338167"

COH <- list(
  GSE145361 = "parkinson", GSE111629 = "parkinson",
  GSE84727  = "schizophrenia", GSE80417 = "schizophrenia",
  GSE338167 = "alzheimer", GSE144858 = "alzheimer",
  GSE198904 = "depression", GSE273970 = "depression",
  GSE87648  = "ibd", GSE112611 = "ibd"
)

anno <- base::readRDS(file.path(PROBE_DIR, "anno450k.rds"))
pid <- base::as.character(base::rownames(anno))
if (!length(pid) || base::mean(base::grepl("^cg", pid)) < 0.5)
  pid <- base::as.character(anno[["Name"]])
chrX <- pid[base::as.character(anno[["chr"]]) == "chrX"]
say("X-chromosome probes in the 450K annotation: ", length(chrX))
say("")

# Translate matrix column names to GSMs: used directly if they match, otherwise looked up in the sample mapping table of 06,
# and finally by scanning the series matrix header. The rds files of the first round hold the raw matrices before mapping,
# with array barcodes rather than GSMs as column names, so a direct intersection would give zero samples.
to_gsm <- function(g, cn) {
  if (base::mean(cn %in% base::unlist(base::lapply(list(cn), function(z) z))) &&
      base::any(base::grepl("^GSM[0-9]+$", cn))) return(cn)
  f <- file.path(MAP_DIR, base::paste0(g, "_sample_map.csv"))
  if (file.exists(f)) {
    m <- read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
    if (base::all(c("gsm", "col_name") %in% base::names(m))) {
      v <- m$gsm[base::match(cn, m$col_name)]
      # Some raw matrices store beta columns and detection P value columns alternately (GSE112611: 804 columns = 402 samples x 2);
      # detection P value columns have no sample in the mapping table. The success rate is computed over beta columns only, detection P value columns
      # are translated to NA and removed by the caller.
      is_det <- base::grepl("^Detection", cn, ignore.case = TRUE)
      if (base::any(!is_det) && base::mean(!base::is.na(v[!is_det])) > 0.95) {
        say("  Column names translated to GSMs via the sample mapping table of 06, success rate ",
            sprintf("%.3f", base::mean(!base::is.na(v[!is_det]))),
            if (base::any(is_det)) base::paste0(" (plus ", base::sum(is_det), " detection P value columns, not used)") else "")
        v[is_det] <- NA
        return(v)
      }
    }
  }
  hd <- tryCatch(mr_header(RAW_DIR, g), error = function(e) NULL)
  if (!is.null(hd)) {
    mp <- mr_map_samples(cn, hd)
    if (mp$frac > 0.95) {
      say("  Column names translated to GSMs via the header field '", mp$name, "', success rate ",
          sprintf("%.3f", mp$frac))
      return(mp$gsm[base::match(cn, mp$key)])
    }
  }
  NULL
}

# Fetch raw beta: first round in rds/, second round in raw/
get_raw_beta <- function(g, want_gsm) {
  f1 <- file.path(RDS_DIR, base::paste0(g, "_beta.rds"))
  if (file.exists(f1)) {
    b <- base::readRDS(f1)
    say("  Source: rds/", base::basename(f1), " (", base::nrow(b), " x ", base::ncol(b), ")")
    gs <- to_gsm(g, base::colnames(b))
    if (is.null(gs)) stop("column names cannot be translated to GSMs")
    base::colnames(b) <- gs
    b <- b[, !base::is.na(gs), drop = FALSE]
    keep <- base::intersect(want_gsm, base::colnames(b))
    if (!length(keep)) stop("no overlap with the GSMs of the design after translation")
    b <- b[, keep, drop = FALSE]
    # What is fetched must be beta values: range between 0 and 1, and the overall mean cannot be close to 0 (that would be detection P values)
    rg <- base::range(b, na.rm = TRUE); mu <- base::mean(b, na.rm = TRUE)
    if (rg[1] < -0.05 || rg[2] > 1.05 || mu < 0.2 || mu > 0.8)
      stop("The columns fetched do not look like beta values: range ", sprintf("%.3f to %.3f", rg[1], rg[2]), ", mean ", sprintf("%.3f", mu))
    return(b)
  }
  f2 <- mr_data_file(RAW_DIR, g)
  say("  Source: raw/", base::basename(f2))
  dt <- mr_read(f2, verbose = FALSE)
  probes <- base::as.character(dt[[1]])
  cn <- base::names(dt)[-1]
  M <- base::as.matrix(dt[, -1, with = FALSE])
  rm(dt); invisible(gc())
  pk <- mr_pick_beta(M, cn, verbose = FALSE)
  M <- M[, pk$idx, drop = FALSE]; cn <- pk$names
  hd <- mr_header(RAW_DIR, g)
  mp <- mr_map_samples(cn, hd)
  if (mp$frac < 0.95) stop(g, ": mapping success rate only ", sprintf("%.3f", mp$frac))
  base::colnames(M) <- mp$gsm[base::match(cn, mp$key)]
  if (g %in% EPICV2) {
    bs <- base::sub("_.*$", "", probes)
    if (base::any(base::duplicated(bs))) {
      s <- base::rowsum(M, group = bs, na.rm = TRUE)
      n <- base::rowsum((!base::is.na(M)) * 1, group = bs, na.rm = TRUE)
      M <- s / n; M[!base::is.finite(M)] <- NA
    } else base::rownames(M) <- bs
  } else base::rownames(M) <- probes
  keep <- base::intersect(want_gsm, base::colnames(M))
  M[, keep, drop = FALSE]
}

rows <- list()
for (g in names(COH)) {
  say("==== ", g, " (", COH[[g]], ") ====")
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  want <- D$gsm
  # Recorded sex: recovered from the design matrix (dummy variables such as sexM)
  sc <- base::grep("^sex", base::colnames(D$design), value = TRUE)
  if (!length(sc)) { say("  No sex term in the design, skipped"); next }
  sex_rec <- base::ifelse(D$design[, sc[1]] == 1,
                          base::sub("^sex", "", sc[1]), "REF")
  tb <- base::table(sex_rec)
  say("  Recorded sex: ", base::paste(sprintf("%s=%d", base::names(tb),
                                            base::as.integer(tb)), collapse = "; "))
  if (length(tb) < 2L) { say("  Only one sex, verification not possible"); next }

  B <- tryCatch(get_raw_beta(g, want), error = function(e) {
    say("  [FAIL] could not fetch the raw matrix: ", conditionMessage(e)); NULL })
  if (is.null(B)) next
  ix <- base::which(base::rownames(B) %in% chrX)
  say("  X probes in the matrix: ", length(ix), "; usable samples ", base::ncol(B), " / ", length(want))
  if (length(ix) < 1000L) { say("  [FAIL] too few X probes, verification not possible"); next }

  MX <- B[ix, , drop = FALSE]
  frac_mid <- base::colMeans(MX > 0.3 & MX < 0.7, na.rm = TRUE)
  rm(B, MX); invisible(gc())
  nd <- length(base::unique(base::round(frac_mid, 6)))
  if (length(frac_mid) < 20L || nd < 3L) {
    say("  [FAIL] usable samples ", length(frac_mid), ", distinct values ", nd,
        "; two clusters cannot be formed")
    rows[[length(rows) + 1]] <- data.frame(gse = g, disease = COH[[g]],
      n = length(frac_mid), separation = NA_real_, concordance = NA_real_,
      concordance_shuffled = NA_real_, verdict = "too few samples", stringsAsFactors = FALSE)
    next
  }
  cl <- stats::kmeans(frac_mid, centers = 2, nstart = 10)
  ctr <- base::sort(cl$centers[, 1])
  sep <- base::diff(ctr)
  say("  Proportion of X probes with intermediate methylation: cluster centres ", sprintf("%.3f and %.3f", ctr[1], ctr[2]),
      ", distance ", sprintf("%.3f", sep))
  if (sep < SEP_MIN) {
    say("  [FAIL] the two clusters are too close; the test has no discriminating power in this cohort and the result is unusable")
    rows[[length(rows) + 1]] <- data.frame(gse = g, disease = COH[[g]],
      n = length(frac_mid), separation = sep, concordance = NA_real_,
      concordance_shuffled = NA_real_, verdict = "no discriminating power", stringsAsFactors = FALSE)
    next
  }
  hi <- base::which.max(cl$centers)
  pred_hi <- cl$cluster == hi          # High intermediate methylation = female
  rec <- sex_rec[base::match(base::names(frac_mid), want)]
  # Which recorded label corresponds to female is decided by the majority match; the concordance is then reported
  lv <- base::names(base::table(rec))
  conc <- base::max(base::sapply(lv, function(L)
    base::mean(pred_hi == (rec == L), na.rm = TRUE)))
  fem <- lv[base::which.max(base::sapply(lv, function(L)
    base::mean(pred_hi == (rec == L), na.rm = TRUE)))]
  shuf <- base::mean(base::sample(pred_hi) == (rec == fem), na.rm = TRUE)
  say("  Inferred female: ", base::sum(pred_hi), "; compared with the recorded '", fem,
      "', concordance ", sprintf("%.4f", conc))
  say("  Concordance after random shuffling: ", sprintf("%.4f", shuf), " (control: if close to the line above, the test is invalid)")
  vd <- if (conc >= CONC_MIN) "pass" else "[FAIL]"
  say("  Decision: ", vd, if (conc < CONC_MIN)
    "  The correspondence between samples and phenotypes is in doubt; all results of this cohort need to be re-checked" else "")
  rows[[length(rows) + 1]] <- data.frame(gse = g, disease = COH[[g]],
    n = length(frac_mid), separation = sep, concordance = conc,
    concordance_shuffled = shuf, verdict = vd, stringsAsFactors = FALSE)
  say("")
}

tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "sex_check.csv"), row.names = FALSE)
say("=== Summary of the sex check ===")
print(tab)
say("")
say("  Passed: ", base::sum(tab$verdict == "pass", na.rm = TRUE), " / ", base::nrow(tab))
if (base::any(tab$verdict != "pass", na.rm = TRUE))
  say("  [ACTION NEEDED] ", base::paste(tab$gse[tab$verdict != "pass"], collapse = ", "))

# --------------------------------------------------------------
say("")
say("=== Validation evidence of each cohort that is independent of this pipeline ===")
say("  'Independent' means that the evidence does not rest on the internal consistency of this pipeline but agrees with external facts.")
EV <- data.frame(rbind(
  c("GSE145361", "agrees with the published summary statistics of Vallerga et al. (rho=0.581)"),
  c("GSE111629", "DNAm age agrees with recorded age"),
  c("GSE87648",  "DNAm age agrees with recorded age; higher neutrophils in cases fit inflammatory bowel disease"),
  c("GSE112611", "DNAm age agrees with recorded age (narrow age range, weaker evidence); higher neutrophils in cases"),
  c("GSE84727",  "DNAm age agrees with recorded age; higher smoking score in cases, in line with the literature"),
  c("GSE80417",  "DNAm age agrees with recorded age; higher smoking score in cases, in line with the literature"),
  c("GSE198904", "DNAm age agrees with recorded age"),
  c("GSE273970", "none: age only in ten-year groups, no external comparison so far"),
  c("GSE338167", "DNAm age agrees with recorded age"),
  c("GSE144858", "DNAm age agrees with recorded age")
), stringsAsFactors = FALSE)
names(EV) <- c("gse", "evidence")
EV <- merge(EV, tab[, c("gse", "concordance", "verdict")], by = "gse", all.x = TRUE)
for (i in base::seq_len(base::nrow(EV)))
  say(sprintf("  %-10s sex check %-8s  %s", EV$gse[i],
              if (base::is.na(EV$concordance[i])) "not finished"
              else sprintf("%.4f", EV$concordance[i]), EV$evidence[i]))
write.csv(EV, file.path(OUT_DIR, "independent_evidence.csv"), row.names = FALSE)

say("")
say("  Note: agreement of DNAm age with recorded age verifies two things at once - the correspondence between probe names and values")
say("  is correct (the clock takes coefficients by probe; a shift would not give ages), and the correspondence between samples and phenotypes is correct")
say("  (the ages compared come from the records). GSE273970 has no external comparison on either side; the sex check of this step")
say("  is its only independent evidence.")

# --------------------------------------------------------------
say("")
say("=== GSE273970: does DNAm age agree with the ten-year age groups? ===")
say("  This cohort has no numeric age but ten-year groups, which allow a test of the rank order.")
g <- "GSE273970"
hd <- mr_header(RAW_DIR, g)
gsm <- mr_split(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
band <- NULL
for (one in base::grep("^!Sample_characteristics_ch1", hd, value = TRUE)) {
  v <- mr_split(one)[-1]
  vv <- v[nzchar(v) & base::grepl(":", v, fixed = TRUE)]
  if (!length(vv)) next
  if (base::grepl("age", base::trimws(base::sub(":.*$", "", vv[1])), ignore.case = TRUE)) {
    band <- base::trimws(base::sub("^[^:]*:", "", v)); break
  }
}
if (is.null(band)) { say("  Age group field not found") } else {
  cov <- read.csv(file.path(COV_DIR, base::paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  bb <- band[base::match(cov$gsm, gsm)]
  ok <- !base::is.na(bb) & base::is.finite(cov$dnam_age_horvath)
  tb <- base::table(bb[ok])
  say("  Distribution of groups: ", base::paste(sprintf("%s=%d", base::names(tb),
                                          base::as.integer(tb)), collapse = "; "))
  mid <- base::sapply(base::strsplit(base::gsub("[^0-9-]", "", bb[ok]), "-"),
                      function(z) base::mean(base::as.numeric(z)))
  r <- stats::cor(mid, cov$dnam_age_horvath[ok], method = "spearman")
  say("  Spearman correlation between group midpoints and DNAm age = ", sprintf("%.3f", r))
  for (L in base::names(tb))
    say(sprintf("    %-10s n=%3d  median DNAm age %.1f", L, tb[[L]],
                stats::median(cov$dnam_age_horvath[ok][bb[ok] == L], na.rm = TRUE)))
  say("  Decision: ", if (r >= 0.6) "pass, rank order agrees"
      else "[FAIL] rank order does not agree; the sample correspondence of this cohort must be re-checked")
}
log_stop()
