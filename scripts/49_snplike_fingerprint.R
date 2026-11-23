# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 49: sample fingerprints from SNP-like probes, to cover cohorts without rs probes
#
# Only two of the ten sets kept the rs probes of the array, and the two pairs most in need of checking are not among them:
# GSE84727 vs GSE80417 for schizophrenia and GSE145361 vs GSE111629 for Parkinson's disease.
# If these pairs shared individuals, the within-disease positive control would not measure "reproducible" but "the same people
# counted twice".
#
# This can be checked without downloading IDATs. A considerable number of probes on the array are affected by a single-nucleotide polymorphism
# underneath them, and their beta values form three clusters (near 0, 0.5 and 1) corresponding to the three genotypes. Such
# probes can be recognised directly in the beta matrix and used as a sample fingerprint.
#
# The rule is fixed beforehand: a probe is SNP-like if, in the cohort, the proportions of samples in the three intervals (below 0.25, 0.25 to
# 0.75, above 0.75) are each at least 5%.
# For a comparison of two sets the probes common to both are taken, and a decision requires at least 200.
#
# The key is calibration rather than direct use: GSE87648 and GSE80417 have rs probes and
# known answers (2 duplicate pairs within the former, none in the latter). This step first verifies in these two sets
# whether SNP-like probes reproduce the same answers - if calibration fails the method is not used further and IDATs are downloaded instead.
#
# Calibration is by rank, not by absolute threshold. rs probes and SNP-like probes are two different sets of sites
# with different allele frequency spectra, so the absolute level of concordance differs, and transferring a threshold calibrated on one set
# to the other would be a category error. Two criteria are used instead, both fixed beforehand:
#   (1) the concordances given by the two kinds of probes correlate at a Spearman coefficient of at least 0.80 over all sample pairs;
#   (2) the duplicate pairs identified by rs probes must occupy the top ranks under the SNP-like probes.
# The number of probes required is not preset either; it follows from the calibration: separating the known answers means sufficient.
#
# These probes are removed from the formal analysis by rule D2 (minor allele frequency at CpG or SBE
# above 0.01), so they are taken from the raw matrices, not from the analysis sets.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/49_snplike_fingerprint.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "49_snplike_fingerprint"

RAW_DIR <- file.path(DATA_ROOT, "raw")
RDS_DIR <- file.path(DATA_ROOT, "rds")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
ID_DIR  <- file.path(PROJ_ROOT, "output", "48_identity")
OUT_DIR <- file.path(PROJ_ROOT, "output", "49_fingerprint")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(data.table))

MIN_FRAC  <- 0.05   # Minimum proportion of samples in each of the three intervals
END_LO    <- 0.15   # Upper limit of the mean of the low cluster (homozygotes should be near 0)
END_HI    <- 0.85   # Lower limit of the mean of the high cluster (homozygotes should be near 1)
MID_LO    <- 0.40   # Lower limit of the mean of the middle cluster (heterozygotes should be near 0.5)
MID_HI    <- 0.60   # Upper limit of the mean of the middle cluster
CL_SD     <- 0.06   # Upper limit of the within-cluster SD of each of the three clusters (genotype-driven clusters should be tight)
MIN_SHARE <- 50L    # Lower limit of SNP-like probes shared by two sets (usability is decided by calibration; this only guards against extreme scarcity)
CAL_RHO   <- 0.80   # Calibration requirement 1: correlation of the concordances of the two kinds of probes over all sample pairs
TOP_MULT  <- 1.15   # Calibration requirement 2: the concordance of known duplicate pairs must reach this multiple of the highest remaining value

# Combinations to check: the two most critical within-disease pairs, plus two cohorts with known answers for calibration
CAL  <- c("GSE87648", "GSE80417")
WANT <- list(c("GSE84727", "GSE80417"), c("GSE145361", "GSE111629"))
NEED <- base::unique(c(CAL, base::unlist(WANT)))

to_gsm <- function(g, cn) {
  if (base::mean(base::grepl("^GSM[0-9]+$", cn)) > 0.95) return(cn)
  f <- file.path(MAP_DIR, base::paste0(g, "_sample_map.csv"))
  if (file.exists(f)) {
    m <- read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
    if (base::all(c("gsm", "col_name") %in% base::names(m))) {
      v <- m$gsm[base::match(cn, m$col_name)]
      if (base::mean(!base::is.na(v)) > 0.95) return(v)
    }
  }
  hd <- tryCatch(mr_header(RAW_DIR, g), error = function(e) NULL)
  if (!is.null(hd)) {
    mp <- mr_map_samples(cn, hd)
    if (mp$frac > 0.95) return(mp$gsm[base::match(cn, mp$key)])
  }
  NULL
}
get_beta <- function(g) {
  f1 <- file.path(RDS_DIR, base::paste0(g, "_beta.rds"))
  if (file.exists(f1)) {
    b <- base::readRDS(f1)
    gs <- to_gsm(g, base::colnames(b))
    if (is.null(gs)) stop("column names cannot be translated to sample IDs")
    base::colnames(b) <- gs
    return(b[, !base::is.na(gs), drop = FALSE])
  }
  f2 <- mr_data_file(RAW_DIR, g)
  dt <- mr_read(f2, verbose = FALSE)
  probes <- base::as.character(dt[[1]])
  cn <- base::names(dt)[-1]
  M <- base::as.matrix(dt[, -1, with = FALSE]); rm(dt); invisible(gc())
  pk <- mr_pick_beta(M, cn, verbose = FALSE)
  M <- M[, pk$idx, drop = FALSE]; cn <- pk$names
  base::rownames(M) <- base::sub("_.*$", "", probes)
  hd <- mr_header(RAW_DIR, g)
  mp <- mr_map_samples(cn, hd)
  if (mp$frac < 0.95) stop("mapping success rate only ", sprintf("%.3f", mp$frac))
  base::colnames(M) <- mp$gsm[base::match(cn, mp$key)]
  M[, !base::is.na(base::colnames(M)), drop = FALSE]
}
# Requiring "three modes" alone is not enough: technical batch, cell composition and X inactivation can all create three modes,
# and these are shared between samples processed together, which would raise the concordance of samples from the same batch.
# Genotype-driven trimodality has an additional shape signature - each of the three clusters must be tight, the outer ones near 0 and 1,
# and the middle one near 0.5. The four conditions below are required together.
snp_like <- function(B) {
  lo <- B < 0.25; mi <- B >= 0.25 & B <= 0.75; hi <- B > 0.75
  f_lo <- base::rowMeans(lo, na.rm = TRUE)
  f_mi <- base::rowMeans(mi, na.rm = TRUE)
  f_hi <- base::rowMeans(hi, na.rm = TRUE)
  cluster_stat <- function(mask) {
    X <- B; X[!mask] <- NA_real_
    m <- base::rowMeans(X, na.rm = TRUE)
    v <- base::rowMeans((X - m)^2, na.rm = TRUE)
    list(m = m, sd = base::sqrt(v))
  }
  L <- cluster_stat(lo); M <- cluster_stat(mi); H <- cluster_stat(hi)
  ok <- f_lo >= MIN_FRAC & f_mi >= MIN_FRAC & f_hi >= MIN_FRAC &
        base::is.finite(L$m) & base::is.finite(M$m) & base::is.finite(H$m) &
        L$m < END_LO & H$m > END_HI &
        M$m > MID_LO & M$m < MID_HI &
        L$sd < CL_SD & M$sd < CL_SD & H$sd < CL_SD
  ok[base::is.na(ok)] <- FALSE
  ok
}
call_gt <- function(B) {
  G <- base::matrix(NA_integer_, base::nrow(B), base::ncol(B), dimnames = base::dimnames(B))
  G[B < 0.25] <- 0L; G[B >= 0.25 & B <= 0.75] <- 1L; G[B > 0.75] <- 2L
  G
}
concord <- function(A, B, selfpair) {
  cp <- base::intersect(base::rownames(A), base::rownames(B))
  if (length(cp) < MIN_SHARE) return(NULL)
  A <- A[cp, , drop = FALSE]; B <- B[cp, , drop = FALSE]
  same <- base::matrix(0, base::ncol(A), base::ncol(B))
  for (k in 0:2) {
    Ak <- (A == k) * 1; Ak[base::is.na(Ak)] <- 0
    Bk <- (B == k) * 1; Bk[base::is.na(Bk)] <- 0
    same <- same + base::crossprod(Ak, Bk)
  }
  cnt <- base::crossprod((!base::is.na(A)) * 1, (!base::is.na(B)) * 1)
  C <- same / cnt
  base::dimnames(C) <- list(base::colnames(A), base::colnames(B))
  if (selfpair) C[base::upper.tri(C, diag = TRUE)] <- NA_real_
  list(C = C, n_probe = length(cp))
}

GT <- list(); NP <- integer(0)
say("=== 1. Identifying SNP-like probes set by set ===")
for (g in NEED) {
  say("---- ", g, " ----")
  B <- tryCatch(get_beta(g), error = function(e) { say("  Could not be fetched: ", conditionMessage(e)); NULL })
  if (is.null(B)) next
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  B <- B[, base::intersect(D$gsm, base::colnames(B)), drop = FALSE]
  sel <- snp_like(B)
  say("  Full matrix ", base::nrow(B), " rows x ", base::ncol(B), " samples; probes with genotype shape: ",
      base::sum(sel, na.rm = TRUE), " ( ",
      sprintf("%.3f%%", 100 * base::mean(sel, na.rm = TRUE)), ")")
  GT[[g]] <- call_gt(B[base::which(sel), , drop = FALSE]); NP[g] <- base::sum(sel, na.rm = TRUE)
  rm(B); invisible(gc())
}

say("")
say("=== 2. Calibration on the two sets with known answers ===")
say("  Calibration is by rank: do the two kinds of probes give concordances in the same order over all sample pairs, ")
say("  and do the duplicate pairs identified by rs probes rank at the top under the SNP-like probes?")
known <- tryCatch(read.csv(file.path(ID_DIR, "identity_matches.csv"),
                           stringsAsFactors = FALSE), error = function(e) NULL)
rs_gt <- function(g) {
  B <- tryCatch(get_beta(g), error = function(e) NULL)
  if (is.null(B)) return(NULL)
  ix <- base::grep("^rs", base::rownames(B))
  if (!length(ix)) return(NULL)
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  B <- B[ix, base::intersect(D$gsm, base::colnames(B)), drop = FALSE]
  call_gt(B)
}
cal_ok <- TRUE
for (g in CAL) {
  if (is.null(GT[[g]])) { say("  ", g, ": no SNP-like probes obtained, calibration not possible"); cal_ok <- FALSE; next }
  R <- rs_gt(g)
  if (is.null(R)) { say("  ", g, ": no rs probes obtained, calibration not possible"); cal_ok <- FALSE; next }
  a1 <- concord(GT[[g]], GT[[g]], TRUE)
  a2 <- concord(R, R, TRUE)
  if (is.null(a1) || is.null(a2)) {
    say("  ", g, ": number of probes below ", MIN_SHARE, ", calibration not possible"); cal_ok <- FALSE; next
  }
  cm <- base::intersect(base::rownames(a1$C), base::rownames(a2$C))
  X <- a1$C[cm, cm]; Y <- a2$C[cm, cm]
  ok <- base::is.finite(X) & base::is.finite(Y)
  rr <- stats::cor(X[ok], Y[ok], method = "spearman")
  say(sprintf("  %-10s SNP-like probes %5d, rs probes %3d, comparable sample pairs %7d",
              g, a1$n_probe, a2$n_probe, base::sum(ok)))
  say(sprintf("    Spearman correlation of the concordances given by the two kinds of probes = %.3f %s",
              rr, if (rr >= CAL_RHO) "" else "  [FAIL]"))
  if (rr < CAL_RHO) cal_ok <- FALSE
  kn <- if (is.null(known)) NULL else known[known$a == g & known$b == g, , drop = FALSE]
  if (!is.null(kn) && base::nrow(kn)) {
    vals <- base::vapply(base::seq_len(base::nrow(kn)), function(i) {
      x <- kn$sample_a[i]; y <- kn$sample_b[i]
      v <- c(X[x, y], X[y, x]); v <- v[base::is.finite(v)]
      if (length(v)) base::max(v) else NA_real_
    }, numeric(1))
    others <- X[ok]
    for (i in base::seq_len(base::nrow(kn)))
      others <- others[others != vals[i]]
    top_other <- base::max(others, na.rm = TRUE)
    say(sprintf("    the %d duplicate pairs identified by rs probes have concordances of %s under the SNP-like probes; highest of the rest %.3f",
                base::nrow(kn), base::paste(sprintf("%.3f", vals), collapse = ", "), top_other))
    pass <- base::all(base::is.finite(vals)) && base::all(vals >= TOP_MULT * top_other)
    say("    clearly above the rest (reaching ", TOP_MULT, " times): ", if (pass) "yes" else "no  [FAIL]")
    if (!pass) cal_ok <- FALSE
  } else say("    no known duplicates in this cohort; calibrated by correlation only")
}

say("")
if (!cal_ok) {
  say("  [CALIBRATION FAILED] the SNP-like probes did not reproduce the known answers; the method is not used to judge the other cohorts.")
  say("  To continue, the IDAT raw files of GSE84727, GSE145361 and GSE111629 would have to be downloaded.")
  log_stop()
} else {
  say("  Calibration passed: the method gives consistent results on the known answers and can be used for the other cohorts.")

say("")
say("=== 3. Checking the two critical within-disease pairs ===")
rows <- list()
for (pr in WANT) {
  a <- pr[1]; b <- pr[2]
  if (is.null(GT[[a]]) || is.null(GT[[b]])) {
    say(sprintf("  %-10s vs %-10s one side not obtained, skipped", a, b)); next
  }
  r <- concord(GT[[a]], GT[[b]], FALSE)
  if (is.null(r)) {
    say(sprintf("  %-10s vs %-10s fewer than %d shared SNP-like probes, no decision possible", a, b, MIN_SHARE)); next
  }
  v <- r$C[base::is.finite(r$C)]
  # The criterion follows the distribution of the pair itself: only pairs beyond the 99.99th percentile and reaching TOP_MULT times that value
  # are considered the same person. Absolute thresholds calibrated on another probe set are not carried over.
  base_hi <- stats::quantile(v, 0.9999)
  cut <- TOP_MULT * base_hi
  hit <- base::which(r$C >= cut, arr.ind = TRUE)
  say(sprintf("  %-10s vs %-10s shared probes %5d  compared %8d pairs  median %.3f  99.99th percentile %.3f  maximum %.3f",
              a, b, r$n_probe, length(v), stats::median(v), base_hi, base::max(v)))
  say(sprintf("    bound for calling the same person %.3f; pairs meeting it %d", cut, base::nrow(hit)))
  if (base::nrow(hit)) {
    say("    [ACTION NEEDED] this within-disease pair shares individuals; its positive control must be recomputed after removing duplicates: ")
    for (i in base::seq_len(base::min(20L, base::nrow(hit))))
      say(sprintf("      %-12s and %-12s concordance %.3f",
                  base::rownames(r$C)[hit[i, 1]], base::colnames(r$C)[hit[i, 2]],
                  r$C[hit[i, 1], hit[i, 2]]))
    write.csv(data.frame(a = a, b = b,
                         sample_a = base::rownames(r$C)[hit[, 1]],
                         sample_b = base::colnames(r$C)[hit[, 2]],
                         concordance = r$C[hit], stringsAsFactors = FALSE),
              file.path(OUT_DIR, base::paste0("shared_", a, "_", b, ".csv")), row.names = FALSE)
  } else say("    No shared individuals found; the positive control of this pair is unaffected.")
  rows[[length(rows) + 1]] <- data.frame(
    a = a, b = b, n_probe = r$n_probe, n_compare = length(v),
    median = stats::median(v), max = base::max(v), n_same = base::nrow(hit),
    stringsAsFactors = FALSE)
}
if (length(rows)) {
  tb <- base::do.call(rbind, rows)
  write.csv(tb, file.path(OUT_DIR, "key_pairs_identity.csv"), row.names = FALSE)
  say("")
  say("=== Summary ===")
  say("  Critical within-disease pairs checked: ", base::nrow(tb), "; with shared individuals: ",
      base::sum(tb$n_same), "")
  if (base::sum(tb$n_same) == 0)
    say("  The two most critical within-disease positive controls are not affected by shared individuals; their conclusions stand.")
}
log_stop()
}
