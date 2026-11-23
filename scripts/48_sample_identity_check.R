# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 48: check of sample identity across cohorts
#
# The gap to be closed: two "mutually independent" cohorts may contain the same person. The two schizophrenia
# sets belong to the same UK consortium and the same submission, and the two Parkinson's disease sets are both disease-specific cohorts.
# If they shared individuals, the within-disease positive control would not measure "reproducible" but "the same people
# counted twice" - and the positive controls are the pivot on which everything else in the paper is read.
#
# The method uses the SNP probes carried by the array (names starting with rs). These probes are designed for sample
# fingerprinting: the genotypes of one person agree between two arrays, those of different people do not.
# ewastools uses the same probes (Heiss and Just, Clinical Epigenetics,
# 2018;10:73); they applied this check to 80 GEO data sets with 8,327 samples and
# found 940 samples flagged by at least one control metric and 133 samples in 20 data sets
# with wrongly annotated sex - such problems are common in public data.
#
# Criteria fixed beforehand:
#   genotype call     beta < 0.25 is 0; 0.25 to 0.75 is 1; > 0.75 is 2
#   comparable        at least 30 rs probes shared between two sets
#   same person       genotype concordance of at least 0.90
# The complete distribution of concordance is printed, so that the position of the threshold is visible rather than a bare yes or no.
#
# The step first counts how many rs probes remain in the ten matrices. Processed matrices do not necessarily keep them;
# cohorts with too few are flagged as such, and no conclusion is forced from a handful of probes.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/48_sample_identity_check.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "48_sample_identity_check"

RAW_DIR <- file.path(DATA_ROOT, "raw")
RDS_DIR <- file.path(DATA_ROOT, "rds")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "48_identity")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(data.table))

MIN_SNP  <- 30L    # Number of shared rs probes required for comparison
SAME_CUT <- 0.90   # Genotype concordance at which two samples are called the same person
EPICV2   <- "GSE338167"

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
ALL <- base::unlist(DIS, use.names = FALSE)

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

get_snp <- function(g) {
  f1 <- file.path(RDS_DIR, base::paste0(g, "_beta.rds"))
  if (file.exists(f1)) {
    b <- base::readRDS(f1)
    ix <- base::grep("^rs", base::rownames(b))
    say("  Source rds/: full matrix ", base::nrow(b), " rows, of which rs probes ", length(ix))
    if (!length(ix)) return(NULL)
    S <- b[ix, , drop = FALSE]; rm(b); invisible(gc())
    gs <- to_gsm(g, base::colnames(S))
    if (is.null(gs)) { say("  [FAIL] column names cannot be translated to sample IDs"); return(NULL) }
    base::colnames(S) <- gs
    return(S[, !base::is.na(gs), drop = FALSE])
  }
  f2 <- mr_data_file(RAW_DIR, g)
  dt <- mr_read(f2, verbose = FALSE)
  probes <- base::as.character(dt[[1]])
  ix <- base::grep("^rs", probes)
  say("  Source raw/: full matrix ", base::nrow(dt), " rows, of which rs probes ", length(ix))
  if (!length(ix)) { rm(dt); invisible(gc()); return(NULL) }
  cn <- base::names(dt)[-1]
  M <- base::as.matrix(dt[ix, -1, with = FALSE])
  rm(dt); invisible(gc())
  pk <- mr_pick_beta(M, cn, verbose = FALSE)
  M <- M[, pk$idx, drop = FALSE]; cn <- pk$names
  base::rownames(M) <- base::sub("_.*$", "", probes[ix])
  hd <- mr_header(RAW_DIR, g)
  mp <- mr_map_samples(cn, hd)
  if (mp$frac < 0.95) { say("  [FAIL] mapping success rate ", sprintf("%.3f", mp$frac)); return(NULL) }
  base::colnames(M) <- mp$gsm[base::match(cn, mp$key)]
  M[, !base::is.na(base::colnames(M)), drop = FALSE]
}

# --------------------------------------------------------------
say("=== 1. Scanning the ten matrices for rs probes ===")
SNP <- list(); inv <- list()
for (g in ALL) {
  say("---- ", g, " ----")
  S <- tryCatch(get_snp(g), error = function(e) { say("  Could not be fetched: ", conditionMessage(e)); NULL })
  n <- if (is.null(S)) 0L else base::nrow(S)
  if (!is.null(S)) {
    D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
    S <- S[, base::intersect(D$gsm, base::colnames(S)), drop = FALSE]
    say("  After alignment with the analysis samples: ", base::nrow(S), " probes x ", base::ncol(S), " samples")
    SNP[[g]] <- S
  }
  inv[[length(inv) + 1]] <- data.frame(gse = g, n_snp = n,
    n_sample = if (is.null(S)) 0L else base::ncol(S), stringsAsFactors = FALSE)
}
iv <- base::do.call(rbind, inv)
write.csv(iv, file.path(OUT_DIR, "snp_probe_inventory.csv"), row.names = FALSE)
say("")
say("=== Inventory of rs probes ===")
for (i in base::seq_len(base::nrow(iv)))
  say(sprintf("  %-10s rs probes %3d  samples %4d  %s", iv$gse[i], iv$n_snp[i], iv$n_sample[i],
              if (iv$n_snp[i] >= MIN_SNP) "available" else "[not usable] too few probes"))
usable <- iv$gse[iv$n_snp >= MIN_SNP]
say("  Usable cohorts: ", length(usable), " / ", base::nrow(iv))
if (length(usable) < 2L) {
  say("")
  say("  Fewer than two usable cohorts; no cross-cohort comparison possible.")
  say("  To continue, the IDAT raw files of the cohorts concerned would have to be downloaded to recover the SNP probes.")
  log_stop()
} else {

# --------------------------------------------------------------
call_gt <- function(B) {
  G <- base::matrix(NA_integer_, base::nrow(B), base::ncol(B),
                    dimnames = base::dimnames(B))
  G[B < 0.25] <- 0L; G[B >= 0.25 & B <= 0.75] <- 1L; G[B > 0.75] <- 2L
  G
}
say("")
say("=== 2. Pairwise comparison of sample genotypes ===")
say("  Criterion for the same person: at least ", MIN_SNP, " shared rs probes and a genotype concordance of at least ", SAME_CUT)
GT <- base::lapply(SNP[usable], call_gt)
pairs <- list()
for (i in 1:(length(usable) - 1)) for (j in (i + 1):length(usable))
  pairs[[length(pairs) + 1]] <- c(usable[i], usable[j])
for (g in usable) pairs[[length(pairs) + 1]] <- c(g, g)

hits <- list(); summ <- list()
for (pr in pairs) {
  a <- pr[1]; b <- pr[2]
  A <- GT[[a]]; B <- GT[[b]]
  cp <- base::intersect(base::rownames(A), base::rownames(B))
  if (length(cp) < MIN_SNP) {
    say(sprintf("  %-10s vs %-10s only %d shared rs probes, skipped", a, b, length(cp))); next
  }
  A <- A[cp, , drop = FALSE]; B <- B[cp, , drop = FALSE]
  same <- base::matrix(0, base::ncol(A), base::ncol(B))
  cnt  <- base::matrix(0, base::ncol(A), base::ncol(B))
  for (k in 0:2) {
    Ak <- (A == k) * 1; Ak[base::is.na(Ak)] <- 0
    Bk <- (B == k) * 1; Bk[base::is.na(Bk)] <- 0
    same <- same + base::crossprod(Ak, Bk)
  }
  Av <- (!base::is.na(A)) * 1; Bv <- (!base::is.na(B)) * 1
  cnt <- base::crossprod(Av, Bv)
  conc <- same / cnt
  base::dimnames(conc) <- list(base::colnames(A), base::colnames(B))
  if (a == b) conc[base::upper.tri(conc, diag = TRUE)] <- NA_real_
  v <- conc[base::is.finite(conc)]
  n_hit <- base::sum(v >= SAME_CUT)
  say(sprintf("  %-10s vs %-10s shared probes %3d  compared %8d pairs  median concordance %.3f  upper quartile %.3f  maximum %.3f  meeting criterion %d",
              a, b, length(cp), length(v), stats::median(v),
              stats::quantile(v, 0.75), base::max(v), n_hit))
  summ[[length(summ) + 1]] <- data.frame(
    a = a, b = b, kind = if (a == b) "within cohort" else "cross-cohort",
    n_probe = length(cp), n_compare = length(v), median = stats::median(v),
    q75 = stats::quantile(v, 0.75), max = base::max(v), n_same = n_hit,
    stringsAsFactors = FALSE)
  if (n_hit > 0) {
    w <- base::which(conc >= SAME_CUT, arr.ind = TRUE)
    hits[[length(hits) + 1]] <- data.frame(
      a = a, b = b, sample_a = base::rownames(conc)[w[, 1]],
      sample_b = base::colnames(conc)[w[, 2]],
      concordance = conc[w], stringsAsFactors = FALSE)
  }
}

st <- base::do.call(rbind, summ)
write.csv(st, file.path(OUT_DIR, "pairwise_identity_summary.csv"), row.names = FALSE)
say("")
say("=== Summary ===")
say("  Cohort combinations compared: ", base::nrow(st), "; range of median concordance ",
    sprintf("%.3f to %.3f", base::min(st$median), base::max(st$median)))
tot <- base::sum(st$n_same)
if (length(hits)) {
  hh <- base::do.call(rbind, hits)
  write.csv(hh, file.path(OUT_DIR, "identity_matches.csv"), row.names = FALSE)
  say("")
  say("  [ACTION NEEDED] sample pairs meeting the same-person criterion: ", base::nrow(hh), "")
  wd <- hh[hh$a != hh$b, , drop = FALSE]
  wi <- hh[hh$a == hh$b, , drop = FALSE]
  say("    duplicates within cohorts: ", base::nrow(wi), "; duplicates across cohorts: ", base::nrow(wd), "")
  if (base::nrow(wd)) {
    say("    combinations involved in cross-cohort duplicates: ")
    t2 <- base::table(base::paste(wd$a, "pairs", wd$b))
    for (k in base::names(t2)) say("      ", k, ": ", t2[[k]], "")
    say("    The within-disease concordance of these combinations may be inflated by shared individuals and must be recomputed after removing duplicates.")
  }
  for (i in base::seq_len(base::min(20L, base::nrow(hh))))
    say(sprintf("      %-10s %-12s | %-10s %-12s | concordance %.3f",
                hh$a[i], hh$sample_a[i], hh$b[i], hh$sample_b[i], hh$concordance[i]))
} else say("  No duplicate samples meeting the criterion: no individuals are shared between the ten sets and there are no duplicates within cohorts.")

say("")
say("=== Plotting ===")
p <- ggplot(st, aes(median, max, colour = kind)) +
  geom_hline(yintercept = SAME_CUT, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_manual(values = c("within cohort" = PAL[2], "cross-cohort" = PAL[1]),
                      labels = c("Within cohort", "Between cohorts"), name = NULL) +
  labs(title = "Genotype concordance from array SNP probes",
       subtitle = "Each point is one cohort combination; the dashed line marks the threshold for calling two samples the same person",
       x = "Median concordance over all sample pairs",
       y = "Maximum concordance observed") +
  theme_paper()
save_fig(p, "fig51_sample_identity", width = 8, height = 5.5)
log_stop()
}
