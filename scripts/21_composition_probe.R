# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 21: interrogate the positive result of axis two
#
# Result of 20: PD and IBD agree in direction for all twelve cell types, permutation P = 0.011;
# both negative controls are exactly 6/12, i.e. chance level. This is a positive result that axis one (CpG effect sizes)
# did not give, so it must be questioned as hard as possible before deciding whether to believe it.
#
# Three questions must be answered:
#
#   1   Is this merely "neutrophils up, lymphocytes down"?
#       In the six data sets the pattern of PD and IBD is almost always neutrophils up and the rest down,
#       which is essentially a raised neutrophil/lymphocyte ratio, present in any systemic inflammatory state.
#       If the concordance disappears after removing neutrophils, the sharing is non-specific
#       inflammatory background and the wording must be tightened accordingly.
#
#   2   Is it caused by medication? Glucocorticoids cause exactly neutrophilia with fewer lymphocytes and
#       eosinophils, matching the observed pattern. The same applies to levodopa in PD.
#       GSE72774 carries the daily levodopa dose, which was kept for this step.
#
#   3   Does it hold with another reference matrix? The twelve fractions are estimated by EpiDISH from the same
#       methylation data, and the estimation error may be related to disease. Re-estimating with the seven-cell-type reference
#       and finding a reversed conclusion would mean that the result depends on the choice of reference matrix.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/21_composition_probe.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "21_composition_probe"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RDS_DIR <- file.path(DATA_ROOT, "rds")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "21_probe")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 10000L
RUN_REF7 <- TRUE          # Re-estimation with the seven-cell-type reference requires reading the beta matrices again and is slow
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

ROLE <- c(GSE145361 = "PD discovery", GSE111629 = "PD replication",
          GSE72774 = "PD sensitivity", GSE87648 = "IBD main",
          GSE112611 = "IBD paediatric", GSE152027 = "Negative control")
AGE_SOURCE <- c(GSE145361 = "dnam", GSE152027 = "dnam", GSE111629 = "chrono",
                GSE72774 = "chrono", GSE87648 = "chrono", GSE112611 = "chrono")
USE_SMOKING <- c(GSE145361 = TRUE, GSE152027 = TRUE, GSE111629 = TRUE,
                 GSE72774 = TRUE, GSE87648 = TRUE, GSE112611 = FALSE)
PAIRS <- list(
  list(id = "PC1", kind = "positive control", a = "GSE145361", b = "GSE111629"),
  list(id = "PC2", kind = "positive control", a = "GSE87648",  b = "GSE112611"),
  list(id = "M1",  kind = "primary criterion",   a = "GSE145361", b = "GSE87648"),
  list(id = "M2",  kind = "secondary",     a = "GSE111629", b = "GSE87648"),
  list(id = "NC1", kind = "negative control", a = "GSE145361", b = "GSE152027"),
  list(id = "NC2", kind = "negative control", a = "GSE87648",  b = "GSE152027"))

lg <- function(p) { p <- pmin(pmax(p, 1e-4), 1 - 1e-4); log(p / (1 - p)) }
slide_of <- function(g, gsm) {
  mf <- file.path(MAP_DIR, paste0(g, "_sample_map.csv"))
  if (!file.exists(mf)) return(rep(NA_character_, length(gsm)))
  m <- read.csv(mf, stringsAsFactors = FALSE, colClasses = "character")
  bc <- rep(NA_character_, nrow(m))
  hit <- regexpr(BARCODE, m$col_name) > 0
  bc[hit] <- regmatches(m$col_name, regexpr(BARCODE, m$col_name))
  sub("_R[0-9]{2}C[0-9]{2}$", "", bc)[match(gsm, m$gsm)]
}

# Build the design matrix and cell fractions of one cohort
build <- function(g, prop = NULL, extra = NULL) {
  f <- file.path(SET_DIR, paste0(g, "_analysis_filtered.rds"))
  obj <- readRDS(f); ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(colnames(obj$beta), cov$gsm), , drop = FALSE]
  d <- data.frame(group = factor(as.character(ph$group),
                                 levels = c("control", "case")))
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sx)) d$sex <- factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  if (AGE_SOURCE[g] == "dnam") d$age <- cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    d$age <- suppressWarnings(as.numeric(ph[[ag]]))
  }
  if (isTRUE(USE_SMOKING[g])) d$smoking <- cov$dnam_smoking
  if ("scan_date" %in% names(ph)) d$batch <- factor(as.character(ph$scan_date))
  else {
    sl <- slide_of(g, colnames(obj$beta))
    if (sum(!is.na(sl)) > 0.95 * length(sl) &&
        length(unique(sl)) <= 0.25 * length(sl)) d$batch <- factor(sl)
  }
  if (!is.null(extra)) for (k in names(extra)) d[[k]] <- extra[[k]]
  if (is.null(prop)) {
    ct <- grep("^ct_", names(cov), value = TRUE)
    Y <- as.matrix(cov[, ct, drop = FALSE]); colnames(Y) <- sub("^ct_", "", ct)
  } else Y <- prop
  ok <- stats::complete.cases(d) & stats::complete.cases(Y)
  list(d = d[ok, , drop = FALSE], Y = Y[ok, , drop = FALSE],
       pheno = ph[ok, , drop = FALSE], gsm = colnames(obj$beta)[ok])
}
coef_vec <- function(B, perm = FALSE) {
  d <- B$d
  if (perm) d$group <- sample(d$group)
  X <- stats::model.matrix(~ ., data = d)
  j <- which(colnames(X) == "groupcase")
  out <- vapply(seq_len(ncol(B$Y)), function(k)
    stats::lm.fit(X, lg(B$Y[, k]))$coefficients[j], numeric(1))
  names(out) <- colnames(B$Y); out
}
test_pair <- function(Ba, Bb, ea, eb, keep) {
  a <- ea[keep]; b <- eb[keep]
  rho <- stats::cor(a, b, method = "spearman")
  dir <- mean(sign(a) == sign(b))
  pr <- numeric(N_PERM); pd <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    bp <- coef_vec(Bb, perm = TRUE)[keep]
    pr[i] <- stats::cor(a, bp, method = "spearman")
    pd[i] <- mean(sign(a) == sign(bp))
  }
  c(rho = rho, dir = dir,
    p_rho = (1 + sum(pr >= rho)) / (N_PERM + 1),
    p_dir = (1 + sum(pd >= dir)) / (N_PERM + 1),
    z = (rho - mean(pr)) / stats::sd(pr))
}

# --------------------------------------------------------------
say("=== 1. Is the sharing merely neutrophils up, lymphocytes down? ===")
B <- list(); E <- list()
for (g in names(ROLE)) { B[[g]] <- build(g); E[[g]] <- coef_vec(B[[g]]) }
cells <- names(E[[1]])
say("Cell types: ", paste(cells, collapse = ", "))
say("Neutrophil effect in each cohort: ",
    paste(sprintf("%s %+.2f", names(E), vapply(E, function(x) x[["Neu"]], 0)),
          collapse = "; "))
say("")

SETS <- list(
  "all twelve types"      = cells,
  "without neutrophils"      = setdiff(cells, "Neu"),
  "without neutrophils and eosinophils" = setdiff(cells, c("Neu", "Eos")),
  "lymphocyte subsets only"  = intersect(cells, c("CD4Tnv", "CD4Tmem", "CD8Tnv",
                                         "CD8Tmem", "Treg", "Bnv", "Bmem", "NK")))
rows <- list()
for (P in PAIRS) for (nm in names(SETS)) {
  keep <- SETS[[nm]]
  r <- test_pair(B[[P$a]], B[[P$b]], E[[P$a]], E[[P$b]], keep)
  rows[[length(rows) + 1]] <- data.frame(
    pair = P$id, kind = P$kind, subset = nm, n_ct = length(keep),
    rho = r[["rho"]], dir = r[["dir"]], p_rho = r[["p_rho"]],
    p_dir = r[["p_dir"]], z = r[["z"]], stringsAsFactors = FALSE)
}
d1 <- do.call(rbind, rows)
write.csv(d1, file.path(OUT_DIR, "subset_tests.csv"), row.names = FALSE)
for (P in PAIRS) {
  s <- d1[d1$pair == P$id, ]
  say("  ", P$id, " ", P$kind, "  ", P$a, " vs ", P$b)
  for (i in seq_len(nrow(s)))
    say(sprintf("     %-18s n=%2d  rho=%+.3f  sign %.2f  P(sign)=%.4f",
                s$subset[i], s$n_ct[i], s$rho[i], s$dir[i], s$p_dir[i]))
}
say("")

# --------------------------------------------------------------
say("=== 2. Levodopa sensitivity (GSE72774) ===")
g <- "GSE72774"
ph <- readRDS(file.path(SET_DIR, paste0(g, "_analysis_filtered.rds")))$pheno
ld <- suppressWarnings(as.numeric(ph$levodopaperdayincludingagonists))
if (all(is.na(ld))) say("  This cohort has no usable levodopa dose field") else {
  say("  Daily levodopa dose: median ", stats::median(ld, na.rm = TRUE),
      ", maximum ", max(ld, na.rm = TRUE), ", non-zero ", sum(ld > 0, na.rm = TRUE), " samples")
  B0 <- build(g)
  B1 <- build(g, extra = list(levodopa = ld))
  e0 <- coef_vec(B0); e1 <- coef_vec(B1)
  say("  Effects per cell type before and after adjustment: ")
  for (k in names(e0))
    say(sprintf("     %-9s before %+.3f   after %+.3f   change %+.3f",
                k, e0[[k]], e1[[k]], e1[[k]] - e0[[k]]))
  say("  Correlation between the two effect vectors = ", sprintf("%.3f", stats::cor(e0, e1)))
  r <- test_pair(B[["GSE87648"]], B1, E[["GSE87648"]], e1, cells)
  say("  After adjustment for levodopa, concordance of this cohort with the IBD main set: rho = ",
      sprintf("%.3f", r[["rho"]]), ", sign ", sprintf("%.2f", r[["dir"]]),
      ", P(sign) = ", format.pval(r[["p_dir"]], digits = 3))
}
say("")

# --------------------------------------------------------------
say("=== 3. Re-estimation with the seven-cell-type reference matrix ===")
if (RUN_REF7) {
  suppressPackageStartupMessages(library(EpiDISH))
  e <- new.env(); utils::data("centDHSbloodDMC.m", package = "EpiDISH", envir = e)
  ref7 <- get("centDHSbloodDMC.m", envir = e)
  say("  Reference matrix ", nrow(ref7), " probes x ", ncol(ref7), " types: ",
      paste(colnames(ref7), collapse = ", "))
  B7 <- list(); E7 <- list()
  for (g in names(ROLE)) {
    bf <- file.path(RDS_DIR, paste0(g, "_beta.rds"))
    mf <- file.path(MAP_DIR, paste0(g, "_sample_map.csv"))
    if (!file.exists(bf) || !file.exists(mf)) { say("  ", g, ": raw matrix missing, skipped"); next }
    m <- read.csv(mf, stringsAsFactors = FALSE, colClasses = "character")
    bm <- readRDS(bf)
    bm <- bm[, as.integer(m$col_pos_in_data) - 1L, drop = FALSE]
    colnames(bm) <- m$gsm
    est <- EpiDISH::epidish(beta.m = bm, ref.m = ref7, method = "RPC")$estF
    rm(bm); invisible(gc())
    B7[[g]] <- build(g, prop = est[match(build(g)$gsm, rownames(est)), , drop = FALSE])
    E7[[g]] <- coef_vec(B7[[g]])
    say("  ", g, " finished, median row sum ", sprintf("%.3f", stats::median(rowSums(est))))
  }
  c7 <- names(E7[[1]])
  rows <- list()
  for (P in PAIRS) {
    if (is.null(E7[[P$a]]) || is.null(E7[[P$b]])) next
    r <- test_pair(B7[[P$a]], B7[[P$b]], E7[[P$a]], E7[[P$b]], c7)
    rows[[length(rows) + 1]] <- data.frame(
      pair = P$id, kind = P$kind, rho = r[["rho"]], dir = r[["dir"]],
      p_rho = r[["p_rho"]], p_dir = r[["p_dir"]], stringsAsFactors = FALSE)
    say(sprintf("  %-4s %-8s rho=%+.3f  sign %.2f  P(sign)=%.4f",
                P$id, P$kind, r[["rho"]], r[["dir"]], r[["p_dir"]]))
  }
  d3 <- do.call(rbind, rows)
  write.csv(d3, file.path(OUT_DIR, "ref7_tests.csv"), row.names = FALSE)
} else say("  Skipped (RUN_REF7 = FALSE)")

# --------------------------------------------------------------
say("")
say("=== Plotting ===")
d1$panel <- factor(paste0(d1$pair, " ", d1$kind), levels = unique(paste0(d1$pair, " ", d1$kind)))
d1$subset <- factor(d1$subset, levels = names(SETS))
p <- ggplot(d1, aes(subset, dir, fill = kind)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("P=%.3f", p_dir)), vjust = -0.4, size = 2.3) +
  facet_wrap(~ panel, ncol = 3) +
  scale_fill_manual(values = c("positive control" = PAL[3], "primary criterion" = PAL[2],
                               "secondary" = PAL[5], "negative control" = PAL[8]), name = NULL) +
  labs(title = "Is the shared composition signal more than a neutrophil shift?",
       subtitle = "Direction concordance after removing neutrophils and eosinophils from the comparison",
       x = NULL, y = "Proportion of cell types concordant in direction") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
save_fig(p, "fig27_composition_subsets", width = 10, height = 6)

log_stop()
