# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 44: data quality as a covariate; recompute the five within-disease positive controls
#
# Background: 43 found that sample-level bimodality (proportion of autosomal probes between 0.3 and 0.7)
# is significantly related to case-control status in three cohorts -
#   GSE145361  controls 8.0% vs cases 2.4%    p = 2.0e-07
#   GSE111629  cases 4.8% vs controls 1.3%    p = 0.030
#   GSE87648   controls 10.2% vs cases 2.5%   p = 0.002
# The two Parkinson's disease sets go in opposite directions, so this is a cohort-specific technical phenomenon rather than a property of the disease,
# but it weakens the concordance of effect sizes between the two cohorts.
#
# A covariate rather than removal is chosen for two reasons:
#   (1) removal in these three sets would change the composition of the two groups and create selection bias, which is worse than keeping the samples;
#   (2) a continuous covariate uses the actual value of each sample, needs no threshold and avoids the argument about
#       "how high is bad" - the relative criterion of 43 is suitable for description only, not for selection.
#
# This step recomputes only the five within-disease positive controls and leaves the forty cross-disease comparisons untouched. The positive controls decide
# how the whole paper is read: if the Z of Parkinson's disease rose clearly from 1.81, its earlier weakness would be technical noise,
# the statement "within-disease replication fails for three of five disorders" would have to be rewritten and the forty
# cross-disease comparisons re-run under the new model; if it does not move, the present model is kept and this step is reported as a sensitivity
# analysis. One hour is spent on the decision first, rather than three hours on recomputation.
#
# An additional diagnostic: correlation of the quality index with the methylation principal components already in the model. If the components already
# absorb the quality variation, adding the term naturally changes nothing - this must be shown with data, not guessed.
#
# Environment: Windows + RStudio. About one hour
# Location: <project root>/scripts/44_quality_covariate_sensitivity.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "44_quality_covariate_sensitivity"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
QC_DIR  <- file.path(PROJ_ROOT, "output", "43_qc")
CON_DIR <- file.path(PROJ_ROOT, "output", "39_concordance")
OUT_DIR <- file.path(PROJ_ROOT, "output", "44_quality")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260919)
N_PERM    <- 1000L
BETA_CLIP <- 0.001

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
LAB_EN <- c(parkinson = "Parkinson", schizophrenia = "Schizophrenia",
            alzheimer = "Alzheimer", depression = "Depression", ibd = "IBD")

QS <- read.csv(file.path(QC_DIR, "bimodality_per_sample.csv"), stringsAsFactors = FALSE)
say("Sample quality index of 43 loaded: ", base::nrow(QS), " rows, covering ",
    length(base::unique(QS$gse)), " cohorts")
OLD <- tryCatch(read.csv(file.path(CON_DIR, "concordance_all_pairs.csv"),
                         stringsAsFactors = FALSE), error = function(e) NULL)

# --------------------------------------------------------------
# Per set: prepare the M matrix and the two versions of the design matrix
# --------------------------------------------------------------
prep <- function(g) {
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())

  q <- QS[QS$gse == g, , drop = FALSE]
  fq <- q$frac_mid[base::match(D$gsm, q$gsm)]
  if (base::mean(base::is.na(fq)) > 0.02)
    stop(g, ": quality index missing for ", sprintf("%.1f%%", 100 * base::mean(base::is.na(fq))))
  fq[base::is.na(fq)] <- stats::median(fq, na.rm = TRUE)
  # Standardisation only makes the coefficient readable; it does not affect residuals or tests
  fqz <- (fq - base::mean(fq)) / stats::sd(fq)

  des0 <- D$design
  mpc <- base::grep("^mpc", base::colnames(des0), value = TRUE)
  if (length(mpc)) {
    rr <- base::sapply(mpc, function(k) stats::cor(des0[, k], fqz))
    say("  Correlation of the quality index with the ", length(mpc), " methylation principal components in the model: largest absolute value ",
        sprintf("%.3f", base::max(base::abs(rr))),
        " (", base::names(rr)[base::which.max(base::abs(rr))], ")")
  } else say("  This model contains no methylation principal components")
  gcol <- base::which(base::colnames(des0) == "groupcase")
  des1 <- base::cbind(des0, quality = fqz)
  qq <- qr(des1)
  dropped <- character(0)
  if (qq$rank < base::ncol(des1)) {
    keepc <- base::sort(qq$pivot[base::seq_len(qq$rank)])
    dropped <- base::colnames(des1)[base::setdiff(base::seq_len(base::ncol(des1)), keepc)]
    des1 <- des1[, keepc, drop = FALSE]
  }
  has_q <- "quality" %in% base::colnames(des1)
  say("  After adding the quality covariate: parameters ", base::ncol(des0), " -> ", base::ncol(des1),
      if (length(dropped)) base::paste0("; collinear columns removed ",
        base::paste(utils::head(dropped, 4), collapse = ", ")) else "",
      if (!has_q) "   [NOTE] the quality term itself was found collinear and removed, i.e. it is already covered by the existing covariates" else "")
  list(M = M, des0 = des0, des1 = des1, gsm = D$gsm, group = D$group,
       probes = D$probes, has_q = has_q)
}

eff_of <- function(M, des) {
  fit <- limma::eBayes(limma::lmFit(M, des))
  base::as.numeric(limma::topTable(fit, coef = "groupcase", number = Inf,
                                   sort.by = "none")$logFC)
}
state_of <- function(M, des) {
  jg <- base::which(base::colnames(des) == "groupcase")
  X0 <- des[, -jg, drop = FALSE]
  q0 <- qr(X0)
  list(q0 = q0, R = qr.resid(q0, base::t(M)))
}
coef_fwl <- function(S, gvec) {
  gr <- qr.resid(S$q0, gvec)
  base::as.numeric(base::crossprod(gr, S$R)) / base::sum(gr^2)
}

rows <- list()
for (d in names(DIS)) {
  a <- DIS[[d]][1]; bcoh <- DIS[[d]][2]
  say("==== ", LAB_EN[[d]], ": ", a, " vs ", bcoh, " ====")

  say("  [fixed side] ", a)
  PA <- prep(a)
  ea0 <- eff_of(PA$M, PA$des0); ea1 <- eff_of(PA$M, PA$des1)
  say("  Correlation of effect sizes between the two models: r = ", sprintf("%.4f", stats::cor(ea0, ea1)))
  pa <- PA$probes; has_qa <- PA$has_q
  rm(PA); invisible(gc())

  say("  [permuted side] ", bcoh)
  PB <- prep(bcoh)
  gnum <- base::as.numeric(PB$group == "case")
  S0 <- state_of(PB$M, PB$des0); S1 <- state_of(PB$M, PB$des1)
  eb0 <- coef_fwl(S0, gnum); eb1 <- coef_fwl(S1, gnum)
  say("  Correlation of effect sizes between the two models: r = ", sprintf("%.4f", stats::cor(eb0, eb1)))
  pb <- PB$probes; has_qb <- PB$has_q
  rm(PB); invisible(gc())

  cm <- base::intersect(pa, pb)
  ia <- base::match(cm, pa); ib <- base::match(cm, pb)
  for (v in c("original model", "with quality covariate")) {
    ea <- if (v == "original model") ea0[ia] else ea1[ia]
    eb <- if (v == "original model") eb0[ib] else eb1[ib]
    S  <- if (v == "original model") S0 else S1
    rho <- stats::cor(ea, eb, method = "spearman")
    dir <- base::mean(base::sign(ea) == base::sign(eb))
    ra <- base::rank(ea)
    pr <- numeric(N_PERM)
    for (i in base::seq_len(N_PERM))
      pr[i] <- stats::cor(ra, base::rank(coef_fwl(S, base::sample(gnum))[ib]))
    pv <- (1 + base::sum(pr >= rho)) / (N_PERM + 1)
    z <- (rho - base::mean(pr)) / stats::sd(pr)
    say(sprintf("  %-14s probes %6d  rho=%+.4f  sign %.4f  Z=%+5.2f  P=%.4f",
                v, length(cm), rho, dir, z, pv))
    rows[[length(rows) + 1]] <- data.frame(
      disease = LAB_EN[[d]], model = v, n_probe = length(cm), rho = rho,
      dir = dir, z = z, p = pv, quality_kept_fixed = has_qa,
      quality_kept_permuted = has_qb, stringsAsFactors = FALSE)
  }
  rm(S0, S1); invisible(gc())
  say("")
}

tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "within_disease_quality_sensitivity.csv"), row.names = FALSE)
say("=== Five within-disease positive controls: the two models side by side ===")
for (d in base::unique(tab$disease)) {
  s <- tab[tab$disease == d, ]
  say(sprintf("  %-14s original rho=%+.4f Z=%+5.2f P=%.4f | with quality rho=%+.4f Z=%+5.2f P=%.4f",
              d, s$rho[1], s$z[1], s$p[1], s$rho[2], s$z[2], s$p[2]))
}
say("")
say("  Original model reaching P<0.05: ", base::sum(tab$p[tab$model == "original model"] < 0.05), " / 5; ",
    "with quality covariate: ", base::sum(tab$p[tab$model == "with quality covariate"] < 0.05), " / 5")
dz <- tab$z[tab$model == "with quality covariate"] - tab$z[tab$model == "original model"]
say("  Change in Z: ", base::paste(sprintf("%s %+.2f", base::unique(tab$disease), dz),
                                collapse = "; "))
say("")
if (base::max(base::abs(dz)) < 0.5) {
  say("  Interpretation: the change in Z is below 0.5 for all five positive controls; quality confounding does not materially affect the conclusions.")
  say("  The present model is kept and this step is reported as a sensitivity analysis; the forty cross-disease comparisons need not be re-run.")
} else {
  say("  Interpretation: the change in Z reaches or exceeds 0.5 for at least one positive control; quality confounding affects the conclusions.")
  say("  The ten EWAS need to be rebuilt under the model with the quality covariate, and the forty cross-disease comparisons and the composition matrix re-run.")
}

say("")
say("=== Plotting ===")
tab$model <- factor(tab$model, levels = c("original model", "with quality covariate"))
p <- ggplot(tab, aes(disease, z, fill = model)) +
  geom_hline(yintercept = stats::qnorm(0.95), linetype = 2, colour = "grey55",
             linewidth = 0.4) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  scale_fill_manual(values = c("original model" = PAL[1], "with quality covariate" = PAL[3]),
                    labels = c("Original model", "Plus quality covariate"), name = NULL) +
  labs(title = "Within-disease replication before and after adjusting for sample quality",
       subtitle = "Dashed line marks the one-sided 5% threshold of the permutation null",
       x = NULL, y = "Z against the permutation null") +
  theme_paper()
save_fig(p, "fig47_quality_sensitivity", width = 8, height = 5)
log_stop()
