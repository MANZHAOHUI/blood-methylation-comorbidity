# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 16: investigate the residual inflation in GSE145361 and GSE152027
#
# Conclusion of 15: with cell fractions, four of the six sets reach lambda between 0.95 and 1.10,
# but GSE145361 reaches 1.385 at best (K=30) and GSE152027 1.187 at best (K=50),
# and further principal components bring no further decrease. All VIFs are below 3.2, so collinearity is not the cause.
#
# This step tests four sources not yet excluded:
#   (1) sample-wide distributions are not aligned. If the deposited matrix was not normalised between samples,
#       the overall M value level differs between samples and shifts all probes in the same direction,
#       producing very strong correlation structure. This is the most common cause of large lambda inflation
#       and the most easily overlooked, because each data set looks "normal" on its own.
#   (2) DNAm age as a covariate. It is a linear combination computed from the same methylation data,
#       which amounts to putting a direction of the data themselves into the model.
#   (3) array as a fixed effect with 168 levels, which may add noise rather than remove it.
#   (4) a few extreme samples.
#
# Sample-level descriptions are produced for all six sets first (fast), then variants are compared for the two problematic sets.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/16_residual_inflation.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "16_residual_inflation"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "16_inflation")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260914)
N_PROBE <- 50000L
BETA_CLIP <- 0.001
DROP_CT <- "ct_Baso"
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
PROBLEM <- c("GSE145361", "GSE152027")
K_BASE  <- c(GSE145361 = 30L, GSE152027 = 50L)

ROLE <- c(GSE145361 = "PD discovery", GSE111629 = "PD replication",
          GSE72774 = "PD sensitivity", GSE87648 = "IBD main",
          GSE112611 = "IBD paediatric", GSE152027 = "Negative control")
ORDER <- names(ROLE)
AGE_SOURCE <- c(GSE145361 = "dnam", GSE152027 = "dnam", GSE111629 = "chrono",
                GSE72774 = "chrono", GSE87648 = "chrono", GSE112611 = "chrono")
USE_SMOKING <- c(GSE145361 = TRUE, GSE152027 = TRUE, GSE111629 = TRUE,
                 GSE72774 = TRUE, GSE87648 = TRUE, GSE112611 = FALSE)

slide_of <- function(g, gsm) {
  mf <- file.path(MAP_DIR, paste0(g, "_sample_map.csv"))
  if (!file.exists(mf)) return(rep(NA_character_, length(gsm)))
  m <- read.csv(mf, stringsAsFactors = FALSE, colClasses = "character")
  bc <- rep(NA_character_, nrow(m))
  hit <- regexpr(BARCODE, m$col_name) > 0
  bc[hit] <- regmatches(m$col_name, regexpr(BARCODE, m$col_name))
  sub("_R[0-9]{2}C[0-9]{2}$", "", bc)[match(gsm, m$gsm)]
}
lambda_of <- function(p) {
  cs <- stats::qchisq(p, df = 1, lower.tail = FALSE)
  stats::median(cs, na.rm = TRUE) / stats::qchisq(0.5, 1)
}

files <- list.files(SET_DIR, pattern = "_analysis_filtered\\.rds$", full.names = TRUE)
samp_stats <- list(); variants <- list()

for (f in files) {
  obj <- readRDS(f); g <- obj$gse
  say("==== ", g, " (", ROLE[g], ") ====")
  ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(colnames(obj$beta), cov$gsm), , drop = FALSE]
  grp <- factor(as.character(ph$group), levels = c("control", "case"))

  idx <- sample(seq_len(nrow(obj$beta)), min(N_PROBE, nrow(obj$beta)))
  b <- obj$beta[idx, , drop = FALSE]
  b[!is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- log2(b / (1 - b)); M <- M[stats::complete.cases(M), , drop = FALSE]
  rm(b); invisible(gc())

  # ---------- sample-level description ----------
  cm <- colMeans(M); cs <- apply(M, 2, stats::sd)
  pm <- stats::t.test(cm ~ grp)$p.value
  ps <- stats::t.test(cs ~ grp)$p.value
  say("  Sample mean M: overall median ", sprintf("%.4f", median(cm)),
      ", difference between groups ", sprintf("%.4f", diff(tapply(cm, grp, mean))),
      ", p = ", format.pval(pm, digits = 2))
  say("  Sample SD of M: overall median ", sprintf("%.4f", median(cs)),
      ", difference between groups ", sprintf("%.4f", diff(tapply(cs, grp, mean))),
      ", p = ", format.pval(ps, digits = 2))
  say("  Dispersion of sample mean M (SD / absolute mean): ",
      sprintf("%.4f", stats::sd(cm) / abs(mean(cm))))
  samp_stats[[g]] <- data.frame(gse = g, gsm = colnames(M),
                                group = as.character(grp),
                                mean_M = cm, sd_M = cs,
                                stringsAsFactors = FALSE)

  if (!(g %in% PROBLEM)) { say(""); rm(obj, M); invisible(gc()); next }

  # ---------- comparison of variants (only for the two problematic sets) ----------
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  sex <- if (is.na(sx)) NULL else factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  age <- if (AGE_SOURCE[g] == "dnam") cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    suppressWarnings(as.numeric(ph[[ag]]))
  }
  cts <- cov[, setdiff(grep("^ct_", names(cov), value = TRUE), DROP_CT), drop = FALSE]
  smk <- if (isTRUE(USE_SMOKING[g])) cov$dnam_smoking else NULL
  bvar <- if ("scan_date" %in% names(ph)) as.character(ph$scan_date) else {
    sl <- slide_of(g, colnames(obj$beta))
    if (sum(!is.na(sl)) > 0.95 * length(sl)) sl else NULL
  }

  build <- function(Mx, use_age = TRUE, use_batch = TRUE, add_mean = FALSE,
                    keep_idx = NULL, K = K_BASE[[g]]) {
    if (is.null(keep_idx)) keep_idx <- seq_len(ncol(Mx))
    Mx <- Mx[, keep_idx, drop = FALSE]
    gg <- droplevels(grp[keep_idx])
    d <- data.frame(group = gg)
    if (!is.null(sex)) d$sex <- droplevels(sex[keep_idx])
    if (use_age) d$age <- age[keep_idx]
    if (!is.null(smk)) d$smoking <- smk[keep_idx]
    if (use_batch && !is.null(bvar)) d$batch <- droplevels(factor(bvar[keep_idx]))
    d <- cbind(d, cts[keep_idx, , drop = FALSE])
    if (add_mean) d$meanM <- colMeans(Mx)
    if (K > 0) {
      pr <- stats::prcomp(t(Mx), center = TRUE, scale. = FALSE)
      mp <- pr$x[, seq_len(min(K + 15L, ncol(pr$x))), drop = FALSE]
      pg <- apply(mp, 2, function(v) stats::anova(stats::lm(v ~ gg))$`Pr(>F)`[1])
      take <- utils::head(which(pg >= 0.01), K)
      if (length(take)) {
        pcd <- as.data.frame(mp[, take, drop = FALSE])
        names(pcd) <- paste0("mpc", seq_along(take))
        d <- cbind(d, pcd)
      }
    }
    ok <- stats::complete.cases(d)
    des <- stats::model.matrix(~ ., data = d[ok, , drop = FALSE])
    if (qr(des)$rank < ncol(des)) return(NA_real_)
    fit <- limma::eBayes(limma::lmFit(Mx[, ok, drop = FALSE], des))
    lambda_of(fit$p.value[, "groupcase"])
  }

  # Between-sample quantile normalisation used by variant 2
  Mq <- limma::normalizeQuantiles(M)

  # Variant 5: samples deviating by four standard deviations on the first five principal components
  pr0 <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  z <- scale(pr0$x[, 1:5, drop = FALSE])
  keep_in <- which(apply(abs(z), 1, max) <= 4)
  say("  Remaining after removal of extreme samples: ", length(keep_in), " / ", ncol(M))

  vs <- list(
    V0 = list(desc = "baseline (with cells, K as optimal in 15)", f = function() build(M)),
    V1 = list(desc = "sample mean M added as covariate", f = function() build(M, add_mean = TRUE)),
    V2 = list(desc = "after between-sample quantile normalisation", f = function() build(Mq)),
    V3 = list(desc = "age term removed", f = function() build(M, use_age = FALSE)),
    V4 = list(desc = "batch term removed", f = function() build(M, use_batch = FALSE)),
    V5 = list(desc = "extreme samples removed", f = function() build(M, keep_idx = keep_in)))
  for (nm in names(vs)) {
    lam <- tryCatch(vs[[nm]]$f(), error = function(e) NA_real_)
    say("  ", nm, "  ", vs[[nm]]$desc, "  lambda = ",
        if (is.na(lam)) "failed" else sprintf("%.3f", lam))
    variants[[length(variants) + 1]] <- data.frame(
      gse = g, variant = nm, desc = vs[[nm]]$desc, lambda = lam,
      stringsAsFactors = FALSE)
  }
  say("")
  rm(obj, M, Mq, pr0); invisible(gc())
}

ss <- do.call(rbind, samp_stats)
write.csv(ss, file.path(OUT_DIR, "sample_level_stats.csv"), row.names = FALSE)
if (length(variants)) {
  vv <- do.call(rbind, variants)
  write.csv(vv, file.path(OUT_DIR, "inflation_variants.csv"), row.names = FALSE)
  say("=== Summary of variants ===")
  print(vv[, c("gse", "variant", "desc", "lambda")])
}

say("")
say("=== Plotting ===")
ss$panel <- factor(paste0(ss$gse, "\n", ROLE[ss$gse]),
                   levels = paste0(ORDER, "\n", ROLE[ORDER]))
ss$group <- factor(ss$group, levels = c("control", "case"))
p <- ggplot(ss, aes(mean_M, sd_M, colour = group)) +
  geom_point_big(size = 0.7, alpha = 0.6) +
  facet_wrap(~ panel, ncol = 3, scales = "free") +
  scale_colour_manual(values = c(control = PAL[1], case = PAL[2]),
                      name = NULL, labels = c("Control", "Case")) +
  labs(title = "Sample-level distribution of M-values",
       subtitle = "Samples spread apart on these axes indicate the deposited matrices were not normalised across samples",
       x = "Mean M-value per sample", y = "SD of M-values per sample") +
  theme_paper() +
  guides(colour = guide_legend(override.aes = list(size = 2)))
save_fig(p, "fig19_sample_level_M", width = 10, height = 6)

if (length(variants)) {
  vv$panel <- factor(paste0(vv$gse, "\n", ROLE[vv$gse]),
                     levels = paste0(ORDER, "\n", ROLE[ORDER]))
  p <- ggplot(vv, aes(variant, lambda)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
    geom_col(fill = PAL[1], width = 0.65) +
    geom_text(aes(label = sprintf("%.2f", lambda)), vjust = -0.5, size = 2.8) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    labs(title = "Residual inflation under model variants",
         subtitle = "V0 baseline; V1 +mean M; V2 quantile-normalised; V3 no age; V4 no batch; V5 outliers removed",
         x = NULL, y = expression(lambda)) +
    theme_paper()
  save_fig(p, "fig20_inflation_variants", width = 9, height = 4.4)
}

log_stop()
