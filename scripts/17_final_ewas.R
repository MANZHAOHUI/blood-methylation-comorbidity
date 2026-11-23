# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 17: formal within-cohort EWAS (frozen settings of decisions D4, D12 and D13)
#
# Differences from 13:
#   - methylation principal components are computed on the full filtered matrix, no longer on a subset of fifty thousand probes
#   - each data set uses its own K according to Table 4
#   - genomic control is applied to cohorts with an inflation factor above 1.10; effect sizes are not rescaled
#   - the design matrices are saved for reuse in the permutation tests of 18, avoiding recomputation of the principal components
#
# Note on the principal component algorithm: principal components in the sample direction are needed. A direct SVD
# of a matrix of 400 thousand by more than a thousand is wasteful; instead the Gram matrix between samples is computed
# (inner products over the probe dimension, only n x n) and then eigen-decomposed.
# This is mathematically equivalent to an SVD of the centred matrix but saves most of the computation.
#
# Environment: Windows + RStudio (peak memory about 12 GB, for the largest cohort)
# Location: <project root>/scripts/17_final_ewas.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "17_final_ewas"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "17_ewas")
for (d in c(RES_DIR, OUT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260914)
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
MPC_GROUP_P <- 0.01
GC_CUT      <- 1.10          # Genomic control is applied above this value
BARCODE     <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

# Frozen settings of Table 4
K_BY <- c(GSE145361 = 30L, GSE111629 = 0L, GSE72774 = 10L,
          GSE87648 = 5L, GSE112611 = 0L, GSE152027 = 50L)
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

# Principal components in the sample direction: Gram matrix + eigen-decomposition
sample_pcs <- function(M, npc) {
  Mc <- M - rowMeans(M, na.rm = TRUE)
  G <- crossprod(Mc)                    # n x n
  rm(Mc); invisible(gc())
  ev <- eigen(G, symmetric = TRUE)
  k <- min(npc, sum(ev$values > 1e-8))
  sc <- ev$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(pmax(ev$values[seq_len(k)], 0)), nrow = k)
  rownames(sc) <- colnames(M)
  colnames(sc) <- paste0("PC", seq_len(k))
  list(scores = sc, var = ev$values / sum(ev$values))
}

files <- list.files(SET_DIR, pattern = "_analysis_filtered\\.rds$", full.names = TRUE)
summ <- list(); qqd <- list()

for (f in files) {
  obj <- readRDS(f); g <- obj$gse
  K <- K_BY[[g]]
  say("==== ", g, " (", ROLE[g], ")  K = ", K, " ====")
  ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(colnames(obj$beta), cov$gsm), , drop = FALSE]
  if (anyNA(cov$gsm)) stop(g, ": covariates and matrix samples do not match")

  b <- obj$beta
  b[!is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- log2(b / (1 - b))
  rm(b); invisible(gc())
  keep_probe <- stats::complete.cases(M)
  M <- M[keep_probe, , drop = FALSE]
  say("  Probes ", nrow(M), " (removed for missing values: ", sum(!keep_probe), "), samples ", ncol(M))

  grp <- factor(as.character(ph$group), levels = c("control", "case"))

  d <- data.frame(group = grp)
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sx)) d$sex <- factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  if (AGE_SOURCE[g] == "dnam") d$age <- cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    d$age <- suppressWarnings(as.numeric(ph[[ag]]))
  }
  if (isTRUE(USE_SMOKING[g])) d$smoking <- cov$dnam_smoking
  bv <- if ("scan_date" %in% names(ph)) as.character(ph$scan_date) else {
    sl <- slide_of(g, colnames(M))
    if (sum(!is.na(sl)) > 0.95 * length(sl)) sl else NULL
  }
  if (!is.null(bv)) d$batch <- factor(bv)
  cts <- cov[, setdiff(grep("^ct_", names(cov), value = TRUE), DROP_CT), drop = FALSE]
  d <- cbind(d, cts)

  if (K > 0) {
    say("  Computing sample principal components on the full matrix ...")
    t0 <- Sys.time()
    pcs <- sample_pcs(M, K + 25L)
    say("  Time: ", sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))),
        " min; the first ", ncol(pcs$scores), " components explain cumulatively ",
        sprintf("%.1f%%", 100 * sum(pcs$var[seq_len(ncol(pcs$scores))])))
    pg <- apply(pcs$scores, 2, function(v)
      stats::anova(stats::lm(v ~ grp))$`Pr(>F)`[1])
    avail <- which(pg >= MPC_GROUP_P)
    take <- utils::head(avail, K)
    say("  Components associated with group: ", sum(pg < MPC_GROUP_P), " excluded; used: ",
        length(take), " (indices ", paste(range(take), collapse = " to "), ")")
    if (length(take) < K) say("  [NOTE] fewer usable components than K; actually used: ", length(take), "")
    pcd <- as.data.frame(pcs$scores[, take, drop = FALSE])
    names(pcd) <- paste0("mpc", seq_along(take))
    d <- cbind(d, pcd)
    rm(pcs); invisible(gc())
  }

  ok <- stats::complete.cases(d)
  if (sum(!ok)) say("  Removed because of missing covariates: ", sum(!ok), " samples")
  d <- d[ok, , drop = FALSE]; M <- M[, ok, drop = FALSE]
  des <- stats::model.matrix(~ ., data = d)
  if (qr(des)$rank < ncol(des)) stop(g, ": design matrix rank-deficient")
  tb <- table(d$group)
  say("  In the model: controls ", tb[["control"]], ", cases ", tb[["case"]],
      "; parameters ", ncol(des), ", residual degrees of freedom ", nrow(des) - ncol(des))

  fit <- limma::eBayes(limma::lmFit(M, des))
  tt <- limma::topTable(fit, coef = "groupcase", number = Inf, sort.by = "none")
  res <- data.frame(probe = rownames(tt), effect = tt$logFC,
                    se = fit$stdev.unscaled[, "groupcase"] * sqrt(fit$s2.post),
                    t = tt$t, p = tt$P.Value, stringsAsFactors = FALSE)
  lam <- lambda_of(res$p)
  if (lam > GC_CUT) {
    cs <- stats::qchisq(res$p, 1, lower.tail = FALSE) / lam
    res$p_adj_scale <- stats::pchisq(cs, 1, lower.tail = FALSE)
    gc_flag <- TRUE
  } else {
    res$p_adj_scale <- res$p; gc_flag <- FALSE
  }
  res$fdr <- stats::p.adjust(res$p_adj_scale, method = "BH")
  say("  lambda = ", sprintf("%.3f", lam),
      "; genomic control: ", if (gc_flag) "applied" else "not applied",
      "; probes with FDR < 0.05: ", sum(res$fdr < 0.05, na.rm = TRUE))

  saveRDS(res, file.path(RES_DIR, paste0(g, "_ewas_final.rds")))
  saveRDS(list(design = des, gsm = colnames(M), probes = res$probe,
               group = d$group, gse = g),
          file.path(RES_DIR, paste0(g, "_design.rds")))
  write.csv(utils::head(res[order(res$p), ], 200),
            file.path(OUT_DIR, paste0(g, "_top200.csv")), row.names = FALSE)

  summ[[g]] <- data.frame(gse = g, role = ROLE[g], disease = obj$disease,
                          K = K, n_control = tb[["control"]], n_case = tb[["case"]],
                          n_probe = nrow(res), lambda = lam, gc = gc_flag,
                          n_fdr05 = sum(res$fdr < 0.05, na.rm = TRUE),
                          stringsAsFactors = FALSE)
  o <- order(res$p)
  qqd[[g]] <- data.frame(gse = g, lambda = lam,
                         expected = -log10(stats::ppoints(nrow(res))),
                         observed = -log10(res$p[o]), stringsAsFactors = FALSE)
  say("")
  rm(obj, M, fit, des); invisible(gc())
}

s_all <- do.call(rbind, summ)
write.csv(s_all, file.path(OUT_DIR, "ewas_final_summary.csv"), row.names = FALSE)
say("=== Summary ===")
print(s_all[, c("gse", "role", "K", "n_control", "n_case", "lambda", "gc", "n_fdr05")])

say("")
say("=== Plotting ===")
qd <- do.call(rbind, qqd)
qd$panel <- factor(paste0(qd$gse, "\n", ROLE[qd$gse]),
                   levels = paste0(ORDER, "\n", ROLE[ORDER]))
lt <- unique(qd[, c("panel", "lambda")])
lt$lab <- paste0("lambda = ", sprintf("%.3f", lt$lambda))
p <- ggplot(qd, aes(expected, observed)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55",
              linetype = 2, linewidth = 0.4) +
  geom_point_big(size = 0.5, alpha = 0.55, colour = PAL[1]) +
  geom_text(data = lt, aes(x = -Inf, y = Inf, label = lab), hjust = -0.12,
            vjust = 1.6, size = 2.9, inherit.aes = FALSE) +
  facet_wrap(~ panel, ncol = 3, scales = "free") +
  labs(title = "Quantile-quantile plots under the frozen model",
       subtitle = "Cohort-specific number of methylation PCs as pre-specified in decision D12",
       x = expression(Expected~-log[10]*(P)),
       y = expression(Observed~-log[10]*(P))) +
  theme_paper()
save_fig(p, "fig21_final_qq", width = 10, height = 6)

log_stop()
