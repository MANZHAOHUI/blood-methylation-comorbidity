# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 15: settle the final model
#
# Step 14 left three problems that must be solved:
#   (1) GSE145361 still has lambda = 1.70 after adding methylation principal components;
#       none of the measurable covariates (cells, smoking, batch, genetic PCs) brings it down.
#   (2) GSE152027 hardly moves along the whole ladder (1.64 to 1.35); the negative control itself
#       carries unadjusted structure and would contaminate the criterion.
#   (3) GSE87648 drops to lambda = 0.48 after adding cell fractions, which is over-adjustment.
#       In IBD the change in leukocyte composition is itself part of the disease phenotype;
#       it behaves more like a mediator than a confounder, and removing it removes signal as well.
#
# This step does three things:
#   A  scan the number K of methylation principal components (0/5/10/20/30/50) to see whether lambda
#      returns to about 1 as K increases, and how many are needed.
#   B  run both "with cell fractions" and "without cell fractions". This is not just tuning:
#      the two answer different questions - without cell fractions is the overall peripheral signature,
#      with cell fractions is the signature within cell types. Both are reported in the end.
#   C  compute the variance inflation factor of the case-control term to check whether the low lambda of GSE87648
#      comes from collinearity between group and cell fractions.
#
# This step still yields no conclusion about sharing.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/15_model_selection.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "15_model_selection"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "15_model")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260914)
N_PROBE     <- 50000L
K_GRID      <- c(0L, 5L, 10L, 20L, 30L, 50L)
N_MPC_MAX   <- 60L
MPC_GROUP_P <- 0.01
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
BARCODE     <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

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
# Variance inflation factor of the case-control column: regress it on the other columns, VIF = 1/(1-R^2)
vif_group <- function(des) {
  j <- which(colnames(des) == "groupcase")
  if (!length(j) || ncol(des) < 3) return(NA_real_)
  x <- des[, setdiff(seq_len(ncol(des)), c(1L, j)), drop = FALSE]
  if (!ncol(x)) return(1)
  r2 <- summary(stats::lm(des[, j] ~ x))$r.squared
  1 / (1 - r2)
}

files <- list.files(SET_DIR, pattern = "_analysis_filtered\\.rds$", full.names = TRUE)
res <- list()

for (f in files) {
  obj <- readRDS(f); g <- obj$gse
  say("==== ", g, " (", ROLE[g], ") ====")
  ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(colnames(obj$beta), cov$gsm), , drop = FALSE]

  grp <- factor(as.character(ph$group), levels = c("control", "case"))
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  sex <- if (is.na(sx)) NULL else factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  if (AGE_SOURCE[g] == "dnam") age <- cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    age <- suppressWarnings(as.numeric(ph[[ag]]))
  }
  cts <- cov[, setdiff(grep("^ct_", names(cov), value = TRUE), DROP_CT), drop = FALSE]
  smk <- if (isTRUE(USE_SMOKING[g])) cov$dnam_smoking else NULL
  bvar <- if ("scan_date" %in% names(ph)) as.character(ph$scan_date) else {
    sl <- slide_of(g, colnames(obj$beta))
    if (sum(!is.na(sl)) > 0.95 * length(sl)) sl else NULL
  }

  idx <- sample(seq_len(nrow(obj$beta)), min(N_PROBE, nrow(obj$beta)))
  b <- obj$beta[idx, , drop = FALSE]
  b[!is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- log2(b / (1 - b)); M <- M[stats::complete.cases(M), , drop = FALSE]
  rm(b); invisible(gc())
  say("  Probes ", nrow(M), ", samples ", ncol(M))

  pr <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  nmp <- min(N_MPC_MAX, ncol(pr$x))
  mpc <- pr$x[, seq_len(nmp), drop = FALSE]
  pg <- apply(mpc, 2, function(v) stats::anova(stats::lm(v ~ grp))$`Pr(>F)`[1])
  ok_pc <- which(pg >= MPC_GROUP_P)
  say("  Candidate methylation principal components: ", nmp, ", associated with group: ", sum(pg < MPC_GROUP_P),
      " removed, usable ", length(ok_pc), "")
  rm(pr); invisible(gc())

  for (usecell in c(FALSE, TRUE)) {
    for (K in K_GRID) {
      d <- data.frame(group = grp)
      if (!is.null(sex)) d$sex <- sex
      d$age <- age
      if (!is.null(smk)) d$smoking <- smk
      if (!is.null(bvar)) d$batch <- factor(bvar)
      if (usecell) d <- cbind(d, cts)
      if (K > 0) {
        take <- utils::head(ok_pc, K)
        if (length(take)) {
          pcd <- as.data.frame(mpc[, take, drop = FALSE])
          names(pcd) <- paste0("mpc", seq_along(take))
          d <- cbind(d, pcd)
        }
      }
      keep <- stats::complete.cases(d)
      dd <- d[keep, , drop = FALSE]
      des <- stats::model.matrix(~ ., data = dd)
      rdf <- nrow(des) - ncol(des)
      if (qr(des)$rank < ncol(des) || rdf < 30L) {
        say("  ", if (usecell) "with cells" else "without cells", " K=", K,
            ": rank-deficient or insufficient residual degrees of freedom (", rdf, "), skipped"); next
      }
      fit <- limma::eBayes(limma::lmFit(M[, keep, drop = FALSE], des))
      lam <- lambda_of(fit$p.value[, "groupcase"])
      vf <- vif_group(des)
      say("  ", if (usecell) "with cells" else "without cells", "  K=", K,
          "  parameters ", ncol(des), "  Residual degrees of freedom: ", rdf,
          "  lambda = ", sprintf("%.3f", lam),
          "  VIF of group = ", sprintf("%.2f", vf))
      res[[length(res) + 1]] <- data.frame(
        gse = g, cells = usecell, K = K, n_par = ncol(des), rdf = rdf,
        lambda = lam, vif = vf, stringsAsFactors = FALSE)
    }
  }
  say("")
  rm(obj, M, mpc); invisible(gc())
}

rr <- do.call(rbind, res)
write.csv(rr, file.path(OUT_DIR, "model_selection.csv"), row.names = FALSE)

say("=== Setting closest to 1 for each data set ===")
for (g in unique(rr$gse)) {
  sub <- rr[rr$gse == g, ]
  for (uc in c(FALSE, TRUE)) {
    s2 <- sub[sub$cells == uc, ]
    if (!nrow(s2)) next
    best <- s2[which.min(abs(s2$lambda - 1)), ]
    say("  ", g, "  ", if (uc) "with cells" else "without cells",
        ": K = ", best$K, "  lambda = ", sprintf("%.3f", best$lambda),
        "  VIF = ", sprintf("%.2f", best$vif))
  }
}

say("")
say("=== Plotting ===")
rr$panel <- factor(paste0(rr$gse, "\n", ROLE[rr$gse]),
                   levels = paste0(ORDER, "\n", ROLE[ORDER]))
rr$model <- factor(ifelse(rr$cells, "with cell composition", "without cell composition"),
                   levels = c("without cell composition", "with cell composition"))
p <- ggplot(rr, aes(K, lambda, colour = model, group = model)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.5) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("without cell composition" = PAL[1],
                                 "with cell composition" = PAL[2]), name = NULL) +
  scale_x_continuous(breaks = K_GRID) +
  labs(title = "Genomic inflation against the number of methylation PCs",
       subtitle = "All models include group, sex, age, smoking and batch where available",
       x = "Number of methylation principal components included",
       y = expression(lambda)) +
  theme_paper()
save_fig(p, "fig17_lambda_vs_npc", width = 10, height = 6)

p <- ggplot(rr, aes(K, vif, colour = model, group = model)) +
  geom_hline(yintercept = 5, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.55) + geom_point(size = 1.5) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("without cell composition" = PAL[1],
                                 "with cell composition" = PAL[2]), name = NULL) +
  scale_x_continuous(breaks = K_GRID) +
  labs(title = "Variance inflation of the case-control term",
       subtitle = "High values mean the covariates are collinear with disease status, which costs power",
       x = "Number of methylation principal components included",
       y = "VIF of the case-control term") +
  theme_paper()
save_fig(p, "fig18_vif_group", width = 10, height = 6)

log_stop()
