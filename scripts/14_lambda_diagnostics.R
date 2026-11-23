# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 14: diagnosis of the inflation factor - find the source of the runaway lambda
#
# The results of 13 cannot be used as they are. Lambda ranges from 0.633 to 3.197 across the six data sets;
# in GSE145361 (PD discovery set) lambda = 3.197 with 56,773 probes at FDR < 0.05.
# This number conflicts directly with published results: in a larger meta-analysis that includes these samples,
# Vallerga et al. reported only a single-digit number of significant sites. The 56,773 are therefore almost certainly inflation,
# not discovery.
#
# This directly threatens the main criterion. The criterion compares effect sizes across cohorts;
# if one set carries unadjusted structure, the structure inflates the test statistics and also leaves its imprint on
# the effect sizes themselves; if two data sets share the same kind of technical structure, the comparison yields a
# completely spurious "shared signal".
#
# This step is diagnostic only and changes no conclusion:
#   Part A  balance between batch and case-control status. If a cohort put cases and controls on different
#           arrays, covariates cannot rescue it, and this must be known first.
#   Part B  model ladder. Covariates are added one at a time to see how lambda moves,
#           locating the term that causes the inflation and the term that removes it.
#
# To limit run time the ladder estimates lambda on a random subset of probes. Lambda is a median
# statistic and fifty thousand probes estimate it accurately.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/14_lambda_diagnostics.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "14_lambda_diagnostics"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "14_lambda")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260914)
N_PROBE_LAMBDA <- 50000L   # Number of random probes used to estimate lambda
N_MPC          <- 20L      # Candidate number of methylation principal components
MPC_GROUP_P    <- 0.01     # Principal components significantly associated with case-control status are not modelled
BETA_CLIP      <- 0.001
DROP_CT        <- "ct_Baso"
BARCODE        <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

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
bal <- list(); lad <- list()

for (f in files) {
  obj <- readRDS(f)
  g <- obj$gse
  say("==== ", g, " (", ROLE[g], ") ====")
  ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(colnames(obj$beta), cov$gsm), , drop = FALSE]

  grp <- factor(as.character(ph$group), levels = c("control", "case"))

  # ---------- Part A: balance between batch and group ----------
  bvar <- NULL; bname <- NA_character_
  if ("scan_date" %in% names(ph)) {
    bvar <- as.character(ph$scan_date); bname <- "scan date"
  } else {
    sl <- slide_of(g, colnames(obj$beta))
    if (sum(!is.na(sl)) > 0.95 * length(sl)) { bvar <- sl; bname <- "slide" }
  }
  if (!is.null(bvar)) {
    tb <- table(bvar, grp)
    pure <- sum(rowSums(tb > 0) == 1L)
    ct <- suppressWarnings(stats::chisq.test(tb))
    say("  Batch type: ", bname, ", levels ", nrow(tb),
        ", batches with a single group ", pure, " (", sprintf("%.0f%%", 100 * pure / nrow(tb)), ")")
    say("  Test of independence between batch and group: p = ", format.pval(ct$p.value, digits = 3))
    frac_case <- tb[, "case"] / rowSums(tb)
    say("  Proportion of cases per batch: median ", sprintf("%.2f", median(frac_case)),
        ", quartiles ", sprintf("%.2f", quantile(frac_case, 0.25)),
        " to ", sprintf("%.2f", quantile(frac_case, 0.75)))
    bal[[g]] <- data.frame(gse = g, batch = rownames(tb),
                           n = as.integer(rowSums(tb)),
                           frac_case = as.numeric(frac_case),
                           stringsAsFactors = FALSE)
  } else {
    say("  No batch information")
  }

  # ---------- covariate blocks ----------
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  sex <- if (is.na(sx)) NULL else factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  if (AGE_SOURCE[g] == "dnam") age <- cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    age <- suppressWarnings(as.numeric(ph[[ag]]))
  }
  cts <- cov[, setdiff(grep("^ct_", names(cov), value = TRUE), DROP_CT), drop = FALSE]
  smk <- if (isTRUE(USE_SMOKING[g])) cov$dnam_smoking else NULL
  # Note the use of base::all. EpiSmokEr exports a function named all;
  # once that package has been attached, unqualified all() reaches it.
  anc <- if (base::all(c("PC1", "PC2") %in% names(ph)))
    data.frame(anc1 = as.numeric(ph$PC1), anc2 = as.numeric(ph$PC2)) else NULL

  # ---------- methylation principal components ----------
  idx <- sample(seq_len(nrow(obj$beta)), min(N_PROBE_LAMBDA, nrow(obj$beta)))
  b <- obj$beta[idx, , drop = FALSE]
  b[!is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- log2(b / (1 - b))
  M <- M[stats::complete.cases(M), , drop = FALSE]
  say("  Probes used to estimate lambda: ", nrow(M), "")

  pr <- stats::prcomp(t(M), center = TRUE, scale. = FALSE)
  mpc <- pr$x[, seq_len(min(N_MPC, ncol(pr$x))), drop = FALSE]
  pg <- apply(mpc, 2, function(v)
    stats::anova(stats::lm(v ~ grp))$`Pr(>F)`[1])
  keep_pc <- which(pg >= MPC_GROUP_P)
  say("  Methylation principal components: ", ncol(mpc), ", of which associated with group (p < ", MPC_GROUP_P,
      "): ", sum(pg < MPC_GROUP_P), " not modelled, remaining ", length(keep_pc), "")
  mpc_keep <- if (length(keep_pc)) mpc[, keep_pc, drop = FALSE] else NULL

  # ---------- model ladder ----------
  blocks <- list(
    M0 = list(),
    M1 = list(cells = cts),
    M2 = list(cells = cts, smoking = smk),
    M3 = list(cells = cts, smoking = smk, batch = bvar),
    M4 = list(cells = cts, smoking = smk, batch = bvar, ancestry = anc),
    M5 = list(cells = cts, smoking = smk, batch = bvar, ancestry = anc,
              mpc = mpc_keep))
  desc <- c(M0 = "group+sex+age", M1 = "+cell fractions", M2 = "+smoking",
            M3 = "+batch", M4 = "+genetic PCs", M5 = "+methylation PCs")

  for (nm in names(blocks)) {
    d <- data.frame(group = grp)
    if (!is.null(sex)) d$sex <- sex
    d$age <- age
    bl <- blocks[[nm]]
    for (k in names(bl)) {
      v <- bl[[k]]
      if (is.null(v)) next
      if (is.data.frame(v) || is.matrix(v)) {
        v <- as.data.frame(v); names(v) <- paste0(k, "_", seq_len(ncol(v)))
        d <- cbind(d, v)
      } else {
        d[[k]] <- if (is.character(v)) factor(v) else v
      }
    }
    ok <- stats::complete.cases(d)
    dd <- d[ok, , drop = FALSE]
    des <- stats::model.matrix(~ ., data = dd)
    if (qr(des)$rank < ncol(des)) {
      say("  ", nm, " ", desc[nm], ": design matrix rank-deficient, skipped"); next
    }
    fit <- limma::eBayes(limma::lmFit(M[, ok, drop = FALSE], des))
    pv <- fit$p.value[, "groupcase"]
    lam <- lambda_of(pv)
    say("  ", nm, "  ", desc[nm], "  parameters ", ncol(des),
        "  Residual degrees of freedom: ", nrow(des) - ncol(des),
        "  lambda = ", sprintf("%.3f", lam))
    lad[[length(lad) + 1]] <- data.frame(
      gse = g, model = nm, label = unname(desc[nm]),
      n_par = ncol(des), lambda = lam, stringsAsFactors = FALSE)
  }
  say("")
  rm(obj, M, b, pr); invisible(gc())
}

ld <- do.call(rbind, lad)
write.csv(ld, file.path(OUT_DIR, "lambda_ladder.csv"), row.names = FALSE)
if (length(bal)) write.csv(do.call(rbind, bal),
                           file.path(OUT_DIR, "batch_balance.csv"), row.names = FALSE)

say("=== Ladder summary ===")
print(reshape(ld[, c("gse", "model", "lambda")], idvar = "gse",
              timevar = "model", direction = "wide"))

say("")
say("=== Plotting ===")
ld$panel <- factor(paste0(ld$gse, "\n", ROLE[ld$gse]),
                   levels = paste0(ORDER, "\n", ROLE[ORDER]))
ld$model <- factor(ld$model, levels = paste0("M", 0:5))
p <- ggplot(ld, aes(model, lambda, group = gse)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_line(colour = PAL[1], linewidth = 0.5) +
  geom_point(colour = PAL[1], size = 1.6) +
  geom_text(aes(label = sprintf("%.2f", lambda)), vjust = -0.9, size = 2.5) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  labs(title = "Genomic inflation across a ladder of models",
       subtitle = "M0 group+sex+age; M1 +cell composition; M2 +smoking; M3 +batch; M4 +genetic PCs; M5 +methylation PCs",
       x = NULL, y = expression(lambda)) +
  theme_paper()
save_fig(p, "fig15_lambda_ladder", width = 10, height = 6)

if (length(bal)) {
  bb <- do.call(rbind, bal)
  bb$panel <- factor(paste0(bb$gse, "\n", ROLE[bb$gse]),
                     levels = paste0(ORDER, "\n", ROLE[ORDER]))
  p <- ggplot(bb, aes(frac_case)) +
    geom_histogram(bins = 20, fill = PAL[1], colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey45", linewidth = 0.4) +
    facet_wrap(~ panel, ncol = 3, scales = "free_y") +
    labs(title = "Case fraction within each processing batch",
         subtitle = "Batches concentrated near 0 or 1 would mean batch is confounded with disease status",
         x = "Proportion of cases in the batch", y = "Number of batches") +
    theme_paper()
  save_fig(p, "fig16_batch_balance", width = 10, height = 5.5)
}

log_stop()
