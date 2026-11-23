# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 54: rebuild the model of GSE87648 after removing duplicate samples
#
# Basis: with the rs probes of the array, step 48 found two pairs of samples within GSE87648 with concordant genotypes
# (concordance 0.968 and 0.921), i.e. the same person included twice. Duplicate samples
# overstate the effective sample size of the cohort and understate the variance; only one sample per pair is kept.
#
# Retention rule (written down beforehand, without looking at results): of each pair, keep the sample with the lower quality index of step 43
# (proportion of beta values between 0.3 and 0.7 among fifty thousand random autosomal probes), i.e. the one with better data
# quality; if equal, keep the column that comes first in the matrix.
#
# This step does three things:
#   (1) write the exclusion list output/54_dedup/excluded_samples.csv;
#   (2) move the existing GSE87648_ewas_final.rds and GSE87648_design.rds to
#       ewas_final/archive_before_dedup/ for reference;
#   (3) re-determine K on the de-duplicated samples by the rule of 37 and rebuild the model of
#       GSE87648 as in 38, writing it back under the same file names. The covariates are built exactly as in 38.
#
# Downstream steps 39-53 all take their samples from gsm of *_design.rds and use the de-duplicated
# sample set automatically on a re-run; 43 reads the analysis set directly and now reads the exclusion list of this step.
#
# Environment: Windows + RStudio. About 2 minutes. Safe to run repeatedly: the archive is not overwritten and the design is read from the archive.
# Location: <project root>/scripts/54_dedup_gse87648.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "54_dedup_gse87648"

RAW_DIR <- file.path(DATA_ROOT, "raw")
SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
ARC_DIR <- file.path(RES_DIR, "archive_before_dedup")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
QC_DIR  <- file.path(PROJ_ROOT, "output", "43_qc")
ID_DIR  <- file.path(PROJ_ROOT, "output", "48_identity")
TAB_DIR <- file.path(PROJ_ROOT, "output", "38_rebuild")
OUT_DIR <- file.path(PROJ_ROOT, "output", "54_dedup")
for (d in c(ARC_DIR, OUT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
}

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260917)
GSE         <- "GSE87648"
K_GRID      <- c(0L, 5L, 10L, 20L, 30L, 50L)
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
MPC_GROUP_P <- 0.01
LAM_LO      <- 0.95
LAM_HI      <- 1.10
GC_CUT      <- 1.10
AGE_SRC     <- "chrono"     # Same settings as 37/38 for GSE87648
USE_SMOKING <- TRUE

lambda_of <- function(p) {
  cs <- stats::qchisq(p, df = 1, lower.tail = FALSE)
  stats::median(cs, na.rm = TRUE) / stats::qchisq(0.5, 1)
}
sample_pcs <- function(M, npc) {
  Mc <- M - base::rowMeans(M, na.rm = TRUE)
  G <- base::crossprod(Mc); rm(Mc); invisible(gc())
  ev <- base::eigen(G, symmetric = TRUE)
  k <- base::min(npc, base::sum(ev$values > 1e-8))
  sc <- ev$vectors[, base::seq_len(k), drop = FALSE] %*%
    base::diag(base::sqrt(base::pmax(ev$values[base::seq_len(k)], 0)), nrow = k)
  base::rownames(sc) <- base::colnames(M)
  base::colnames(sc) <- base::paste0("PC", base::seq_len(k))
  sc
}
fit_model <- function(M, d) {
  keep <- stats::complete.cases(d)
  des <- stats::model.matrix(~ ., data = d[keep, , drop = FALSE])
  q <- qr(des)
  if (q$rank < base::ncol(des))
    des <- des[, base::sort(q$pivot[base::seq_len(q$rank)]), drop = FALSE]
  fit <- limma::eBayes(limma::lmFit(M[, keep, drop = FALSE], des))
  list(fit = fit, des = des, keep = keep,
       lam = lambda_of(fit$p.value[, "groupcase"]),
       npar = base::ncol(des), rdf = base::nrow(des) - base::ncol(des))
}
slide_from_map <- function(g, gsm_want) {
  f <- file.path(MAP_DIR, base::paste0(g, "_sample_map.csv"))
  if (!file.exists(f)) return(NULL)
  m <- read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
  bc <- base::rep(NA_character_, base::nrow(m))
  hit <- base::regexpr(MR_BARCODE, m$col_name) > 0
  if (!base::any(hit)) return(NULL)
  bc[hit] <- base::regmatches(m$col_name, base::regexpr(MR_BARCODE, m$col_name))
  base::sub("_R[0-9]{2}C[0-9]{2}$", "", bc)[base::match(gsm_want, m$gsm)]
}

# --------------------------------------------------------------
say("=== 1. Determine the samples to exclude ===")
ID <- read.csv(file.path(ID_DIR, "identity_matches.csv"), stringsAsFactors = FALSE)
ID <- ID[ID$a == GSE & ID$b == GSE, , drop = FALSE]
if (!base::nrow(ID)) stop("The match table of step 48 contains no within-cohort duplicate pairs for ", GSE, "")
QS <- read.csv(file.path(QC_DIR, "bimodality_per_sample.csv"), stringsAsFactors = FALSE)
QS <- QS[QS$gse == GSE, , drop = FALSE]
# Design matrix before de-duplication: in ewas_final/ on the first run, already in the archive directory on later runs
arc_des <- file.path(ARC_DIR, base::paste0(GSE, "_design.rds"))
D_old <- base::readRDS(if (file.exists(arc_des)) arc_des else
                         file.path(RES_DIR, base::paste0(GSE, "_design.rds")))
obj <- base::readRDS(file.path(SET_DIR, base::paste0(GSE, "_analysis_filtered.rds")))
col_order <- base::colnames(obj$beta)

EX_FILE <- file.path(OUT_DIR, "excluded_samples.csv")
excl <- list()
if (file.exists(EX_FILE)) {
  EX <- read.csv(EX_FILE, stringsAsFactors = FALSE)
  say("  Exclusion list exists (", base::nrow(EX), " samples); reused, not decided again")
} else for (i in base::seq_len(base::nrow(ID))) {
  pair <- c(ID$sample_a[i], ID$sample_b[i])
  if (!base::all(pair %in% D_old$gsm))
    stop("Duplicate pair ", base::paste(pair, collapse = "/"), " is not entirely in the current design matrix; check step 48 first")
  fm <- QS$frac_mid[base::match(pair, QS$gsm)]
  if (base::any(base::is.na(fm))) stop("Step 43 has no quality index for ", base::paste(pair, collapse = "/"), "")
  drop_i <- if (fm[1] != fm[2]) base::which.max(fm) else
    base::which.max(base::match(pair, col_order))
  keep_i <- 3L - drop_i
  say(sprintf("  Pair %d: %s (quality index %.4f) and %s (%.4f), genotype concordance %.3f; kept %s, excluded %s",
              i, pair[1], fm[1], pair[2], fm[2], ID$concordance[i], pair[keep_i], pair[drop_i]))
  excl[[i]] <- data.frame(gse = GSE, gsm = pair[drop_i], kept_partner = pair[keep_i],
                          frac_mid_dropped = fm[drop_i], frac_mid_kept = fm[keep_i],
                          genotype_concordance = ID$concordance[i],
                          group_dropped = as.character(obj$pheno$group[base::match(pair[drop_i], col_order)]),
                          rule = "keep the member with the lower bimodality index",
                          stringsAsFactors = FALSE)
}
if (!file.exists(EX_FILE)) {
  EX <- base::do.call(rbind, excl)
  write.csv(EX, EX_FILE, row.names = FALSE)
}
say("  removed ", base::nrow(EX), " samples (group: ", base::paste(EX$group_dropped, collapse = ", "),
    "); list written")

# --------------------------------------------------------------
say("")
say("=== 2. Archive the results and design matrix from before de-duplication ===")
for (suf in c("_ewas_final.rds", "_design.rds")) {
  src <- file.path(RES_DIR, base::paste0(GSE, suf))
  dst <- file.path(ARC_DIR, base::paste0(GSE, suf))
  if (file.exists(dst)) {
    say("  ", base::basename(dst), " already in the archive directory, not overwritten")
  } else {
    if (!file.copy(src, dst)) stop("Archiving failed: ", src)
    say("  Archived ", base::basename(src))
  }
}
E_old <- base::readRDS(file.path(ARC_DIR, base::paste0(GSE, "_ewas_final.rds")))

# --------------------------------------------------------------
say("")
say("=== 3. Re-determine K after de-duplication and rebuild the model ===")
ph <- obj$pheno
gsmv <- base::setdiff(col_order, EX$gsm)
ph <- ph[base::match(gsmv, col_order), , drop = FALSE]
cov <- read.csv(file.path(COV_DIR, base::paste0(GSE, "_covariates.csv")), stringsAsFactors = FALSE)
cov <- cov[base::match(gsmv, cov$gsm), , drop = FALSE]
b <- obj$beta[, gsmv, drop = FALSE]; rm(obj); invisible(gc())
b[!base::is.finite(b)] <- NA
b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
M <- M[stats::complete.cases(M), , drop = FALSE]
say("  Probes ", base::nrow(M), ", samples ", base::ncol(M), " (analysis set before de-duplication ", length(col_order), ")")

d <- data.frame(group = factor(as.character(ph$group), levels = c("control", "case")))
sx <- base::grep("^(sex|gender)$", base::names(ph), ignore.case = TRUE, value = TRUE)[1]
d$sex <- factor(base::toupper(base::substr(base::trimws(as.character(ph[[sx]])), 1, 1)))
if (AGE_SRC == "dnam") d$age <- cov$dnam_age_horvath else {
  ac <- if ("age_clean" %in% base::names(ph)) ph$age_clean else
    ph[[base::grep("^age$|ageat", base::names(ph), ignore.case = TRUE, value = TRUE)[1]]]
  d$age <- base::suppressWarnings(base::as.numeric(ac))
}
if (isTRUE(USE_SMOKING) && "dnam_smoking" %in% base::names(cov)) d$smoking <- cov$dnam_smoking
d <- cbind(d, cov[, base::setdiff(base::grep("^ct_", base::names(cov), value = TRUE), DROP_CT), drop = FALSE])
if ("scan_date" %in% base::names(ph)) {
  d$batch <- factor(as.character(ph$scan_date))
} else {
  sl <- slide_from_map(GSE, gsmv)
  if (!is.null(sl) && base::mean(base::is.na(sl)) <= 0.05 &&
      length(base::unique(sl)) <= 0.25 * length(gsmv)) d$batch <- factor(sl)
}
say("  Covariates: ", base::paste(base::setdiff(base::names(d), "group"), collapse = ", "))

# Check against the design before de-duplication: apart from the two samples removed, the model structure must be identical
n_old <- base::nrow(D_old$design)
n_new <- base::sum(stats::complete.cases(d))
hb_old <- base::any(base::grepl("^batch", base::colnames(D_old$design)))
say("  Design before de-duplication: samples ", n_old, ", parameters ", base::ncol(D_old$design), ", batch term ", if (hb_old) "yes" else "none")
if (n_new != n_old - base::nrow(EX))
  stop("The number of samples in the model should be ", n_old - base::nrow(EX), ", obtained ", n_new, "; the covariate sources have changed, stopping")
if (hb_old != ("batch" %in% base::names(d))) stop("The batch term differs from that before de-duplication, stopping")
say("  Model comparison passed: in the model ", n_new, " = before de-duplication ", n_old, " minus ", base::nrow(EX), " samples")

lam <- rep(NA_real_, length(K_GRID))
r0 <- fit_model(M, d); lam[1] <- r0$lam
say(sprintf("  K= 0  parameters %3d  residual df %4d  lambda = %.3f", r0$npar, r0$rdf, r0$lam))
chosen <- NULL
if (!base::is.na(r0$lam) && r0$lam >= LAM_LO && r0$lam <= LAM_HI) {
  say("  K=0 is within the interval; settled by D12"); chosen <- list(K = 0L, r = r0, take = integer(0))
} else {
  say("  K=0 is outside the interval; computing principal components on the full matrix ...")
  t0 <- Sys.time()
  sc <- sample_pcs(M, base::max(K_GRID) + 25L)
  say("  Time: ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  pg <- base::apply(sc, 2, function(v) stats::anova(stats::lm(v ~ d$group))$`Pr(>F)`[1])
  avail <- base::which(pg >= MPC_GROUP_P)
  say("  Components associated with group: ", base::sum(pg < MPC_GROUP_P), " removed, usable ", length(avail), "")
  fits <- list()
  for (i in 2:length(K_GRID)) {
    take <- utils::head(avail, K_GRID[i])
    pcd <- base::as.data.frame(sc[, take, drop = FALSE])
    base::names(pcd) <- base::paste0("mpc", base::seq_along(take))
    r <- fit_model(M, cbind(d, pcd)); lam[i] <- r$lam; fits[[i]] <- list(r = r, take = take, d = cbind(d, pcd))
    say(sprintf("  K=%2d  parameters %3d  residual df %4d  lambda = %.3f", K_GRID[i], r$npar, r$rdf, r$lam))
    if (!base::is.na(r$lam) && r$lam >= LAM_LO && r$lam <= LAM_HI) {
      say("  within the interval; by D12 the smallest K is taken and no larger K is tried"); break
    }
  }
  inb <- base::which(!base::is.na(lam) & lam >= LAM_LO & lam <= LAM_HI)
  j <- if (length(inb)) inb[1] else base::which.min(base::abs(lam - 1))
  if (j == 1L) chosen <- list(K = 0L, r = r0, take = integer(0), d = d) else
    chosen <- list(K = K_GRID[j], r = fits[[j]]$r, take = fits[[j]]$take, d = fits[[j]]$d)
  rm(sc); invisible(gc())
}
if (is.null(chosen$d)) chosen$d <- d
newK <- chosen$K; r <- chosen$r
gcf <- !(r$lam >= LAM_LO && r$lam <= LAM_HI)
if (gcf && r$lam <= GC_CUT) gcf <- FALSE
say("  Determined: K = ", newK, ", lambda = ", sprintf("%.3f", r$lam),
    ", genomic control ", if (gcf) "applied" else "not applied",
    " (before de-duplication K = 5, lambda = 1.086, not applied)")

des <- r$des; fit <- r$fit; dd <- chosen$d[r$keep, , drop = FALSE]
tb <- base::table(dd$group)
say("  In the model: controls ", tb[["control"]], ", cases ", tb[["case"]],
    "; parameters ", base::ncol(des), ", residual degrees of freedom ", base::nrow(des) - base::ncol(des))
tt <- limma::topTable(fit, coef = "groupcase", number = Inf, sort.by = "none")
res <- data.frame(probe = base::rownames(tt), effect = tt$logFC,
                  se = fit$stdev.unscaled[, "groupcase"] * base::sqrt(fit$s2.post),
                  t = tt$t, p = tt$P.Value, stringsAsFactors = FALSE)
if (gcf) {
  cs <- stats::qchisq(res$p, 1, lower.tail = FALSE) / r$lam
  res$p_adj_scale <- stats::pchisq(cs, 1, lower.tail = FALSE)
} else res$p_adj_scale <- res$p
res$fdr <- stats::p.adjust(res$p_adj_scale, method = "BH")
say("  Probes with FDR < 0.05: ", base::sum(res$fdr < 0.05, na.rm = TRUE),
    " (before de-duplication ", base::sum(E_old$fdr < 0.05, na.rm = TRUE), ")")

gsm_keep <- base::colnames(M)[r$keep]
base::saveRDS(res, file.path(RES_DIR, base::paste0(GSE, "_ewas_final.rds")))
base::saveRDS(list(design = des, gsm = gsm_keep, probes = res$probe,
                   group = dd$group, gse = GSE),
              file.path(RES_DIR, base::paste0(GSE, "_design.rds")))
write.csv(utils::head(res[base::order(res$p), ], 200),
          file.path(OUT_DIR, base::paste0(GSE, "_top200.csv")), row.names = FALSE)
say("  Written back: ", GSE, " results and design matrix (", length(gsm_keep), " samples)")

# --------------------------------------------------------------
say("")
say("=== 4. Comparison with the results before de-duplication ===")
cm <- base::intersect(res$probe, E_old$probe)
r_eff <- stats::cor(res$effect[base::match(cm, res$probe)], E_old$effect[base::match(cm, E_old$probe)])
say(sprintf("  Pearson correlation of the effect vector with that before de-duplication %.4f (%d probes)", r_eff, length(cm)))
EB <- base::readRDS(file.path(RES_DIR, "GSE112611_ewas_final.rds"))
cm2 <- base::intersect(res$probe, EB$probe)
rho_new <- stats::cor(res$effect[base::match(cm2, res$probe)], EB$effect[base::match(cm2, EB$probe)], method = "spearman")
cm3 <- base::intersect(E_old$probe, EB$probe)
rho_old <- stats::cor(E_old$effect[base::match(cm3, E_old$probe)], EB$effect[base::match(cm3, EB$probe)], method = "spearman")
say(sprintf("  Within-disease concordance with GSE112611 (observed values only): before de-duplication %+.4f, after %+.4f; the formal P value is recomputed by step 45", rho_old, rho_new))

ladder <- data.frame(gse = GSE, K = K_GRID, lambda = lam, stringsAsFactors = FALSE)
write.csv(ladder, file.path(OUT_DIR, "lambda_ladder_gse87648.csv"), row.names = FALSE)
FM <- read.csv(file.path(TAB_DIR, "final_model_table.csv"), stringsAsFactors = FALSE)
i <- base::which(FM$gse == GSE)
FM$K_new[i] <- newK; FM$lambda_new[i] <- r$lam; FM$genomic_control[i] <- gcf
write.csv(FM, file.path(OUT_DIR, "final_model_table.csv"), row.names = FALSE)
cmp <- data.frame(item = c("samples_in_model", "controls", "cases", "K", "lambda", "genomic_control", "fdr05_probes"),
                  before = c(n_old, base::sum(D_old$group == "control"), base::sum(D_old$group == "case"),
                             5, 1.08624329245151, FALSE, base::sum(E_old$fdr < 0.05, na.rm = TRUE)),
                  after = c(length(gsm_keep), tb[["control"]], tb[["case"]], newK, r$lam, gcf,
                            base::sum(res$fdr < 0.05, na.rm = TRUE)), stringsAsFactors = FALSE)
cmp <- rbind(cmp, data.frame(item = "within_IBD_rho_observed", before = rho_old, after = rho_new),
             data.frame(item = "effect_vector_cor_before_vs_after", before = NA, after = r_eff))
write.csv(cmp, file.path(OUT_DIR, "before_after_summary.csv"), row.names = FALSE)
print(cmp)
say("")
say("The final model table of the ten sets was written to output/54_dedup/final_model_table.csv; step 53 reads the frozen inflation factors from it.")
say("Next, re-run steps 39 to 53 in the order given by step 55.")
log_stop()
