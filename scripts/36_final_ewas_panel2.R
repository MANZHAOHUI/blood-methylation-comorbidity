# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 36: formal EWAS of the six additional cohorts
#
# Same method as step 17 of the first round:
#   methylation principal components are recomputed on the full filtered matrix, not on the fifty-thousand-probe subset of 35;
#   each set uses the K selected in 35; genomic control is applied according to D13 when the inflation factor exceeds 1.10,
#   adjusting P values only and not rescaling effect sizes - effect sizes are what is compared across cohorts, and rescaling would destroy comparability.
#
# The principal component algorithm is that of 17: components in the sample direction are needed, a direct singular value decomposition
# of a matrix of four hundred thousand by several hundred is wasteful, so the Gram matrix between samples (n x n) is computed first
# and then eigen-decomposed; the two are mathematically equivalent.
#
# One self-check: the inflation factor recomputed on the full matrix should be close to the one estimated by 35 on the subset.
# A large difference would mean that the subset of fifty thousand probes is not representative and the choice of K must be redone.
#
# The design matrices are saved for reuse by the later permutation tests, avoiding recomputation of the principal components.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/36_final_ewas_panel2.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "36_final_ewas_panel2"

RAW_DIR <- file.path(DATA_ROOT, "raw")
SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
SEL_DIR <- file.path(PROJ_ROOT, "output", "35_model")
OUT_DIR <- file.path(PROJ_ROOT, "output", "36_ewas")
for (d in c(RES_DIR, OUT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260917)
FORCE       <- FALSE
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
MPC_GROUP_P <- 0.01
GC_CUT      <- 1.10
LAM_TOL     <- 0.15     # Allowed difference between the inflation factors estimated on the full matrix and on the subset

COHORTS <- list(
  GSE84727  = list(disease = "schizophrenia", age = "dnam"),
  GSE80417  = list(disease = "schizophrenia", age = "chrono"),
  GSE198904 = list(disease = "depression",    age = "chrono"),
  GSE273970 = list(disease = "depression",    age = "dnam"),
  GSE338167 = list(disease = "alzheimer",     age = "chrono",
                   extra = c("ethnicity", "cohort")),
  GSE144858 = list(disease = "alzheimer",     age = "chrono")
)

SEL <- read.csv(file.path(SEL_DIR, "chosen_K_panel2.csv"), stringsAsFactors = FALSE)
say("Selected settings of 35 loaded: ", base::nrow(SEL), "")

lambda_of <- function(p) {
  cs <- stats::qchisq(p, df = 1, lower.tail = FALSE)
  stats::median(cs, na.rm = TRUE) / stats::qchisq(0.5, 1)
}
sample_pcs <- function(M, npc) {
  Mc <- M - base::rowMeans(M, na.rm = TRUE)
  G <- base::crossprod(Mc)
  rm(Mc); invisible(gc())
  ev <- base::eigen(G, symmetric = TRUE)
  k <- base::min(npc, base::sum(ev$values > 1e-8))
  sc <- ev$vectors[, base::seq_len(k), drop = FALSE] %*%
    base::diag(base::sqrt(base::pmax(ev$values[base::seq_len(k)], 0)), nrow = k)
  base::rownames(sc) <- base::colnames(M)
  base::colnames(sc) <- base::paste0("PC", base::seq_len(k))
  list(scores = sc, var = ev$values / base::sum(ev$values))
}
get_slide <- function(g, gsm_want) {
  hd <- mr_header(RAW_DIR, g)
  gsm <- mr_split(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  for (t in base::unique(base::sub("\t.*$", "", base::grep("^!Sample_", hd, value = TRUE)))) {
    for (one in base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)) {
      v <- mr_split(one)[-1]
      if (length(v) != length(gsm)) next
      hit <- base::regexpr(MR_BARCODE, v) > 0
      if (base::mean(hit) > 0.95) {
        bc <- base::rep(NA_character_, length(v))
        bc[hit] <- base::regmatches(v, base::regexpr(MR_BARCODE, v))
        return(base::sub("_R[0-9]{2}C[0-9]{2}$", "",
                         bc[base::match(gsm_want, gsm)]))
      }
    }
  }
  NULL
}
get_char <- function(g, tag, gsm_want) {
  hd <- mr_header(RAW_DIR, g)
  gsm <- mr_split(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  for (one in base::grep("^!Sample_characteristics_ch1", hd, value = TRUE)) {
    v <- mr_split(one)[-1]
    vv <- v[nzchar(v) & base::grepl(":", v, fixed = TRUE)]
    if (!length(vv)) next
    if (tolower(base::trimws(base::sub(":.*$", "", vv[1]))) == tolower(tag))
      return(base::trimws(base::sub("^[^:]*:", "", v))[base::match(gsm_want, gsm)])
  }
  NULL
}

summ <- list(); qqd <- list()
for (g in names(COHORTS)) {
  C <- COHORTS[[g]]
  out_rds <- file.path(RES_DIR, base::paste0(g, "_ewas_final.rds"))
  sel <- SEL[SEL$gse == g, , drop = FALSE]
  if (!base::nrow(sel)) { say("==== ", g, " ====  no selected setting in 35, skipped"); next }
  K <- sel$K[1]
  say("==== ", g, " (", C$disease, ")  K = ", K, " ====")
  if (file.exists(out_rds) && !FORCE) { say("  Exists, skipped (set FORCE to rebuild)"); next }

  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  ph <- obj$pheno
  gsmv <- base::colnames(obj$beta)
  cov <- read.csv(file.path(COV_DIR, base::paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[base::match(gsmv, cov$gsm), , drop = FALSE]

  b <- obj$beta
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
  keep_p <- stats::complete.cases(M)
  M <- M[keep_p, , drop = FALSE]
  say("  Probes ", base::nrow(M), " (removed for missing values: ", base::sum(!keep_p), "), samples ", base::ncol(M))

  d <- data.frame(group = factor(as.character(ph$group), levels = c("control", "case")))
  d$sex <- factor(base::toupper(base::substr(base::trimws(as.character(ph$sex)), 1, 1)))
  d$age <- if (C$age == "dnam") cov$dnam_age_horvath else
    base::suppressWarnings(base::as.numeric(ph$age_clean))
  d$smoking <- cov$dnam_smoking
  d <- cbind(d, cov[, base::setdiff(base::grep("^ct_", base::names(cov), value = TRUE),
                                    DROP_CT), drop = FALSE])
  sl <- get_slide(g, gsmv)
  if (!is.null(sl) && base::mean(base::is.na(sl)) <= 0.05 &&
      length(base::unique(sl)) <= 0.25 * length(gsmv)) {
    d$batch <- factor(sl)
    say("  Batch: arrays ", length(base::unique(sl)), "")
  } else say("  No batch term")
  if (!is.null(C$extra)) for (tg in C$extra) {
    v <- get_char(g, tg, gsmv)
    if (is.null(v)) next
    tb <- base::table(v)
    if (length(tb) >= 2L && length(tb) <= 12L) {
      d[[tg]] <- factor(v); say("  Additional covariate '", tg, "'", length(tb), " levels")
    }
  }

  if (K > 0) {
    say("  Computing sample principal components on the full matrix ...")
    t0 <- Sys.time()
    pcs <- sample_pcs(M, K + 25L)
    say("  Time: ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))),
        " min; the first ", base::ncol(pcs$scores), " components explain cumulatively ",
        sprintf("%.1f%%", 100 * base::sum(pcs$var[base::seq_len(base::ncol(pcs$scores))])))
    pg <- base::apply(pcs$scores, 2, function(v)
      stats::anova(stats::lm(v ~ d$group))$`Pr(>F)`[1])
    avail <- base::which(pg >= MPC_GROUP_P)
    take <- utils::head(avail, K)
    say("  Components associated with group: ", base::sum(pg < MPC_GROUP_P), " excluded; used: ", length(take), "")
    if (length(take) < K) say("  [NOTE] fewer usable components than K; actually used: ", length(take), "")
    pcd <- base::as.data.frame(pcs$scores[, take, drop = FALSE])
    base::names(pcd) <- base::paste0("mpc", base::seq_along(take))
    d <- cbind(d, pcd)
    rm(pcs); invisible(gc())
  }

  ok <- stats::complete.cases(d)
  if (base::sum(!ok)) say("  Removed because of missing covariates: ", base::sum(!ok), " samples")
  d <- d[ok, , drop = FALSE]; M <- M[, ok, drop = FALSE]
  des <- stats::model.matrix(~ ., data = d)
  q <- qr(des)
  if (q$rank < base::ncol(des)) {
    keepc <- base::sort(q$pivot[base::seq_len(q$rank)])
    say("  Collinear columns removed: ", base::ncol(des) - length(keepc), "",
        base::paste(utils::head(base::colnames(des)[base::setdiff(
          base::seq_len(base::ncol(des)), keepc)], 6), collapse = ", "))
    des <- des[, keepc, drop = FALSE]
  }
  tb <- base::table(d$group)
  say("  In the model: controls ", tb[["control"]], ", cases ", tb[["case"]],
      "; parameters ", base::ncol(des), ", residual degrees of freedom ", base::nrow(des) - base::ncol(des))

  fit <- limma::eBayes(limma::lmFit(M, des))
  tt <- limma::topTable(fit, coef = "groupcase", number = Inf, sort.by = "none")
  res <- data.frame(probe = base::rownames(tt), effect = tt$logFC,
                    se = fit$stdev.unscaled[, "groupcase"] * base::sqrt(fit$s2.post),
                    t = tt$t, p = tt$P.Value, stringsAsFactors = FALSE)
  lam <- lambda_of(res$p)
  say("  lambda on the full matrix = ", sprintf("%.3f", lam),
      "; estimated by 35 on the subset ", sprintf("%.3f", sel$lambda[1]),
      ", difference ", sprintf("%.3f", base::abs(lam - sel$lambda[1])),
      if (base::abs(lam - sel$lambda[1]) > LAM_TOL)
        "  [NOTE] large difference; the representativeness of the subset is doubtful" else "")
  if (lam > GC_CUT) {
    cs <- stats::qchisq(res$p, 1, lower.tail = FALSE) / lam
    res$p_adj_scale <- stats::pchisq(cs, 1, lower.tail = FALSE); gcf <- TRUE
  } else { res$p_adj_scale <- res$p; gcf <- FALSE }
  res$fdr <- stats::p.adjust(res$p_adj_scale, method = "BH")
  say("  Genomic control: ", if (gcf) "applied" else "not applied",
      "; probes with FDR < 0.05: ", base::sum(res$fdr < 0.05, na.rm = TRUE))

  base::saveRDS(res, out_rds)
  base::saveRDS(list(design = des, gsm = base::colnames(M), probes = res$probe,
                     group = d$group, gse = g),
                file.path(RES_DIR, base::paste0(g, "_design.rds")))
  write.csv(utils::head(res[base::order(res$p), ], 200),
            file.path(OUT_DIR, base::paste0(g, "_top200.csv")), row.names = FALSE)

  summ[[length(summ) + 1]] <- data.frame(
    gse = g, disease = C$disease, K = K, n_control = tb[["control"]],
    n_case = tb[["case"]], n_probe = base::nrow(res), lambda_full = lam,
    lambda_subset = sel$lambda[1], gc = gcf,
    n_fdr05 = base::sum(res$fdr < 0.05, na.rm = TRUE), stringsAsFactors = FALSE)
  o <- base::order(res$p)
  qqd[[length(qqd) + 1]] <- data.frame(
    gse = g, disease = C$disease, lambda = lam,
    expected = -base::log10(stats::ppoints(base::nrow(res))),
    observed = -base::log10(res$p[o]), stringsAsFactors = FALSE)
  rm(M, fit, des); invisible(gc())
  say("")
}

if (length(summ)) {
  tab <- base::do.call(rbind, summ)
  write.csv(tab, file.path(OUT_DIR, "ewas_summary_panel2.csv"), row.names = FALSE)
  say("=== Summary ==="); print(tab)
}
if (length(qqd)) {
  qd <- base::do.call(rbind, qqd)
  qd$panel <- factor(base::paste0(qd$gse, "\n", qd$disease),
                     levels = base::unique(base::paste0(qd$gse, "\n", qd$disease)))
  lt <- base::unique(qd[, c("panel", "lambda")])
  lt$lab <- base::paste0("lambda = ", sprintf("%.3f", lt$lambda))
  say("")
  say("=== Plotting ===")
  p <- ggplot(qd, aes(expected, observed)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2, linewidth = 0.4) +
    geom_point_big(size = 0.5, alpha = 0.55, colour = PAL[1]) +
    geom_text(data = lt, aes(x = -Inf, y = Inf, label = lab), hjust = -0.12,
              vjust = 1.6, size = 2.9, inherit.aes = FALSE) +
    facet_wrap(~ panel, ncol = 3, scales = "free") +
    labs(title = "Quantile-quantile plots for the second-round cohorts",
         subtitle = "Cohort-specific number of methylation PCs as pre-specified in decision D12",
         x = expression(Expected~-log[10]*(P)),
         y = expression(Observed~-log[10]*(P))) +
    theme_paper()
  save_fig(p, "fig40_qq_panel2", width = 10, height = 6)
}
log_stop()
