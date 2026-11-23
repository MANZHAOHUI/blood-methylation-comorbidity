# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 38: rebuild the EWAS of cohorts whose K changed
#
# 37 determined K for the ten cohorts on the full matrix and all model comparisons passed (sample number and
# batch term agree with the saved designs). Only GSE145361 needs a change:
# K from 30 to 50, inflation factor from 1.378 to 1.355. Neither is within
# 0.95 to 1.10; by D12 the one closer to 1 is taken and genomic control applied according to D13.
#
# This step re-decides nothing; it only rebuilds the flagged cohorts with the K written by 37
# and produces result files and design matrices in the same format as the other nine sets.
#
# The covariates are built exactly as in 37 - deliberately: the model of 37 has already passed the comparison with
# the original analysis, and copying it guarantees that the rebuilt results are comparable with the other nine sets.
#
# Environment: Windows + RStudio. The principal components of GSE145361 take about seven minutes
# Location: <project root>/scripts/38_rebuild_changed_ewas.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "38_rebuild_changed_ewas"

RAW_DIR <- file.path(DATA_ROOT, "raw")
SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
SEL_DIR <- file.path(PROJ_ROOT, "output", "37_reselect")
OUT_DIR <- file.path(PROJ_ROOT, "output", "38_rebuild")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260917)
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
MPC_GROUP_P <- 0.01
GC_CUT      <- 1.10

AGE_SRC <- list(GSE145361 = "dnam", GSE111629 = "chrono", GSE87648 = "chrono",
                GSE112611 = "chrono", GSE84727 = "dnam", GSE80417 = "chrono",
                GSE198904 = "chrono", GSE273970 = "dnam", GSE338167 = "chrono",
                GSE144858 = "chrono")
USE_SMOKING <- list(GSE145361 = TRUE, GSE111629 = TRUE, GSE87648 = TRUE,
                    GSE112611 = FALSE, GSE84727 = TRUE, GSE80417 = TRUE,
                    GSE198904 = TRUE, GSE273970 = TRUE, GSE338167 = TRUE,
                    GSE144858 = TRUE)
EXTRA <- list(GSE338167 = c("ethnicity", "cohort"))
DISEASE <- list(GSE145361 = "parkinson", GSE111629 = "parkinson",
                GSE87648 = "ibd", GSE112611 = "ibd",
                GSE84727 = "schizophrenia", GSE80417 = "schizophrenia",
                GSE198904 = "depression", GSE273970 = "depression",
                GSE338167 = "alzheimer", GSE144858 = "alzheimer")

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
slide_from_hdr <- function(g, gsm_want) {
  hd <- tryCatch(mr_header(RAW_DIR, g), error = function(e) NULL)
  if (is.null(hd)) return(NULL)
  gsm <- mr_split(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  for (t in base::unique(base::sub("\t.*$", "", base::grep("^!Sample_", hd, value = TRUE)))) {
    for (one in base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)) {
      v <- mr_split(one)[-1]
      if (length(v) != length(gsm)) next
      hit <- base::regexpr(MR_BARCODE, v) > 0
      if (base::mean(hit) > 0.95) {
        bc <- base::rep(NA_character_, length(v))
        bc[hit] <- base::regmatches(v, base::regexpr(MR_BARCODE, v))
        return(base::sub("_R[0-9]{2}C[0-9]{2}$", "", bc[base::match(gsm_want, gsm)]))
      }
    }
  }
  NULL
}
get_char <- function(g, tag, gsm_want) {
  hd <- tryCatch(mr_header(RAW_DIR, g), error = function(e) NULL)
  if (is.null(hd)) return(NULL)
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

SEL <- read.csv(file.path(SEL_DIR, "K_reselect.csv"), stringsAsFactors = FALSE)
todo <- SEL[SEL$changed %in% c(TRUE, "TRUE"), , drop = FALSE]
say("Of the ten sets determined by 37, to rebuild: ", base::nrow(todo), ": ",
    base::paste(todo$gse, collapse = ", "))
if (!base::nrow(todo)) { say("No cohort needs to be rebuilt"); log_stop(); }

for (i in base::seq_len(base::nrow(todo))) {
  g <- todo$gse[i]; K <- todo$K_new[i]
  say("")
  say("==== ", g, "  K from ", todo$K_current[i], " to ", K, " ====")
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  ph <- obj$pheno; gsmv <- base::colnames(obj$beta)
  cov <- read.csv(file.path(COV_DIR, base::paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[base::match(gsmv, cov$gsm), , drop = FALSE]
  b <- obj$beta; rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
  M <- M[stats::complete.cases(M), , drop = FALSE]
  say("  Probes ", base::nrow(M), ", samples ", base::ncol(M))

  d <- data.frame(group = factor(as.character(ph$group), levels = c("control", "case")))
  sx <- base::grep("^(sex|gender)$", base::names(ph), ignore.case = TRUE, value = TRUE)[1]
  d$sex <- factor(base::toupper(base::substr(base::trimws(as.character(ph[[sx]])), 1, 1)))
  if (AGE_SRC[[g]] == "dnam") d$age <- cov$dnam_age_horvath else {
    ac <- if ("age_clean" %in% base::names(ph)) ph$age_clean else
      ph[[base::grep("^age$|ageat", base::names(ph), ignore.case = TRUE, value = TRUE)[1]]]
    d$age <- base::suppressWarnings(base::as.numeric(ac))
  }
  if (isTRUE(USE_SMOKING[[g]]) && "dnam_smoking" %in% base::names(cov))
    d$smoking <- cov$dnam_smoking
  d <- cbind(d, cov[, base::setdiff(base::grep("^ct_", base::names(cov), value = TRUE),
                                    DROP_CT), drop = FALSE])
  if ("scan_date" %in% base::names(ph)) {
    d$batch <- factor(as.character(ph$scan_date))
  } else {
    sl <- slide_from_map(g, gsmv)
    if (is.null(sl)) sl <- slide_from_hdr(g, gsmv)
    if (!is.null(sl) && base::mean(base::is.na(sl)) <= 0.05 &&
        length(base::unique(sl)) <= 0.25 * length(gsmv)) d$batch <- factor(sl)
  }
  if (!is.null(EXTRA[[g]])) for (tg in EXTRA[[g]]) {
    v <- get_char(g, tg, gsmv)
    if (!is.null(v)) { tb <- base::table(v)
      if (length(tb) >= 2L && length(tb) <= 12L) d[[tg]] <- factor(v) }
  }
  say("  Covariates: ", base::paste(base::setdiff(base::names(d), "group"), collapse = ", "))

  # Check against the saved design: the rebuilt model must be the one used by 37
  df <- file.path(RES_DIR, base::paste0(g, "_design.rds"))
  if (file.exists(df)) {
    od <- base::readRDS(df)$design
    n_old <- base::nrow(od); n_new <- base::sum(stats::complete.cases(d))
    hb_old <- base::any(base::grepl("^batch", base::colnames(od)))
    say("  Original design: samples ", n_old, ", parameters ", base::ncol(od),
        ", batch term ", if (hb_old) "yes" else "none")
    if (n_old != n_new || hb_old != ("batch" %in% base::names(d)))
      say("  [FAIL] differs from the original design; check the covariate sources before rebuilding")
  }

  say("  Computing principal components on the full matrix ...")
  t0 <- Sys.time()
  sc <- sample_pcs(M, K + 25L)
  say("  Time: ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  pg <- base::apply(sc, 2, function(v) stats::anova(stats::lm(v ~ d$group))$`Pr(>F)`[1])
  take <- utils::head(base::which(pg >= MPC_GROUP_P), K)
  say("  Components associated with group: ", base::sum(pg < MPC_GROUP_P), " removed; used: ", length(take), "")
  pcd <- base::as.data.frame(sc[, take, drop = FALSE])
  base::names(pcd) <- base::paste0("mpc", base::seq_along(take))
  d <- cbind(d, pcd); rm(sc); invisible(gc())

  ok <- stats::complete.cases(d)
  d <- d[ok, , drop = FALSE]; M <- M[, ok, drop = FALSE]
  des <- stats::model.matrix(~ ., data = d)
  q <- qr(des)
  if (q$rank < base::ncol(des))
    des <- des[, base::sort(q$pivot[base::seq_len(q$rank)]), drop = FALSE]
  tb <- base::table(d$group)
  say("  In the model: controls ", tb[["control"]], ", cases ", tb[["case"]],
      "; parameters ", base::ncol(des), ", residual degrees of freedom ", base::nrow(des) - base::ncol(des))

  fit <- limma::eBayes(limma::lmFit(M, des))
  tt <- limma::topTable(fit, coef = "groupcase", number = Inf, sort.by = "none")
  res <- data.frame(probe = base::rownames(tt), effect = tt$logFC,
                    se = fit$stdev.unscaled[, "groupcase"] * base::sqrt(fit$s2.post),
                    t = tt$t, p = tt$P.Value, stringsAsFactors = FALSE)
  lam <- lambda_of(res$p)
  say("  lambda = ", sprintf("%.3f", lam), "; determined by 37: ",
      sprintf("%.3f", todo$lambda_new[i]),
      ", difference ", sprintf("%.4f", base::abs(lam - todo$lambda_new[i])))
  if (base::abs(lam - todo$lambda_new[i]) > 0.02)
    say("  [NOTE] differs from the value determined by 37; the models of the two steps may differ")
  if (lam > GC_CUT) {
    cs <- stats::qchisq(res$p, 1, lower.tail = FALSE) / lam
    res$p_adj_scale <- stats::pchisq(cs, 1, lower.tail = FALSE); gcf <- TRUE
  } else { res$p_adj_scale <- res$p; gcf <- FALSE }
  res$fdr <- stats::p.adjust(res$p_adj_scale, method = "BH")
  say("  Genomic control: ", if (gcf) "applied" else "not applied",
      "; probes with FDR < 0.05: ", base::sum(res$fdr < 0.05, na.rm = TRUE))

  base::saveRDS(res, file.path(RES_DIR, base::paste0(g, "_ewas_final.rds")))
  base::saveRDS(list(design = des, gsm = base::colnames(M), probes = res$probe,
                     group = d$group, gse = g),
                file.path(RES_DIR, base::paste0(g, "_design.rds")))
  write.csv(utils::head(res[base::order(res$p), ], 200),
            file.path(OUT_DIR, base::paste0(g, "_top200.csv")), row.names = FALSE)
  say("  Rebuilt and overwritten: results and design matrix of ", g, "")
  rm(M, fit, des); invisible(gc())
}

say("")
say("=== Final settings of the ten cohorts ===")
fin <- SEL[, c("gse", "round", "K_new", "lambda_new", "genomic_control")]
fin$disease <- base::unlist(DISEASE[fin$gse])
print(fin[base::order(fin$disease, fin$gse), ])
write.csv(fin, file.path(OUT_DIR, "final_model_table.csv"), row.names = FALSE)
say("")
say("The per-cohort analyses of the ten sets are final; pairwise comparisons can start.")
log_stop()
