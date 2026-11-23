# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 37: re-determine K on the full matrix
#
# Background: the self-check of 36 exposed a flaw in the practice of "choosing K on a subset and analysing on the full matrix".
# Listing the inflation factors of the two by K makes the pattern obvious:
#     two sets with K= 0    difference 0.013 and 0.013
#     two sets with K= 5    difference 0.013 and 0.011
#     one set with K=20     difference 0.105
#     one set with K=50     difference 0.112
# At K=0 the two are the same model and differ only by sampling noise; at K>0 the principal components computed on the fifty-thousand-probe
# subset and on the full matrix are not the same vectors, so the models differ. The larger K,
# the less reliable the lambda estimated on the subset, and the cohorts needing a large K depend on this judgement the most.
#
# Concrete consequence: on the subset GSE84727 was judged "none within 0.95 to 1.10" and K=20 was taken;
# recomputed on the full matrix, lambda at K=20 is 1.099, within the interval. Lambda at K=0 on the full matrix
# had never been computed - if it is also within the interval, D12 requires K=0 and 20 degrees of freedom were spent for nothing.
#
# This step redoes the choice as D12 intends: estimate the inflation factor on the full matrix for each K and take the smallest K
# within the interval. K=0 is computed first; if it is within the interval the choice is settled and no principal component is needed;
# otherwise the principal components are computed and K is increased step by step.
#
# The six sets of the first round are checked as well: their subset and full-matrix values differ by at most 0.031 and the choice is expected to stand,
# but all ten sets should go through the same procedure rather than two rounds with two standards.
#
# This step only determines K and changes no existing result. Cohorts with unchanged K need no new EWAS;
# changed ones are listed explicitly in the log and rebuilt separately.
#
# Environment: Windows + RStudio. The slowest is GSE145361 (1638 samples)
# Location: <project root>/scripts/37_K_reselect_full.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "37_K_reselect_full"

RAW_DIR <- file.path(DATA_ROOT, "raw")
SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "37_reselect")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260917)
K_GRID      <- c(0L, 5L, 10L, 20L, 30L, 50L)
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
MPC_GROUP_P <- 0.01
LAM_LO      <- 0.95
LAM_HI      <- 1.10

# Current settings: first round from steps 15/17, second round from steps 35/36
CUR <- list(
  GSE145361 = list(K = 30, age = "dnam",   round = 1),
  GSE111629 = list(K = 0,  age = "chrono", round = 1),
  GSE87648  = list(K = 5,  age = "chrono", round = 1),
  GSE112611 = list(K = 0,  age = "chrono", round = 1),
  GSE84727  = list(K = 20, age = "dnam",   round = 2),
  GSE80417  = list(K = 5,  age = "chrono", round = 2),
  GSE198904 = list(K = 50, age = "chrono", round = 2),
  GSE273970 = list(K = 0,  age = "dnam",   round = 2),
  GSE338167 = list(K = 0,  age = "chrono", round = 2, extra = c("ethnicity", "cohort")),
  GSE144858 = list(K = 5,  age = "chrono", round = 2)
)

# Whether the smoking score is included (as in the first round)
USE_SMOKING <- list(
  GSE145361 = TRUE,  GSE111629 = TRUE,  GSE87648  = TRUE,  GSE112611 = FALSE,
  GSE84727  = TRUE,  GSE80417  = TRUE,  GSE198904 = TRUE,  GSE273970 = TRUE,
  GSE338167 = TRUE,  GSE144858 = TRUE)

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
  base::rownames(sc) <- base::colnames(M); base::colnames(sc) <- base::paste0("PC", base::seq_len(k))
  sc
}
fit_lambda <- function(M, d) {
  keep <- stats::complete.cases(d)
  des <- stats::model.matrix(~ ., data = d[keep, , drop = FALSE])
  q <- qr(des)
  if (q$rank < base::ncol(des))
    des <- des[, base::sort(q$pivot[base::seq_len(q$rank)]), drop = FALSE]
  rdf <- base::nrow(des) - base::ncol(des)
  if (rdf < 30L) return(list(lam = NA_real_, npar = base::ncol(des), rdf = rdf))
  fit <- limma::eBayes(limma::lmFit(M[, keep, drop = FALSE], des))
  list(lam = lambda_of(fit$p.value[, "groupcase"]), npar = base::ncol(des), rdf = rdf)
}
# Barcode source of the first round: the sample mapping table generated by 06
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
get_slide <- function(g, gsm_want) {
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

ONLY <- character(0)   # Empty = all; cohort names = only these sets
rows <- list(); ladder <- list()
for (g in names(CUR)) {
  if (length(ONLY) && !(g %in% ONLY)) next
  C <- CUR[[g]]
  say("==== ", g, " (round ", C$round, ", current K = ", C$K, ") ====")
  f <- file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds"))
  if (!file.exists(f)) { say("  Analysis set not found, skipped"); next }
  obj <- base::readRDS(f); ph <- obj$pheno; gsmv <- base::colnames(obj$beta)
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
  if (C$age == "dnam") d$age <- cov$dnam_age_horvath else {
    ac <- if ("age_clean" %in% base::names(ph)) ph$age_clean else
      ph[[base::grep("^age$|ageat", base::names(ph), ignore.case = TRUE, value = TRUE)[1]]]
    d$age <- base::suppressWarnings(base::as.numeric(ac))
  }
  # The smoking score is decided per cohort: GSE112611 is a paediatric cohort and the score cannot be computed for 52 samples;
  # the term was already excluded in the first round. Without excluding it here, those 52 samples would be dropped for missingness,
  # and both model and sample set would differ from the original analysis.
  if (isTRUE(USE_SMOKING[[g]]) && "dnam_smoking" %in% base::names(cov))
    d$smoking <- cov$dnam_smoking else say("  Smoking score: not included for this cohort")
  d <- cbind(d, cov[, base::setdiff(base::grep("^ct_", base::names(cov), value = TRUE),
                                    DROP_CT), drop = FALSE])
  # Batch has two sources, and the order must not be reversed:
  #   the barcodes of the six first-round sets are in the sample mapping table generated by 06 (parsed from the matrix column names),
  #   those of the six second-round sets are in a field of the series matrix header.
  # Looking only in the header would lose the batch term for three first-round sets, and the model would no longer match 17.
  if ("scan_date" %in% base::names(ph)) {
    d$batch <- factor(as.character(ph$scan_date))
  } else {
    sl <- slide_from_map(g, gsmv)
    if (is.null(sl)) sl <- get_slide(g, gsmv)
    if (!is.null(sl) && base::mean(base::is.na(sl)) <= 0.05 &&
        length(base::unique(sl)) <= 0.25 * length(gsmv)) d$batch <- factor(sl)
  }
  if (!is.null(C$extra)) for (tg in C$extra) {
    v <- get_char(g, tg, gsmv)
    if (!is.null(v)) { tb <- base::table(v)
      if (length(tb) >= 2L && length(tb) <= 12L) d[[tg]] <- factor(v) }
  }
  say("  Covariates: ", base::paste(base::setdiff(base::names(d), "group"), collapse = ", "))
  # Comparison with the saved design matrix: the model rebuilt here must be the one used when the EWAS was run,
  # otherwise a change in K would be due to a different model rather than to the re-determination.
  df <- file.path(RES_DIR, base::paste0(g, "_design.rds"))
  if (file.exists(df)) {
    old_des <- base::readRDS(df)$design
    base_old <- base::unique(base::sub("^(batch|sex|ethnicity|cohort)[^ ]*$", "\\1",
                                       base::colnames(old_des)))
    has_batch_old <- base::any(base::grepl("^batch", base::colnames(old_des)))
    has_batch_new <- "batch" %in% base::names(d)
    n_mpc_old <- base::sum(base::grepl("^mpc", base::colnames(old_des)))
    n_old <- base::nrow(old_des)
    say("  Saved design: samples ", n_old, ", parameters ", base::ncol(old_des),
        ", of which principal components ", n_mpc_old, ", batch term ", if (has_batch_old) "yes" else "none")
    # All three must agree: presence of the batch term, number of samples in the model, number of parameters without principal components.
    # Checking batch alone would miss differences such as samples dropped for a missing smoking score.
    n_new <- base::sum(stats::complete.cases(d))
    npar_old_base <- base::ncol(old_des) - n_mpc_old
    if (has_batch_old != has_batch_new)
      say("  [FAIL] batch term differs: original ", has_batch_old, ", this step ", has_batch_new)
    if (n_old != n_new)
      say("  [FAIL] number of samples in the model differs: original ", n_old, ", this step ", n_new,
          ". Most likely samples were dropped because a covariate is missing")
    if (n_old == n_new && has_batch_old == has_batch_new)
      say("  Model comparison passed (sample number and batch term agree; expected parameters at K=0: ",
          npar_old_base, ")")
  }

  lam <- rep(NA_real_, length(K_GRID))
  r0 <- fit_lambda(M, d)
  lam[1] <- r0$lam
  say(sprintf("  K= 0  parameters %3d  residual df %4d  lambda = %.3f", r0$npar, r0$rdf, r0$lam))
  if (!base::is.na(r0$lam) && r0$lam >= LAM_LO && r0$lam <= LAM_HI) {
    say("  K=0 is within the interval; settled by D12, no principal components needed")
  } else {
    say("  K=0 is outside the interval; computing principal components on the full matrix ...")
    t0 <- Sys.time()
    sc <- sample_pcs(M, base::max(K_GRID) + 25L)
    say("  Time: ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
    pg <- base::apply(sc, 2, function(v) stats::anova(stats::lm(v ~ d$group))$`Pr(>F)`[1])
    avail <- base::which(pg >= MPC_GROUP_P)
    say("  Components associated with group: ", base::sum(pg < MPC_GROUP_P), " removed, usable ", length(avail), "")
    for (i in 2:length(K_GRID)) {
      take <- utils::head(avail, K_GRID[i])
      pcd <- base::as.data.frame(sc[, take, drop = FALSE])
      base::names(pcd) <- base::paste0("mpc", base::seq_along(take))
      r <- fit_lambda(M, cbind(d, pcd))
      lam[i] <- r$lam
      say(sprintf("  K=%2d  parameters %3d  residual df %4d  lambda = %.3f",
                  K_GRID[i], r$npar, r$rdf, r$lam))
      if (!base::is.na(r$lam) && r$lam >= LAM_LO && r$lam <= LAM_HI) {
        say("  within the interval; by D12 the smallest K is taken and no larger K is tried"); break
      }
    }
    rm(sc); invisible(gc())
  }

  inb <- base::which(!base::is.na(lam) & lam >= LAM_LO & lam <= LAM_HI)
  if (length(inb)) { newK <- K_GRID[inb[1]]; newL <- lam[inb[1]]; gcf <- FALSE
  } else { j <- base::which.min(base::abs(lam - 1)); newK <- K_GRID[j]; newL <- lam[j]; gcf <- TRUE }
  chg <- newK != C$K
  say("  Determined: K = ", newK, ", lambda = ", sprintf("%.3f", newL),
      ", genomic control ", if (gcf) "applied" else "not applied",
      if (chg) base::paste0("   [differs from the current K=", C$K, "; the EWAS of this cohort must be rebuilt]")
      else "   same as current, no rebuild needed")
  rows[[length(rows) + 1]] <- data.frame(
    gse = g, round = C$round, K_current = C$K, K_new = newK,
    lambda_new = newL, genomic_control = gcf, changed = chg, stringsAsFactors = FALSE)
  for (i in base::seq_along(K_GRID))
    ladder[[length(ladder) + 1]] <- data.frame(gse = g, K = K_GRID[i],
                                               lambda = lam[i], stringsAsFactors = FALSE)
  rm(M); invisible(gc())
  say("")
}

tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "K_reselect.csv"), row.names = FALSE)
write.csv(base::do.call(rbind, ladder), file.path(OUT_DIR, "lambda_ladder_full.csv"),
          row.names = FALSE)
say("=== Determined settings ===")
print(tab)
say("")
say("Cohorts whose EWAS must be rebuilt: ", base::sum(tab$changed), "",
    if (base::sum(tab$changed)) base::paste0(" (",
      base::paste(tab$gse[tab$changed], collapse = ", "), ")") else "")
log_stop()
