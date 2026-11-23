# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 35: model selection for the six additional cohorts (choice of K and inflation factor)
#
# Same method as step 15 of the first round; the rule follows decision D12, fixed beforehand and not changed after seeing results:
#   among K = 0, 5, 10, 20, 30, 50 take the smallest value that brings the inflation factor into 0.95 to 1.10;
#   if none does, take the value that brings it closest to 1 and apply genomic control in the formal analysis
#   according to D13. Methylation principal components associated with case-control status at p < 0.01
#   are never modelled.
#
# Covariates follow decision D4 together with what was measured in this round:
#   sex, age (source as in the table of 34b; three sets use DNAm age),
#   eleven cell fractions (Baso dropped), DNAm smoking score, batch.
#   Batch is the array ID; the field holding the barcode differs between cohorts and is found by the shared scanner;
#   cohorts without barcodes have no batch term, which is recorded as it is (as for GSE72774 in the first round).
#
# One point specific to this round: GSE338167 spans five population groups (Hispanic, African American, non-Hispanic White,
# Peruvian, Cuban); population group and sub-cohort must be in the model, otherwise population stratification would be taken for a disease effect.
#
# This step estimates the inflation factor on a random subset of fifty thousand probes only to choose K; the formal analysis recomputes
# on the full matrix (as in the first round).
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/35_model_selection_panel2.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "35_model_selection_panel2"

RAW_DIR <- file.path(DATA_ROOT, "raw")
SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
OUT_DIR <- file.path(PROJ_ROOT, "output", "35_model")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(limma))

set.seed(20260917)
N_PROBE     <- 50000L
K_GRID      <- c(0L, 5L, 10L, 20L, 30L, 50L)
N_MPC_MAX   <- 70L
MPC_GROUP_P <- 0.01
BETA_CLIP   <- 0.001
DROP_CT     <- "ct_Baso"
LAM_LO      <- 0.95
LAM_HI      <- 1.10

COHORTS <- list(
  GSE84727  = list(disease = "schizophrenia", age = "dnam"),
  GSE80417  = list(disease = "schizophrenia", age = "chrono"),
  GSE198904 = list(disease = "depression",    age = "chrono"),
  GSE273970 = list(disease = "depression",    age = "dnam"),
  GSE338167 = list(disease = "alzheimer",     age = "chrono",
                   extra = c("ethnicity", "cohort")),
  GSE144858 = list(disease = "alzheimer",     age = "chrono")
)

lambda_of <- function(p) {
  cs <- stats::qchisq(p, df = 1, lower.tail = FALSE)
  stats::median(cs, na.rm = TRUE) / stats::qchisq(0.5, 1)
}
# Find the header field carrying barcodes and take the array ID as batch
get_slide <- function(g, gsm_want) {
  hd <- mr_header(RAW_DIR, g)
  gsm <- mr_split(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  best <- NULL
  for (t in base::unique(base::sub("\t.*$", "", base::grep("^!Sample_", hd, value = TRUE)))) {
    for (one in base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)) {
      v <- mr_split(one)[-1]
      if (length(v) != length(gsm)) next
      hit <- base::regexpr(MR_BARCODE, v) > 0
      if (base::mean(hit) > 0.95) {
        bc <- base::rep(NA_character_, length(v))
        bc[hit] <- base::regmatches(v, base::regexpr(MR_BARCODE, v))
        best <- list(field = base::sub("^!Sample_", "", t),
                     slide = base::sub("_R[0-9]{2}C[0-9]{2}$", "", bc))
        break
      }
    }
    if (!is.null(best)) break
  }
  if (is.null(best)) return(NULL)
  list(field = best$field, slide = best$slide[base::match(gsm_want, gsm)])
}
# Take any characteristic field from the header
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

res <- list(); chosen <- list()
for (g in names(COHORTS)) {
  C <- COHORTS[[g]]
  say("==== ", g, " (", C$disease, ") ====")
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, base::paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[base::match(base::colnames(obj$beta), cov$gsm), , drop = FALSE]
  gsmv <- base::colnames(obj$beta)

  d <- data.frame(group = factor(as.character(ph$group), levels = c("control", "case")))
  d$sex <- factor(base::toupper(base::substr(base::trimws(as.character(ph$sex)), 1, 1)))
  d$age <- if (C$age == "dnam") cov$dnam_age_horvath else
    base::suppressWarnings(base::as.numeric(ph$age_clean))
  say("  Age source: ", if (C$age == "dnam") "DNAm age" else "chronological age",
      ", missing ", base::sum(base::is.na(d$age)), " samples")
  d$smoking <- cov$dnam_smoking
  cts <- cov[, base::setdiff(base::grep("^ct_", base::names(cov), value = TRUE), DROP_CT),
             drop = FALSE]
  d <- cbind(d, cts)

  sl <- get_slide(g, gsmv)
  if (is.null(sl) || base::mean(base::is.na(sl$slide)) > 0.05) {
    say("  Batch: no array barcode found in the header; this cohort has no batch term")
  } else {
    nlev <- length(base::unique(sl$slide))
    if (nlev > 0.25 * length(gsmv)) {
      say("  Batch: arrays ", nlev, "; fewer than 4 samples per array on average, a fixed effect would exhaust the degrees of freedom, not included")
    } else {
      d$batch <- factor(sl$slide)
      say("  Batch: taken from '", sl$field, "', arrays ", nlev, "")
    }
  }
  if (!is.null(C$extra)) {
    for (tg in C$extra) {
      v <- get_char(g, tg, gsmv)
      if (is.null(v)) { say("  [NOTE] field not found: '", tg, "'"); next }
      tb <- base::table(v)
      if (length(tb) < 2L || length(tb) > 12L) { say("  Field '", tg, "': number of levels unsuitable for modelling, skipped"); next }
      d[[tg]] <- factor(v)
      say("  Additional covariate '", tg, "': ",
          base::paste(sprintf("%s=%d", base::names(tb), base::as.integer(tb)), collapse = "; "))
    }
  }

  idx <- base::sample(base::seq_len(base::nrow(obj$beta)), base::min(N_PROBE, base::nrow(obj$beta)))
  b <- obj$beta[idx, , drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); M <- M[stats::complete.cases(M), , drop = FALSE]
  rm(b); invisible(gc())
  say("  Probes used to estimate the inflation factor: ", base::nrow(M), ", samples ", base::ncol(M))

  pr <- stats::prcomp(base::t(M), center = TRUE, scale. = FALSE)
  nmp <- base::min(N_MPC_MAX, base::ncol(pr$x))
  mpc <- pr$x[, base::seq_len(nmp), drop = FALSE]
  pg <- base::apply(mpc, 2, function(v) stats::anova(stats::lm(v ~ d$group))$`Pr(>F)`[1])
  ok_pc <- base::which(pg >= MPC_GROUP_P)
  say("  Candidate principal components: ", nmp, ", associated with group: ", base::sum(pg < MPC_GROUP_P),
      " removed, usable ", length(ok_pc), "")
  rm(pr); invisible(gc())

  lam <- rep(NA_real_, length(K_GRID))
  for (i in base::seq_along(K_GRID)) {
    K <- K_GRID[i]
    dd <- d
    if (K > 0) {
      take <- utils::head(ok_pc, K)
      if (length(take)) {
        pcd <- base::as.data.frame(mpc[, take, drop = FALSE])
        base::names(pcd) <- base::paste0("mpc", base::seq_along(take))
        dd <- cbind(dd, pcd)
      }
    }
    keep <- stats::complete.cases(dd)
    des <- stats::model.matrix(~ ., data = dd[keep, , drop = FALSE])
    # Collinear columns are removed automatically: in some cohorts sub-cohort, population group and array are nested in each other
    # and the model is necessarily rank-deficient with all of them. The pivot order of the QR identifies redundant columns, which are removed,
    # and what is removed is printed - it must not be dropped silently, because it is information about the model specification.
    q <- qr(des)
    if (q$rank < base::ncol(des)) {
      keepc <- base::sort(q$pivot[base::seq_len(q$rank)])
      dropped <- base::colnames(des)[base::setdiff(base::seq_len(base::ncol(des)), keepc)]
      if (K == K_GRID[1])
        say("  Collinear columns removed: ", length(dropped), "",
            base::paste(utils::head(dropped, 6), collapse = ", "),
            if (length(dropped) > 6) " ..." else "")
      des <- des[, keepc, drop = FALSE]
    }
    rdf <- base::nrow(des) - base::ncol(des)
    if (rdf < 30L) {
      say(sprintf("  K=%2d  parameters %3d  residual df %4d  insufficient degrees of freedom, skipped",
                  K, base::ncol(des), rdf)); next
    }
    fit <- limma::eBayes(limma::lmFit(M[, keep, drop = FALSE], des))
    lam[i] <- lambda_of(fit$p.value[, "groupcase"])
    say(sprintf("  K=%2d  parameters %3d  residual df %4d  lambda = %.3f",
                K, base::ncol(des), rdf, lam[i]))
  }

  if (base::all(base::is.na(lam))) {
    say("  [FAIL] no K can be fitted; this cohort is skipped and the covariate specification needs manual inspection")
    say("")
    rm(M, mpc); invisible(gc())
    next
  }
  inband <- base::which(!base::is.na(lam) & lam >= LAM_LO & lam <= LAM_HI)
  if (length(inband)) {
    pick <- K_GRID[inband[1]]; gc_flag <- FALSE
    say("  By D12, smallest K falling into ", LAM_LO, " to ", LAM_HI, ": K = ", pick,
        ", lambda = ", sprintf("%.3f", lam[inband[1]]))
  } else {
    j <- base::which.min(base::abs(lam - 1))
    pick <- K_GRID[j]; gc_flag <- TRUE
    say("  None falls into the interval; K closest to 1: K = ", pick,
        ", lambda = ", sprintf("%.3f", lam[j]), "; genomic control is applied in the formal analysis according to D13")
  }
  chosen[[length(chosen) + 1]] <- data.frame(
    gse = g, disease = C$disease, age_source = C$age,
    has_batch = "batch" %in% base::names(d), K = pick,
    lambda = lam[base::match(pick, K_GRID)], genomic_control = gc_flag,
    stringsAsFactors = FALSE)
  for (i in base::seq_along(K_GRID))
    res[[length(res) + 1]] <- data.frame(gse = g, disease = C$disease,
                                         K = K_GRID[i], lambda = lam[i],
                                         stringsAsFactors = FALSE)
  rm(M, mpc); invisible(gc())
  say("")
}

rr <- base::do.call(rbind, res)
cc <- base::do.call(rbind, chosen)
write.csv(rr, file.path(OUT_DIR, "lambda_ladder_panel2.csv"), row.names = FALSE)
write.csv(cc, file.path(OUT_DIR, "chosen_K_panel2.csv"), row.names = FALSE)
say("=== Selected settings ===")
print(cc)

say("")
say("=== Plotting ===")
rr$panel <- factor(base::paste0(rr$gse, "\n", rr$disease),
                   levels = base::unique(base::paste0(rr$gse, "\n", rr$disease)))
p <- ggplot(rr[!base::is.na(rr$lambda), ], aes(K, lambda)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_line(colour = PAL[1], linewidth = 0.55) +
  geom_point(colour = PAL[1], size = 1.6) +
  geom_text(aes(label = sprintf("%.2f", lambda)), vjust = -0.6, size = 2.4) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = K_GRID) +
  labs(title = "Genomic inflation against the number of methylation PCs",
       subtitle = "All models include group, sex, age, cell composition, smoking and batch where available",
       x = "Number of methylation principal components",
       y = expression(lambda)) +
  theme_paper()
save_fig(p, "fig39_lambda_ladder_panel2", width = 10, height = 6)
log_stop()
