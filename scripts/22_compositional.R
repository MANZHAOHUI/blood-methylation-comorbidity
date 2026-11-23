# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 22: redo axis two with the correct approach for compositional data
#
# The first question of 21 was not really answered because the test design had a flaw:
# the twelve fractions always sum to 1. Once neutrophils rise, the sum of the other eleven must
# fall - this is a constraint, not a finding. Modelling each type with a logit treats every type
# as an independent variable and does not remove the constraint at all; "still 12/12 in the same direction
# after removing neutrophils" therefore does not prove sharing beyond the neutrophil effect, because the remaining eleven types
# still carry the shadow of the neutrophil increase.
#
# The evidence is in the numbers of 21: after removing neutrophils and eosinophils the sign concordance is still
# 1.00, but the Spearman correlation drops from 0.483 to 0.103 and is 0.000 with lymphocyte subsets only.
# In other words only the signs agree while the magnitudes are unrelated - exactly the
# artefact produced by the constraint, not evidence of sharing.
#
# The correct approach is compositional data analysis:
#   (1) centred log-ratio (clr): the log of each type minus the mean log of all types
#       of the sample. The effects then sum to zero over the twelve types, "all in the same direction"
#       becomes impossible by construction, and sign concordance is again a criterion that can fail.
#   (2) subcomposition: take the eight lymphocyte types alone and renormalise; the question becomes
#       "within lymphocytes, does the relative composition change in the same direction in both diseases",
#       which is entirely unrelated to neutrophils.
#
# The crash in part 3 of 21 is fixed as well: build() received a fraction matrix aligned to the already
# filtered samples, while the function applied complete.cases to the unfiltered samples,
# so the lengths differed. The fraction matrix is now named by GSM and aligned inside the function.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/22_compositional.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "22_compositional"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RDS_DIR <- file.path(DATA_ROOT, "rds")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "22_compositional")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 10000L
RUN_REF7 <- TRUE
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
LYMPH <- c("CD4Tnv", "CD4Tmem", "CD8Tnv", "CD8Tmem", "Treg", "Bnv", "Bmem", "NK")

ROLE <- c(GSE145361 = "PD discovery", GSE111629 = "PD replication",
          GSE72774 = "PD sensitivity", GSE87648 = "IBD main",
          GSE112611 = "IBD paediatric", GSE152027 = "Negative control")
ORDER <- names(ROLE)
AGE_SOURCE <- c(GSE145361 = "dnam", GSE152027 = "dnam", GSE111629 = "chrono",
                GSE72774 = "chrono", GSE87648 = "chrono", GSE112611 = "chrono")
USE_SMOKING <- c(GSE145361 = TRUE, GSE152027 = TRUE, GSE111629 = TRUE,
                 GSE72774 = TRUE, GSE87648 = TRUE, GSE112611 = FALSE)
PAIRS <- list(
  list(id = "PC1", kind = "positive control", a = "GSE145361", b = "GSE111629"),
  list(id = "PC2", kind = "positive control", a = "GSE87648",  b = "GSE112611"),
  list(id = "M1",  kind = "primary criterion",   a = "GSE145361", b = "GSE87648"),
  list(id = "M2",  kind = "secondary",     a = "GSE111629", b = "GSE87648"),
  list(id = "NC1", kind = "negative control", a = "GSE145361", b = "GSE152027"),
  list(id = "NC2", kind = "negative control", a = "GSE87648",  b = "GSE152027"))

# Centred log-ratio. A very small amount is added to zeros first, otherwise the log diverges.
clr <- function(P) {
  P <- pmax(P, 1e-5)
  L <- log(P)
  L - rowMeans(L)
}
slide_of <- function(g, gsm) {
  mf <- file.path(MAP_DIR, paste0(g, "_sample_map.csv"))
  if (!file.exists(mf)) return(rep(NA_character_, length(gsm)))
  m <- read.csv(mf, stringsAsFactors = FALSE, colClasses = "character")
  bc <- rep(NA_character_, nrow(m))
  hit <- regexpr(BARCODE, m$col_name) > 0
  bc[hit] <- regmatches(m$col_name, regexpr(BARCODE, m$col_name))
  sub("_R[0-9]{2}C[0-9]{2}$", "", bc)[match(gsm, m$gsm)]
}

# prop, if given, must be a matrix with GSM as row names; alignment is done inside the function
build <- function(g, prop = NULL, subset_cells = NULL) {
  obj <- readRDS(file.path(SET_DIR, paste0(g, "_analysis_filtered.rds")))
  gsm <- colnames(obj$beta); ph <- obj$pheno
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(gsm, cov$gsm), , drop = FALSE]
  d <- data.frame(group = factor(as.character(ph$group),
                                 levels = c("control", "case")))
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sx)) d$sex <- factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  if (AGE_SOURCE[g] == "dnam") d$age <- cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    d$age <- suppressWarnings(as.numeric(ph[[ag]]))
  }
  if (isTRUE(USE_SMOKING[g])) d$smoking <- cov$dnam_smoking
  if ("scan_date" %in% names(ph)) d$batch <- factor(as.character(ph$scan_date))
  else {
    sl <- slide_of(g, gsm)
    if (sum(!is.na(sl)) > 0.95 * length(sl) &&
        length(unique(sl)) <= 0.25 * length(sl)) d$batch <- factor(sl)
  }
  if (is.null(prop)) {
    ct <- grep("^ct_", names(cov), value = TRUE)
    P <- as.matrix(cov[, ct, drop = FALSE]); colnames(P) <- sub("^ct_", "", ct)
    rownames(P) <- gsm
  } else {
    P <- prop[match(gsm, rownames(prop)), , drop = FALSE]
    rownames(P) <- gsm
  }
  if (!is.null(subset_cells)) {
    P <- P[, intersect(subset_cells, colnames(P)), drop = FALSE]
    P <- P / rowSums(P)                       # Renormalise the subcomposition
  }
  Y <- clr(P)
  ok <- stats::complete.cases(d) & stats::complete.cases(Y)
  list(d = d[ok, , drop = FALSE], Y = Y[ok, , drop = FALSE], gsm = gsm[ok])
}
coef_vec <- function(B, perm = FALSE) {
  d <- B$d
  if (perm) d$group <- sample(d$group)
  X <- stats::model.matrix(~ ., data = d)
  j <- which(colnames(X) == "groupcase")
  out <- vapply(seq_len(ncol(B$Y)), function(k)
    stats::lm.fit(X, B$Y[, k])$coefficients[j], numeric(1))
  names(out) <- colnames(B$Y); out
}
test_pair <- function(Bb, ea, eb) {
  keep <- intersect(names(ea), names(eb))
  a <- ea[keep]; b <- eb[keep]
  rho <- stats::cor(a, b, method = "spearman")
  dir <- mean(sign(a) == sign(b))
  pr <- numeric(N_PERM); pd <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    bp <- coef_vec(Bb, perm = TRUE)[keep]
    pr[i] <- stats::cor(a, bp, method = "spearman")
    pd[i] <- mean(sign(a) == sign(bp))
  }
  list(rho = rho, dir = dir, n = length(keep),
       p_rho = (1 + sum(pr >= rho)) / (N_PERM + 1),
       p_dir = (1 + sum(pd >= dir)) / (N_PERM + 1),
       z = (rho - mean(pr)) / stats::sd(pr), null = pr)
}

run_block <- function(label, subset_cells = NULL, prop_list = NULL) {
  say("=== ", label, " ===")
  B <- list(); E <- list()
  for (g in names(ROLE)) {
    pr <- if (is.null(prop_list)) NULL else prop_list[[g]]
    if (!is.null(prop_list) && is.null(pr)) next
    B[[g]] <- build(g, prop = pr, subset_cells = subset_cells)
    E[[g]] <- coef_vec(B[[g]])
  }
  say("  Dimensions ", ncol(B[[1]]$Y), ": ", paste(colnames(B[[1]]$Y), collapse = ", "))
  say("  Note: after clr the effects of each data set sum to zero, ")
  say("        so 'all in the same direction' is impossible by construction and sign concordance is again a criterion that can fail")
  rows <- list(); nulls <- list()
  for (P in PAIRS) {
    if (is.null(E[[P$a]]) || is.null(E[[P$b]])) next
    r <- test_pair(B[[P$b]], E[[P$a]], E[[P$b]])
    say(sprintf("  %-4s %-8s %s vs %s  n=%2d  rho=%+.3f  sign %.2f  P(rho)=%.4f  P(sign)=%.4f  Z=%+.2f",
                P$id, P$kind, P$a, P$b, r$n, r$rho, r$dir, r$p_rho, r$p_dir, r$z))
    rows[[length(rows) + 1]] <- data.frame(
      block = label, pair = P$id, kind = P$kind, a = P$a, b = P$b, n = r$n,
      rho = r$rho, dir = r$dir, p_rho = r$p_rho, p_dir = r$p_dir, z = r$z,
      stringsAsFactors = FALSE)
    nulls[[P$id]] <- data.frame(block = label, pair = P$id, kind = P$kind,
                                rho = r$null, obs = r$rho, stringsAsFactors = FALSE)
  }
  say("")
  list(tab = do.call(rbind, rows), eff = E, null = do.call(rbind, nulls))
}

# --------------------------------------------------------------
r1 <- run_block("twelve types clr")
r2 <- run_block("lymphocyte subcomposition clr", subset_cells = LYMPH)

# --------------------------------------------------------------
prop7 <- NULL
if (RUN_REF7) {
  say("=== Re-estimation with the seven-cell-type reference (does the conclusion depend on the reference matrix?) ===")
  suppressPackageStartupMessages(library(EpiDISH))
  e <- new.env(); utils::data("centDHSbloodDMC.m", package = "EpiDISH", envir = e)
  ref7 <- get("centDHSbloodDMC.m", envir = e)
  prop7 <- list()
  for (g in names(ROLE)) {
    bf <- file.path(RDS_DIR, paste0(g, "_beta.rds"))
    mf <- file.path(MAP_DIR, paste0(g, "_sample_map.csv"))
    if (!file.exists(bf) || !file.exists(mf)) { say("  ", g, ": raw matrix missing, skipped"); next }
    m <- read.csv(mf, stringsAsFactors = FALSE, colClasses = "character")
    bm <- readRDS(bf)
    bm <- bm[, as.integer(m$col_pos_in_data) - 1L, drop = FALSE]
    colnames(bm) <- m$gsm
    est <- EpiDISH::epidish(beta.m = bm, ref.m = ref7, method = "RPC")$estF
    rm(bm); invisible(gc())
    prop7[[g]] <- est
    say("  ", g, " finished, ", ncol(est), " types")
  }
  say("")
}
r3 <- if (!is.null(prop7) && length(prop7) >= 2)
  run_block("seven-type reference clr", prop_list = prop7) else NULL

# --------------------------------------------------------------
tabs <- Filter(Negate(is.null), list(r1$tab, r2$tab, if (!is.null(r3)) r3$tab))
alltab <- do.call(rbind, tabs)
write.csv(alltab, file.path(OUT_DIR, "clr_concordance.csv"), row.names = FALSE)

say("=== Primary criterion and controls compared across the three analyses ===")
for (b in unique(alltab$block)) {
  s <- alltab[alltab$block == b, ]
  say("  ", b)
  for (i in seq_len(nrow(s)))
    say(sprintf("     %-4s %-8s rho=%+.3f sign %.2f  P(rho)=%.4f",
                s$pair[i], s$kind[i], s$rho[i], s$dir[i], s$p_rho[i]))
}

# --------------------------------------------------------------
say("")
say("=== Plotting ===")
E1 <- r1$eff
M <- do.call(rbind, lapply(names(E1), function(g)
  data.frame(gse = g, cell = names(E1[[g]]), effect = as.numeric(E1[[g]]),
             stringsAsFactors = FALSE)))
M$panel <- factor(paste0(M$gse, "\n", ROLE[M$gse]),
                  levels = paste0(ORDER, "\n", ROLE[ORDER]))
lim <- max(abs(M$effect))
p <- ggplot(M, aes(cell, panel, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = PAL[1], mid = "grey96", high = PAL[2], midpoint = 0,
                       limits = c(-lim, lim), name = "Case minus control\n(clr scale)") +
  labs(title = "Composition effects on the centred log-ratio scale",
       subtitle = "Effects now sum to zero within each cohort, so a shared direction cannot arise from the unit-sum constraint",
       x = NULL, y = NULL) +
  theme_paper() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
save_fig(p, "fig28_clr_effects", width = 9, height = 4.4)

alltab$panel <- factor(paste0(alltab$pair, " ", alltab$kind),
                       levels = unique(paste0(alltab$pair, " ", alltab$kind)))
p <- ggplot(alltab, aes(panel, rho, fill = kind)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("P=%.3f", p_rho)), vjust = -0.4, size = 2.3) +
  facet_wrap(~ block, ncol = 1) +
  scale_fill_manual(values = c("positive control" = PAL[3], "primary criterion" = PAL[2],
                               "secondary" = PAL[5], "negative control" = PAL[8]), name = NULL) +
  labs(title = "Concordance after proper compositional treatment",
       x = NULL, y = "Spearman correlation of clr effects") +
  theme_paper()
save_fig(p, "fig29_clr_concordance", width = 8.5, height = 7)

log_stop()
