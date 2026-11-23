# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 25: the decisive step - stratify sites by their "correlation with cell composition"
#          and see whether cross-disease sharing lives only in the composition-related stratum
#
# Three facts from 24:
#   (1) The permutation null and the hypergeometric expectation are almost identical (0.6/2.5/59.6 vs
#       0.6/2.5/61.5). The earlier concern that "shared technical structure raises the null"
#       was wrong; the fold enrichments of 23 stand. This is acknowledged.
#   (2) PD and IBD are indeed enriched in the tail: 15 overlapping among the top 1000, P = 0.014;
#       234 among the top 5000, P = 0.008. The global Spearman correlation
#       cannot detect it, because two hundred true signals among four hundred thousand sites hardly change the correlation.
#   (3) But the overlapping sites of all pairs - including true within-disease replication - are extremely enriched for
#       probes correlated with the neutrophil fraction: median |r| of the overlap 0.40 to 0.50,
#       against an overall median of only 0.076. Enrichment for smoking and age is an order of magnitude smaller.
#
# This gives a hypothesis that can be decided: the cross-disease "shared signature" is in fact
# residual cell composition signal - even though eleven cell fractions were removed as covariates in the model,
# linear adjustment leaves a residual, and the tail sites are exactly those with the largest residual.
#
# Decision procedure: divide sites into five strata by the strength of their correlation with the neutrophil fraction and repeat
# the overlap test within each stratum. If the enrichment of cross-disease pairs appears only in the top stratum while within-disease
# replication holds in all strata, "shared signature = composition residual" is confirmed.
#
# A key design point: the composition axis is defined with the control samples of the schizophrenia cohort,
# which is unrelated to both PD and IBD, so that the data being compared do not define their own strata.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/25_composition_strata.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "25_composition_strata"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "25_strata")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 500L
N_STRATA <- 5L
TOP_FRAC <- 0.02          # Within each stratum the top 2% of sites of each side are intersected
BETA_CLIP <- 0.001
AXIS_COHORT <- "GSE152027"    # Its control samples define the composition axis; unrelated to both sides

PAIRS <- list(
  list(id = "PD x IBD",        a = "GSE145361", b = "GSE87648",  kind = "comorbid pair"),
  list(id = "PD x SCZ",        a = "GSE145361", b = "GSE152027", kind = "negative-control pair"),
  list(id = "IBD x SCZ",       a = "GSE87648",  b = "GSE152027", kind = "negative-control pair"),
  list(id = "IBD x IBD(paed)", a = "GSE87648",  b = "GSE112611", kind = "within-disease positive pair"))

ew <- list(); dz <- list()
for (g in c("GSE145361", "GSE87648", "GSE112611", "GSE152027")) {
  f1 <- file.path(RES_DIR, paste0(g, "_ewas_final.rds"))
  f2 <- file.path(RES_DIR, paste0(g, "_design.rds"))
  if (file.exists(f1)) ew[[g]] <- readRDS(f1)
  if (file.exists(f2)) dz[[g]] <- readRDS(f2)
}
say("Loaded ", length(ew), " result sets")

# --------------------------------------------------------------
say("=== Defining the composition axis with the control samples of ", AXIS_COHORT, " ===")
D <- dz[[AXIS_COHORT]]
obj <- readRDS(file.path(SET_DIR, paste0(AXIS_COHORT, "_analysis_filtered.rds")))
cov <- read.csv(file.path(COV_DIR, paste0(AXIS_COHORT, "_covariates.csv")),
                stringsAsFactors = FALSE)
cov <- cov[match(D$gsm, cov$gsm), , drop = FALSE]
ctrl <- which(D$group == "control")
say("  Control samples: ", length(ctrl), " samples")
b <- obj$beta[D$probes, D$gsm[ctrl], drop = FALSE]
rm(obj); invisible(gc())
b[!is.finite(b)] <- NA
b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
M <- log2(b / (1 - b)); rm(b); invisible(gc())
x <- cov$ct_Neu[ctrl]
ok <- is.finite(x)
xs <- scale(x[ok])[, 1]
Mm <- M[, ok, drop = FALSE]
Mm <- Mm - rowMeans(Mm)
rr <- as.numeric(Mm %*% xs) / sqrt(sum(xs^2) * rowSums(Mm^2))
axis_r <- data.frame(probe = D$probes, absr = abs(rr), stringsAsFactors = FALSE)
rm(M, Mm); invisible(gc())
say("  Distribution of |r|: median ", sprintf("%.3f", stats::median(axis_r$absr)),
    ", 80th percentile ", sprintf("%.3f", stats::quantile(axis_r$absr, 0.8)),
    ", maximum ", sprintf("%.3f", max(axis_r$absr)))
write.csv(axis_r, file.path(OUT_DIR, "probe_composition_axis.csv"), row.names = FALSE)
say("")

# --------------------------------------------------------------
make_state <- function(g) {
  Dg <- dz[[g]]
  o <- readRDS(file.path(SET_DIR, paste0(g, "_analysis_filtered.rds")))
  bb <- o$beta[Dg$probes, Dg$gsm, drop = FALSE]
  rm(o); invisible(gc())
  bb[!is.finite(bb)] <- NA
  bb[bb < BETA_CLIP] <- BETA_CLIP; bb[bb > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  MM <- log2(bb / (1 - bb)); rm(bb); invisible(gc())
  des <- Dg$design
  jg <- which(colnames(des) == "groupcase")
  X0 <- des[, -jg, drop = FALSE]
  q0 <- qr(X0)
  R <- qr.resid(q0, t(MM))
  rss0 <- colSums(R^2)
  rm(MM); invisible(gc())
  list(R = R, q0 = q0, rss0 = rss0, probes = Dg$probes,
       group = Dg$group, df = nrow(des) - ncol(des))
}
pvals <- function(S, gvec) {
  gres <- qr.resid(S$q0, gvec)
  den <- sum(gres^2)
  bhat <- as.numeric(crossprod(gres, S$R)) / den
  rss <- S$rss0 - bhat^2 * den
  se <- sqrt(pmax(rss, 1e-12) / S$df / den)
  2 * stats::pt(abs(bhat / se), S$df, lower.tail = FALSE)
}

st <- list(); rows <- list()
for (P in PAIRS) {
  A <- ew[[P$a]]; B <- ew[[P$b]]
  common <- intersect(intersect(A$probe, B$probe), axis_r$probe)
  ia <- match(common, A$probe); ib <- match(common, B$probe)
  pa <- A$p[ia]; pb <- B$p[ib]
  ar <- axis_r$absr[match(common, axis_r$probe)]
  br <- cut(ar, breaks = stats::quantile(ar, seq(0, 1, length.out = N_STRATA + 1)),
            include.lowest = TRUE, labels = paste0("Q", seq_len(N_STRATA)))
  say("==== ", P$id, " (", P$kind, ") common probes ", length(common), " ====")
  say("  Median |r| with the neutrophil fraction per stratum: ",
      paste(sprintf("%s %.3f", levels(br),
                    tapply(ar, br, stats::median)), collapse = "; "))

  if (is.null(st[[P$b]])) { say("  Residualising ", P$b, " ..."); st[[P$b]] <- make_state(P$b) }
  S <- st[[P$b]]
  jb <- match(common, S$probes)
  gnum <- as.numeric(S$group == "case")

  # Observed: within each stratum intersect the top TOP_FRAC of each side
  obs <- integer(N_STRATA); nsel <- integer(N_STRATA)
  idx_by <- split(seq_along(common), br)
  for (k in seq_len(N_STRATA)) {
    ii <- idx_by[[k]]
    n <- max(20L, floor(length(ii) * TOP_FRAC))
    nsel[k] <- n
    sa <- ii[order(pa[ii])][seq_len(n)]
    sb <- ii[order(pb[ii])][seq_len(n)]
    obs[k] <- length(intersect(sa, sb))
  }
  # Permutation: one complete P value vector per permutation, reused for the five strata
  permM <- matrix(0L, N_PERM, N_STRATA)
  for (i in seq_len(N_PERM)) {
    pp <- pvals(S, sample(gnum))[jb]
    for (k in seq_len(N_STRATA)) {
      ii <- idx_by[[k]]
      sa <- ii[order(pa[ii])][seq_len(nsel[k])]
      sb <- ii[order(pp[ii])][seq_len(nsel[k])]
      permM[i, k] <- length(intersect(sa, sb))
    }
  }
  for (k in seq_len(N_STRATA)) {
    mu <- mean(permM[, k])
    pe <- (1 + sum(permM[, k] >= obs[k])) / (N_PERM + 1)
    say(sprintf("  %s (median |r| %.3f, %5d per stratum): overlap %4d, null mean %6.1f, fold %5.2f, P = %.4f",
                levels(br)[k], stats::median(ar[idx_by[[k]]]), nsel[k],
                obs[k], mu, obs[k] / max(mu, 0.5), pe))
    rows[[length(rows) + 1]] <- data.frame(
      pair = P$id, kind = P$kind, stratum = levels(br)[k],
      median_absr = stats::median(ar[idx_by[[k]]]), n_sel = nsel[k],
      observed = obs[k], null_mean = mu, fold = obs[k] / max(mu, 0.5),
      p_perm = pe, stringsAsFactors = FALSE)
  }
  say("")
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "strata_overlap.csv"), row.names = FALSE)

say("=== Interpretation ===")
say("  If the enrichment of cross-disease pairs is concentrated in Q5 (the stratum most related to composition) ")
say("  while within-disease replication holds from Q1 to Q5, 'shared signature = composition residual' holds.")
say("  If cross-disease pairs are also enriched in Q1 and Q2, there is true sharing unrelated to composition.")
say("")

# --------------------------------------------------------------
say("=== Plotting ===")
tab$panel <- factor(tab$pair, levels = vapply(PAIRS, function(x) x$id, ""))
tab$stratum <- factor(tab$stratum, levels = paste0("Q", seq_len(N_STRATA)))
p <- ggplot(tab, aes(stratum, fold, fill = kind)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(p_perm < 0.01, sprintf("%.1fx\nP<0.01", fold),
                               sprintf("%.1fx\nP=%.2f", fold, p_perm))),
            vjust = -0.3, size = 2.1) +
  facet_wrap(~ panel, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("comorbid pair" = PAL[2], "negative-control pair" = PAL[8],
                               "within-disease positive pair" = PAL[3]), name = NULL) +
  labs(title = "Where the apparent sharing lives",
       subtitle = "Probes split into five strata by their correlation with neutrophil fraction, measured in an unrelated cohort's controls; Q5 is the most composition-linked",
       x = "Stratum of correlation with cell composition",
       y = "Overlap divided by the permutation mean") +
  theme_paper()
save_fig(p, "fig34_strata_overlap", width = 11, height = 4.8)

log_stop()
