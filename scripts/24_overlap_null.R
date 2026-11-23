# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 24: give the overlap approach a real null distribution and find out what the overlap contains
#
# Two problems of 23 must be fixed first:
#
#   1   The hypergeometric expectation was used as the null - which is the very assumption that template papers rely on
#       and the thing to be questioned. It assumes that the sites selected in the two diseases are mutually
#       independent, whereas in practice both sides preferentially select the same sites with large variance and
#       technical instability, so the overlap is inflated by nature. Generating the null of the overlap
#       itself within the same permutation framework is the self-consistent approach. Criticising others while using
#       the same assumption was an error and is acknowledged here.
#
#   2   What the overlap contains was not examined. Among the 183 genes output by 23, MPO,
#       PRTN3 and CEACAM8 are neutrophil granule protein genes, FCER1G, CORO1A,
#       SELPLG, ICAM3, LY9 and FYB are immune cell identity genes, and ABCG1,
#       SREBF1 and FKBP5 are loci that appear in almost every blood EWAS.
#       This suggests that the "shared signature" is in fact a cell composition signature. This step quantifies it:
#       in control samples, compute the strength of association of every site with the neutrophil fraction, the smoking score and
#       age, then test whether the overlap is enriched for such sites.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/24_overlap_null.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "24_overlap_null"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "24_overlap")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 1000L
TOPN <- c(500L, 1000L, 5000L)
BETA_CLIP <- 0.001

ROLE <- c(GSE145361 = "PD", GSE87648 = "IBD",
          GSE112611 = "IBD (paediatric)", GSE152027 = "Schizophrenia")
PAIRS <- list(
  list(id = "PD x IBD",        a = "GSE145361", b = "GSE87648",  kind = "comorbid pair"),
  list(id = "PD x SCZ",        a = "GSE145361", b = "GSE152027", kind = "negative-control pair"),
  list(id = "IBD x SCZ",       a = "GSE87648",  b = "GSE152027", kind = "negative-control pair"),
  list(id = "IBD x IBD(paed)", a = "GSE87648",  b = "GSE112611", kind = "within-disease positive pair"))

ew <- list()
for (g in names(ROLE)) {
  f <- file.path(RES_DIR, paste0(g, "_ewas_final.rds"))
  if (file.exists(f)) ew[[g]] <- readRDS(f)
}
dz <- list()
for (g in names(ROLE)) {
  f <- file.path(RES_DIR, paste0(g, "_design.rds"))
  if (file.exists(f)) dz[[g]] <- readRDS(f)
}
say("Loaded ", length(ew), " result sets and ", length(dz), " design matrices")
say("")

# --------------------------------------------------------------
# Residualise, for fast recomputation of coefficients and P values during permutation
# --------------------------------------------------------------
make_state <- function(g) {
  D <- dz[[g]]
  obj <- readRDS(file.path(SET_DIR, paste0(g, "_analysis_filtered.rds")))
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- log2(b / (1 - b)); rm(b); invisible(gc())
  des <- D$design
  jg <- which(colnames(des) == "groupcase")
  X0 <- des[, -jg, drop = FALSE]
  q0 <- qr(X0)
  R <- qr.resid(q0, t(M))                 # samples x probes
  rss0 <- colSums(R^2)
  rm(M); invisible(gc())
  list(R = R, q0 = q0, rss0 = rss0, probes = D$probes,
       group = D$group, df = nrow(des) - ncol(des))
}
pvals <- function(st, gvec) {
  gres <- qr.resid(st$q0, gvec)
  den <- sum(gres^2)
  bhat <- as.numeric(crossprod(gres, st$R)) / den
  rss <- st$rss0 - bhat^2 * den
  se <- sqrt(pmax(rss, 1e-12) / st$df / den)
  2 * stats::pt(abs(bhat / se), st$df, lower.tail = FALSE)
}

st <- list()
rows <- list(); nulls <- list()
for (P in PAIRS) {
  if (is.null(ew[[P$a]]) || is.null(ew[[P$b]])) next
  A <- ew[[P$a]]; B <- ew[[P$b]]
  common <- intersect(A$probe, B$probe)
  ia <- match(common, A$probe); ib <- match(common, B$probe)
  pa <- A$p[ia]; pb <- B$p[ib]
  say("==== ", P$id, " (", P$kind, ")  common probes ", length(common), " ====")

  if (is.null(st[[P$b]])) {
    say("  Residualising ", P$b, "...")
    st[[P$b]] <- make_state(P$b)
  }
  S <- st[[P$b]]
  jb <- match(common, S$probes)
  gnum <- as.numeric(S$group == "case")
  chk <- stats::cor(pvals(S, gnum)[jb], pb, use = "complete.obs")
  say("  Correlation between P values recomputed by residualisation and limma results = ", sprintf("%.4f", chk))

  ordA <- order(pa)
  for (n in TOPN) {
    sa <- common[ordA][seq_len(n)]
    obs <- length(intersect(sa, common[order(pb)][seq_len(n)]))
    perm <- numeric(N_PERM)
    for (i in seq_len(N_PERM)) {
      pp <- pvals(S, sample(gnum))[jb]
      perm[i] <- length(intersect(sa, common[order(pp)][seq_len(n)]))
    }
    pe <- (1 + sum(perm >= obs)) / (N_PERM + 1)
    hyper <- n * n / length(common)
    say(sprintf("  top %5d: observed overlap %4d, permutation null mean %6.1f (hypergeometric expectation %5.1f), fold (vs permutation) %5.2f, P = %.4f",
                n, obs, mean(perm), hyper, obs / max(mean(perm), 0.5), pe))
    rows[[length(rows) + 1]] <- data.frame(
      pair = P$id, kind = P$kind, topn = n, observed = obs,
      null_mean = mean(perm), null_sd = stats::sd(perm),
      hyper = hyper, fold_perm = obs / max(mean(perm), 0.5),
      fold_hyper = obs / hyper, p_perm = pe, stringsAsFactors = FALSE)
    if (n == 5000L) nulls[[P$id]] <- data.frame(
      pair = P$id, kind = P$kind, overlap = perm, obs = obs,
      stringsAsFactors = FALSE)
  }
  say("")
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "overlap_permutation.csv"), row.names = FALSE)

# --------------------------------------------------------------
say("=== What the overlap contains: association with cell composition, smoking and age ===")
annotate_axis <- function(g) {
  D <- dz[[g]]
  obj <- readRDS(file.path(SET_DIR, paste0(g, "_analysis_filtered.rds")))
  cov <- read.csv(file.path(COV_DIR, paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[match(D$gsm, cov$gsm), , drop = FALSE]
  ctrl <- which(D$group == "control")          # Controls only, to avoid circularity
  b <- obj$beta[D$probes, D$gsm[ctrl], drop = FALSE]
  rm(obj); invisible(gc())
  b[!is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- log2(b / (1 - b)); rm(b); invisible(gc())
  ax <- list(neutrophil = cov$ct_Neu[ctrl], smoking = cov$dnam_smoking[ctrl],
             age = cov$dnam_age_horvath[ctrl])
  out <- list()
  for (k in names(ax)) {
    v <- ax[[k]]
    ok <- is.finite(v)
    if (sum(ok) < 30) { out[[k]] <- rep(NA_real_, nrow(M)); next }
    x <- scale(v[ok])[, 1]
    Mm <- M[, ok, drop = FALSE]
    Mm <- Mm - rowMeans(Mm, na.rm = TRUE)
    out[[k]] <- abs(as.numeric(Mm %*% x) / sqrt(sum(x^2) * rowSums(Mm^2, na.rm = TRUE)))
  }
  rm(M); invisible(gc())
  data.frame(probe = D$probes, neutrophil = out$neutrophil,
             smoking = out$smoking, age = out$age, stringsAsFactors = FALSE)
}
ax <- annotate_axis("GSE87648")     # Control samples of the IBD main set serve as the reference
say("  Control samples of the reference cohort GSE87648; per-site strength of correlation with the three axes")

enr <- list()
for (P in PAIRS) {
  if (is.null(ew[[P$a]]) || is.null(ew[[P$b]])) next
  A <- ew[[P$a]]; B <- ew[[P$b]]
  common <- intersect(A$probe, B$probe)
  ia <- match(common, A$probe); ib <- match(common, B$probe)
  n <- 5000L
  sa <- common[order(A$p[ia])][seq_len(n)]
  sb <- common[order(B$p[ib])][seq_len(n)]
  ov <- intersect(sa, sb)
  if (!length(ov)) next
  m <- match(ov, ax$probe); mb <- match(common, ax$probe)
  say("  ", P$id, " (overlap ", length(ov), ")")
  for (k in c("neutrophil", "smoking", "age")) {
    a1 <- ax[[k]][m]; a0 <- ax[[k]][mb]
    a1 <- a1[is.finite(a1)]; a0 <- a0[is.finite(a0)]
    if (!length(a1)) next
    w <- suppressWarnings(stats::wilcox.test(a1, a0))
    say(sprintf("     correlation with %-10s: overlap median %.3f, overall median %.3f, P = %s",
                switch(k, neutrophil = "neutrophil fraction", smoking = "smoking score", age = "age"),
                stats::median(a1), stats::median(a0),
                format.pval(w$p.value, digits = 3)))
    enr[[length(enr) + 1]] <- data.frame(
      pair = P$id, kind = P$kind, axis = k,
      med_overlap = stats::median(a1), med_all = stats::median(a0),
      p = w$p.value, stringsAsFactors = FALSE)
  }
}
en <- do.call(rbind, enr)
write.csv(en, file.path(OUT_DIR, "overlap_axis_enrichment.csv"), row.names = FALSE)

# --------------------------------------------------------------
say("")
say("=== Plotting ===")
tab$panel <- factor(tab$pair, levels = vapply(PAIRS, function(x) x$id, ""))
tab$topn_f <- factor(paste0("top", tab$topn), levels = paste0("top", TOPN))
p <- ggplot(tab, aes(topn_f, fold_perm, fill = kind)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.1fx\nP=%.3f", fold_perm, p_perm)),
            vjust = -0.3, size = 2.2) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c("comorbid pair" = PAL[2], "negative-control pair" = PAL[8],
                               "within-disease positive pair" = PAL[3]), name = NULL) +
  labs(title = "Overlap against a permutation null rather than a hypergeometric one",
       subtitle = "The permutation null preserves each cohort's probe-level variance structure, which the hypergeometric expectation ignores",
       x = NULL, y = "Observed overlap divided by the permutation mean") +
  theme_paper()
save_fig(p, "fig32_overlap_permutation", width = 10, height = 4.8)

if (!is.null(en)) {
  en$panel <- factor(en$pair, levels = vapply(PAIRS, function(x) x$id, ""))
  en$axis <- factor(en$axis, levels = c("neutrophil", "smoking", "age"),
                    labels = c("Neutrophil fraction", "Smoking score", "Age"))
  ln <- rbind(data.frame(en[, c("panel", "axis")], what = "Overlap set",
                         value = en$med_overlap),
              data.frame(en[, c("panel", "axis")], what = "All probes",
                         value = en$med_all))
  p <- ggplot(ln, aes(axis, value, fill = what)) +
    geom_col(position = position_dodge(0.75), width = 0.68) +
    facet_wrap(~ panel, nrow = 1) +
    scale_fill_manual(values = c("Overlap set" = PAL[2], "All probes" = "grey70"),
                      name = NULL) +
    labs(title = "What the overlapping probes are associated with",
         subtitle = "Absolute correlation with each axis, measured in control samples only",
         x = NULL, y = "Median absolute correlation") +
    theme_paper() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
  save_fig(p, "fig33_overlap_axes", width = 10, height = 4.4)
}

log_stop()
