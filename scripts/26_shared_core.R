# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 26: characterise the shared sites unrelated to cell composition and validate them in independent cohorts
#
# The result of 25 contradicts the hypothesis, which must be stated first:
#
#   The enrichment of PD x IBD is strongest in Q1 - the stratum least correlated with the neutrophil fraction
#   (median |r| = 0.018, fold 2.33, P = 0.002) - and not significant in Q5, the stratum most related to
#   composition (1.58, P = 0.142).
#
#   What is truly driven by composition is the negative control IBD x SCZ: the fold rises monotonically from 1.34 in Q1
#   to 3.16 in Q5, with P from 0.044 to 0.002.
#
# The hypothesis "cross-disease sharing = composition residual" therefore does not hold for the comorbid pair but holds for
# the negative control. The stratification also validates the method itself - it can separate the two kinds of sharing.
#
# It follows that PD and IBD share a set of sites unrelated to leukocyte composition. This step
# does three things, and failure of any one overturns that statement:
#
#   1   extract these sites and see whether directions agree and which genes they fall in
#   2   test whether they still agree in direction when both sides are replaced by another independent cohort
#       (GSE111629 and GSE112611) - this is the decisive gate
#   3   re-stratify by the age axis, the smoking axis and the variance of the sites themselves, to confirm that the enrichment
#       in the low strata is not created by a single stratification variable
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/26_shared_core.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "26_shared_core"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
REF_DIR <- file.path(DATA_ROOT, "probe_lists")
STR_DIR <- file.path(PROJ_ROOT, "output", "25_strata")
OUT_DIR <- file.path(PROJ_ROOT, "output", "26_core")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 2000L
N_STRATA <- 5L
TOP_FRAC <- 0.02
BETA_CLIP <- 0.001
AXIS_COHORT <- "GSE152027"
DISC <- c(a = "GSE145361", b = "GSE87648")      # The pair used for discovery
REPL <- c(a = "GSE111629", b = "GSE112611")     # The pair with both sides replaced

ew <- list(); dz <- list()
for (g in c("GSE145361", "GSE87648", "GSE111629", "GSE112611", "GSE152027")) {
  f1 <- file.path(RES_DIR, paste0(g, "_ewas_final.rds"))
  f2 <- file.path(RES_DIR, paste0(g, "_design.rds"))
  if (file.exists(f1)) ew[[g]] <- readRDS(f1)
  if (file.exists(f2)) dz[[g]] <- readRDS(f2)
}
axis_r <- read.csv(file.path(STR_DIR, "probe_composition_axis.csv"),
                   stringsAsFactors = FALSE)
anno <- if (file.exists(file.path(REF_DIR, "anno450k.rds")))
  readRDS(file.path(REF_DIR, "anno450k.rds")) else NULL
say("Loaded ", length(ew), " EWAS result sets; composition axis ", nrow(axis_r), " sites")
say("")

# --------------------------------------------------------------
say("=== 1. Extract the shared sites unrelated to composition ===")
A <- ew[[DISC[["a"]]]]; B <- ew[[DISC[["b"]]]]
common <- intersect(intersect(A$probe, B$probe), axis_r$probe)
ia <- match(common, A$probe); ib <- match(common, B$probe)
pa <- A$p[ia]; pb <- B$p[ib]; ea <- A$effect[ia]; eb <- B$effect[ib]
ar <- axis_r$absr[match(common, axis_r$probe)]
qs <- stats::quantile(ar, seq(0, 1, length.out = N_STRATA + 1))
br <- cut(ar, breaks = qs, include.lowest = TRUE, labels = paste0("Q", 1:N_STRATA))
idx_by <- split(seq_along(common), br)

core <- integer(0)
for (k in 1:2) {                                   # Q1 and Q2
  ii <- idx_by[[k]]
  n <- max(20L, floor(length(ii) * TOP_FRAC))
  core <- c(core, intersect(ii[order(pa[ii])][seq_len(n)],
                            ii[order(pb[ii])][seq_len(n)]))
}
core <- sort(unique(core))
say("  Shared sites in Q1 and Q2 together: ", length(core), "")
dirc <- mean(sign(ea[core]) == sign(eb[core]))
say("  Sign concordance: ", sprintf("%.3f", dirc),
    " (", sum(sign(ea[core]) == sign(eb[core])), " / ", length(core), ")")
say("  Median |r| with the neutrophil fraction: shared sites ", sprintf("%.3f", stats::median(ar[core])),
    ", all ", sprintf("%.3f", stats::median(ar)))
gl <- character(0)
if (!is.null(anno)) {
  g <- anno[intersect(common[core], rownames(anno)), "UCSC_RefGene_Name"]
  gl <- unique(unlist(strsplit(as.character(g), ";", fixed = TRUE)))
  gl <- sort(gl[nzchar(gl)])
  say("  Genes involved: ", length(gl), "")
  say("  First 40: ", paste(utils::head(gl, 40), collapse = ", "))
  writeLines(gl, file.path(OUT_DIR, "core_shared_genes.txt"), useBytes = TRUE)
}
# Genomic clustering: if the 75 sites come from only a few regions, the strength of evidence must be discounted
if (!is.null(anno)) {
  aa <- anno[intersect(common[core], rownames(anno)), c("chr", "pos")]
  if (nrow(aa)) {
    o <- order(aa$chr, aa$pos)
    aa <- aa[o, , drop = FALSE]
    newblock <- c(TRUE, aa$chr[-1] != aa$chr[-nrow(aa)] |
                    diff(aa$pos) > 1000)
    n1k <- sum(newblock)
    newblock10 <- c(TRUE, aa$chr[-1] != aa$chr[-nrow(aa)] |
                      diff(aa$pos) > 10000)
    n10k <- sum(newblock10)
    say("  Genomic clustering: ", nrow(aa), " sites fall into ", n1k,
        " blocks more than 1 kb apart and ", n10k, " blocks more than 10 kb apart")
    say("    mean number of sites per block ", sprintf("%.2f", nrow(aa) / n1k),
        "; a value close to 1 means the sites are independent and the evidence does not rest on a few regions")
    tb <- sort(table(aa$chr), decreasing = TRUE)
    say("    top five chromosomes: ",
        paste(sprintf("%s=%d", names(tb)[1:min(5, length(tb))],
                      as.integer(tb)[1:min(5, length(tb))]), collapse = "; "))
  }
}

write.csv(data.frame(probe = common[core], effect_pd = ea[core], p_pd = pa[core],
                     effect_ibd = eb[core], p_ibd = pb[core],
                     absr_composition = ar[core], stratum = as.character(br[core]),
                     stringsAsFactors = FALSE),
          file.path(OUT_DIR, "core_shared_probes.csv"), row.names = FALSE)
say("")

# --------------------------------------------------------------
say("=== 2. Validation in independent cohorts with both sides replaced ===")
A2 <- ew[[REPL[["a"]]]]; B2 <- ew[[REPL[["b"]]]]
cm2 <- intersect(A2$probe, B2$probe)
tgt <- intersect(common[core], cm2)
say("  Shared sites that can be matched in the replication cohorts: ", length(tgt), " / ", length(core))
if (length(tgt) >= 20L) {
  i2 <- match(tgt, A2$probe); j2 <- match(tgt, B2$probe)
  e1 <- A2$effect[i2]; e2 <- B2$effect[j2]
  obs_dir <- mean(sign(e1) == sign(e2))
  obs_rho <- stats::cor(e1, e2, method = "spearman")
  say("  In the replication cohorts: sign concordance ", sprintf("%.3f", obs_dir),
      ", Spearman rho ", sprintf("%.3f", obs_rho))

  # Null: draw the same number of sites at random from the same composition strata
  ar2 <- axis_r$absr[match(tgt, axis_r$probe)]
  pool <- cm2[cm2 %in% axis_r$probe]
  arp <- axis_r$absr[match(pool, axis_r$probe)]
  bins <- cut(arp, breaks = qs, include.lowest = TRUE)
  tb <- table(cut(ar2, breaks = qs, include.lowest = TRUE))
  pd <- numeric(N_PERM); pr <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    sel <- unlist(lapply(names(tb), function(b) {
      ii <- which(bins == b)
      if (!length(ii) || tb[[b]] == 0) return(integer(0))
      sample(ii, min(tb[[b]], length(ii)))
    }))
    pp <- pool[sel]
    k1 <- match(pp, A2$probe); k2 <- match(pp, B2$probe)
    pd[i] <- mean(sign(A2$effect[k1]) == sign(B2$effect[k2]))
    pr[i] <- stats::cor(A2$effect[k1], B2$effect[k2], method = "spearman")
  }
  say("  Null (same number of sites drawn at random from the same strata): mean sign concordance ",
      sprintf("%.3f", mean(pd)), ", mean rho ", sprintf("%.3f", mean(pr)))
  say("  Empirical P: sign ", format.pval((1 + sum(pd >= obs_dir)) / (N_PERM + 1), digits = 3),
      "; rho ", format.pval((1 + sum(pr >= obs_rho)) / (N_PERM + 1), digits = 3))
  write.csv(data.frame(probe = tgt, effect_pd_rep = e1, effect_ibd_rep = e2,
                       stringsAsFactors = FALSE),
            file.path(OUT_DIR, "core_in_replication.csv"), row.names = FALSE)
} else say("  Too few matched sites; cannot be tested")
say("")

# --------------------------------------------------------------
say("=== 3. Stratify by other variables to confirm that the low-stratum enrichment is not created by the stratification variable ===")
D <- dz[[AXIS_COHORT]]
obj <- readRDS(file.path(SET_DIR, paste0(AXIS_COHORT, "_analysis_filtered.rds")))
cov <- read.csv(file.path(COV_DIR, paste0(AXIS_COHORT, "_covariates.csv")),
                stringsAsFactors = FALSE)
cov <- cov[match(D$gsm, cov$gsm), , drop = FALSE]
ctrl <- which(D$group == "control")
b <- obj$beta[D$probes, D$gsm[ctrl], drop = FALSE]; rm(obj); invisible(gc())
b[!is.finite(b)] <- NA
b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
M <- log2(b / (1 - b)); rm(b); invisible(gc())
axes <- list(age = cov$dnam_age_horvath[ctrl], smoking = cov$dnam_smoking[ctrl])
AX <- list(variance = apply(M, 1, stats::var, na.rm = TRUE))
for (k in names(axes)) {
  v <- axes[[k]]; okk <- is.finite(v)
  xs <- scale(v[okk])[, 1]
  Mm <- M[, okk, drop = FALSE]; Mm <- Mm - rowMeans(Mm, na.rm = TRUE)
  AX[[k]] <- abs(as.numeric(Mm %*% xs) / sqrt(sum(xs^2) * rowSums(Mm^2, na.rm = TRUE)))
}
rm(M); invisible(gc())
for (k in names(AX)) AX[[k]] <- data.frame(probe = D$probes, val = AX[[k]],
                                           stringsAsFactors = FALSE)

make_state <- function(g) {
  Dg <- dz[[g]]
  o <- readRDS(file.path(SET_DIR, paste0(g, "_analysis_filtered.rds")))
  bb <- o$beta[Dg$probes, Dg$gsm, drop = FALSE]; rm(o); invisible(gc())
  bb[!is.finite(bb)] <- NA
  bb[bb < BETA_CLIP] <- BETA_CLIP; bb[bb > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  MM <- log2(bb / (1 - bb)); rm(bb); invisible(gc())
  des <- Dg$design; jg <- which(colnames(des) == "groupcase")
  X0 <- des[, -jg, drop = FALSE]; q0 <- qr(X0)
  R <- qr.resid(q0, t(MM)); rss0 <- colSums(R^2); rm(MM); invisible(gc())
  list(R = R, q0 = q0, rss0 = rss0, probes = Dg$probes,
       group = Dg$group, df = nrow(des) - ncol(des))
}
pvals <- function(S, gvec) {
  gres <- qr.resid(S$q0, gvec); den <- sum(gres^2)
  bh <- as.numeric(crossprod(gres, S$R)) / den
  rss <- S$rss0 - bh^2 * den
  se <- sqrt(pmax(rss, 1e-12) / S$df / den)
  2 * stats::pt(abs(bh / se), S$df, lower.tail = FALSE)
}
S <- make_state(DISC[["b"]])
gnum <- as.numeric(S$group == "case")

rows <- list()
for (k in names(AX)) {
  v <- AX[[k]]$val[match(common, AX[[k]]$probe)]
  qq <- stats::quantile(v, seq(0, 1, length.out = N_STRATA + 1), na.rm = TRUE)
  bb2 <- cut(v, breaks = qq, include.lowest = TRUE, labels = paste0("Q", 1:N_STRATA))
  ib2 <- split(seq_along(common), bb2)
  jb <- match(common, S$probes)
  obs <- integer(N_STRATA); nsel <- integer(N_STRATA)
  for (s in 1:N_STRATA) {
    ii <- ib2[[s]]; n <- max(20L, floor(length(ii) * TOP_FRAC)); nsel[s] <- n
    obs[s] <- length(intersect(ii[order(pa[ii])][seq_len(n)],
                               ii[order(pb[ii])][seq_len(n)]))
  }
  pm <- matrix(0L, 1000L, N_STRATA)
  for (i in seq_len(1000L)) {
    pp <- pvals(S, sample(gnum))[jb]
    for (s in 1:N_STRATA) {
      ii <- ib2[[s]]
      pm[i, s] <- length(intersect(ii[order(pa[ii])][seq_len(nsel[s])],
                                   ii[order(pp[ii])][seq_len(nsel[s])]))
    }
  }
  say("  Stratified by ", k, ": ")
  for (s in 1:N_STRATA) {
    mu <- mean(pm[, s]); pe <- (1 + sum(pm[, s] >= obs[s])) / 1001
    say(sprintf("     Q%d  overlap %4d, null %6.1f, fold %5.2f, P = %.4f",
                s, obs[s], mu, obs[s] / max(mu, 0.5), pe))
    rows[[length(rows) + 1]] <- data.frame(
      axis = k, stratum = paste0("Q", s), observed = obs[s], null_mean = mu,
      fold = obs[s] / max(mu, 0.5), p = pe, stringsAsFactors = FALSE)
  }
}
rb <- do.call(rbind, rows)
write.csv(rb, file.path(OUT_DIR, "alt_axis_strata.csv"), row.names = FALSE)

say("")
say("=== Plotting ===")
rb$stratum <- factor(rb$stratum, levels = paste0("Q", 1:N_STRATA))
p <- ggplot(rb, aes(stratum, fold)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_col(fill = PAL[2], width = 0.7) +
  geom_text(aes(label = sprintf("%.2f\nP=%.3f", fold, p)), vjust = -0.3, size = 2.2) +
  facet_wrap(~ axis, nrow = 1) +
  labs(title = "PD and IBD overlap after stratifying on other probe properties",
       subtitle = "If the low strata stay enriched under every stratifying variable, the shared signal is not an artefact of one axis",
       x = NULL, y = "Overlap divided by the permutation mean") +
  theme_paper()
save_fig(p, "fig35_alt_axis_strata", width = 10, height = 4.4)

log_stop()
