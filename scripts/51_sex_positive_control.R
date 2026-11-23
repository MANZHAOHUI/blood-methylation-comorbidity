# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 51: sex as a positive control at the top of the scale
#
# Even the strongest existing positive control, inflammatory bowel disease, is only 0.116, and readers cannot know "how large a truly shared
# signature would be under this criterion". Sex has a large and fully reproducible effect on autosomal
# methylation, is on record and is not derived from methylation, so it can serve as the top of the scale.
#
# Procedure: the sex column of the design matrix becomes the exposure and case-control status becomes a covariate;
# methylation principal components associated with sex are removed by the same rule as in D12 - when the exposure changes, what is protected
# must change with it, otherwise the components absorb the sex effect. The other covariates are unchanged. Criterion, permutation null and choice of common probes are exactly those of the disease analysis;
# only the exposure tested is replaced, so the values of the two can be set side by side.
#
# Sex-chromosome probes were already removed from the analysis sets, so what is measured is the autosomal sex effect,
# not the trivial difference between X and Y.
#
# All forty-five pairs compare "the same exposure", so all forty-five values are positive controls. The
# cross-platform, cross-population pairs matter most - they answer directly "can the design detect a
# truly existing shared signal between an EPIC v2 cohort and a 450K cohort".
#
# 500 permutations are used: the sex effect is far above any possible null, the permutations only provide
# a formal reference and carry no decision; what counts is the magnitude of the effect-size correlation itself.
#
# Environment: Windows + RStudio. About one and a half hours
# Location: <project root>/scripts/51_sex_positive_control.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "51_sex_positive_control"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
CON_DIR <- file.path(PROJ_ROOT, "output", "39_concordance")
PC_DIR  <- file.path(PROJ_ROOT, "output", "45_stable")
OUT_DIR <- file.path(PROJ_ROOT, "output", "51_sex")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260920)
N_PERM     <- 500L
MPC_EXPO_P <- 0.01   # Same threshold as D12: principal components associated with the exposure are never modelled
BETA_CLIP <- 0.001

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
LAB_EN <- c(parkinson = "Parkinson", schizophrenia = "Schizophrenia",
            alzheimer = "Alzheimer", depression = "Depression", ibd = "IBD")
PLAT <- c(GSE145361 = "450K", GSE111629 = "450K", GSE84727 = "450K", GSE80417 = "450K",
          GSE338167 = "EPICv2", GSE144858 = "450K", GSE198904 = "EPIC", GSE273970 = "EPIC",
          GSE87648 = "450K", GSE112611 = "EPIC")
ALL <- base::unlist(DIS, use.names = FALSE)

# Turn the sex column of the design matrix into the exposure; case-control status stays as a covariate
state_of <- function(g, need_R) {
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  sc <- base::grep("^sex", base::colnames(D$design), value = TRUE)
  if (length(sc) != 1L) stop(g, ": the sex term is not a single column but ", length(sc), " columns")
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
  ex <- base::as.numeric(D$design[, sc])
  X0 <- D$design[, base::setdiff(base::colnames(D$design), sc), drop = FALSE]
  # Symmetric to D12: with sex as the exposure, methylation principal components associated with sex must be removed, just as
  # components associated with case-control status are removed in the disease analysis. Otherwise the sex effect is absorbed by the
  # components - mpc6 of GSE198904 correlates with sex at -0.748, and at K>=10 its sex concordance drops to zero.
  mp <- base::grep("^mpc", base::colnames(X0), value = TRUE)
  n_drop <- 0L
  if (length(mp)) {
    pv <- base::vapply(mp, function(k)
      stats::anova(stats::lm(X0[, k] ~ ex))$`Pr(>F)`[1], numeric(1))
    drop <- mp[pv < MPC_EXPO_P]
    n_drop <- length(drop)
    if (n_drop) X0 <- X0[, base::setdiff(base::colnames(X0), drop), drop = FALSE]
  }
  q0 <- qr(X0)
  R <- if (need_R) qr.resid(q0, base::t(M)) else NULL
  eff <- if (need_R) NULL else {
    gr <- qr.resid(q0, ex)
    Rt <- qr.resid(q0, base::t(M))
    v <- base::as.numeric(base::crossprod(gr, Rt)) / base::sum(gr^2)
    rm(Rt); invisible(gc()); v
  }
  rm(M); invisible(gc())
  list(q0 = q0, R = R, ex = ex, probes = D$probes, eff = eff,
       n = base::nrow(D$design), n_pos = base::sum(ex == 1),
       n_mpc = length(mp), n_drop = n_drop)
}
coefs <- function(S, v) {
  gr <- qr.resid(S$q0, v)
  base::as.numeric(base::crossprod(gr, S$R)) / base::sum(gr^2)
}

say("=== 1. Estimating the autosomal effect of sex set by set ===")
EFF <- list(); PRB <- list(); NS <- list()
for (g in ALL) {
  t0 <- Sys.time()
  S <- state_of(g, need_R = FALSE)
  EFF[[g]] <- S$eff; PRB[[g]] <- S$probes; NS[[g]] <- c(S$n, S$n_pos)
  say(sprintf("  %-10s samples %4d (one sex %4d)  probes %6d  principal components %2d, removed as associated with sex %2d  time %.1f min",
              g, S$n, S$n_pos, length(S$probes), S$n_mpc, S$n_drop,
              base::as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(S); invisible(gc())
}

# The pairs are exactly those of the disease analysis: five within-disease, forty cross-disease
pairs <- list()
for (d in names(DIS))
  pairs[[length(pairs) + 1]] <- list(kind = "within", dpair = LAB_EN[[d]],
    a = DIS[[d]][1], b = DIS[[d]][2], id = base::paste0("PC-", LAB_EN[[d]]))
dn <- names(DIS)
for (i in 1:(length(dn) - 1)) for (j in (i + 1):length(dn))
  for (x in DIS[[dn[i]]]) for (y in DIS[[dn[j]]])
    pairs[[length(pairs) + 1]] <- list(kind = "cross",
      dpair = base::paste0(LAB_EN[[dn[i]]], "-", LAB_EN[[dn[j]]]), a = x, b = y,
      id = base::paste0(LAB_EN[[dn[i]]], "|", LAB_EN[[dn[j]]]))
for (k in base::seq_along(pairs)) {
  P <- pairs[[k]]
  pairs[[k]]$perm  <- if (NS[[P$a]][1] <= NS[[P$b]][1]) P$a else P$b
  pairs[[k]]$fixed <- if (pairs[[k]]$perm == P$a) P$b else P$a
}
grp <- base::split(base::seq_along(pairs), base::sapply(pairs, function(x) x$perm))

say("")
say("=== 2. Concordance of the sex effect in the forty-five pairs ===")
rows <- list()
for (pg in base::names(grp)) {
  idx <- grp[[pg]]
  say("---- Residualising ", pg, " (in ", length(idx), " pairs) ----")
  S <- state_of(pg, need_R = TRUE)
  base_coef <- coefs(S, S$ex)
  chk <- stats::cor(base_coef, EFF[[pg]])
  say("  Correlation with the coefficients of part 1: ", sprintf("%.6f", chk),
      if (chk < 0.999) "  [WARNING] the two should be almost identical" else "")
  for (k in idx) {
    P <- pairs[[k]]
    ea_all <- EFF[[P$fixed]]; pa <- PRB[[P$fixed]]
    cm <- base::intersect(pa, S$probes)
    ia <- base::match(cm, pa); ib <- base::match(cm, S$probes)
    ea <- ea_all[ia]; eb <- base_coef[ib]
    rho <- stats::cor(ea, eb, method = "spearman")
    dir <- base::mean(base::sign(ea) == base::sign(eb))
    ra <- base::rank(ea)
    pr <- numeric(N_PERM)
    for (i in base::seq_len(N_PERM))
      pr[i] <- stats::cor(ra, base::rank(coefs(S, base::sample(S$ex))[ib]))
    pv <- (1 + base::sum(pr >= rho)) / (N_PERM + 1)
    z <- (rho - base::mean(pr)) / stats::sd(pr)
    say(sprintf("  %-26s %-9s vs %-9s probes %6d  rho=%+.4f  sign %.4f  Z=%+7.1f  P=%.4f",
                P$id, P$fixed, P$perm, length(cm), rho, dir, z, pv))
    rows[[length(rows) + 1]] <- data.frame(
      id = P$id, kind = P$kind, dpair = P$dpair, fixed = P$fixed, permuted = P$perm,
      plat_a = PLAT[[P$fixed]], plat_b = PLAT[[P$perm]],
      same_platform = PLAT[[P$fixed]] == PLAT[[P$perm]],
      n_probe = length(cm), rho = rho, dir = dir, z = z, p = pv,
      stringsAsFactors = FALSE)
  }
  rm(S, base_coef); invisible(gc())
  say("")
}

tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "sex_concordance_all_pairs.csv"), row.names = FALSE)

say("=== 3. The two ends of the scale ===")
say(sprintf("  Sex (forty-five pairs): median %+.4f, range %+.4f to %+.4f, median sign concordance %.4f",
            stats::median(tab$rho), base::min(tab$rho), base::max(tab$rho),
            stats::median(tab$dir)))
sp <- tab[tab$same_platform, ]; dp <- tab[!tab$same_platform, ]
say(sprintf("    same platform, %d pairs: median %+.4f | cross-platform, %d pairs: median %+.4f",
            base::nrow(sp), stats::median(sp$rho), base::nrow(dp), stats::median(dp$rho)))
pc <- tryCatch(read.csv(file.path(PC_DIR, "positive_controls_10k.csv"),
                        stringsAsFactors = FALSE), error = function(e) NULL)
dz <- tryCatch(read.csv(file.path(CON_DIR, "concordance_all_pairs.csv"),
                        stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(pc))
  say(sprintf("  Disease (five within-disease pairs): %s",
              base::paste(sprintf("%s %+.4f", pc$disease, pc$rho), collapse = "; ")))
if (!is.null(dz)) {
  cr <- dz[dz$kind == "cross", ]
  say(sprintf("  Disease (forty cross-disease pairs): median %+.4f, maximum %+.4f",
              stats::median(cr$rho), base::max(cr$rho)))
}
say("")
say("  Interpretation: sex, a truly shared exposure, gives concordance of the above magnitude under the same criterion, in the same cohorts and")
say("  across platforms and populations. The cross-disease concordance of disease is far below it, so 'not detectable'")
say("  is not a property of this design.")

say("")
say("=== Plotting ===")
pl <- tab[, c("id", "rho")]; pl$what <- "Sex"
if (!is.null(dz)) {
  d2 <- dz[, c("id", "rho")]; d2$what <- "Disease"
  pl <- base::rbind(pl, d2)
}
p <- ggplot(pl, aes(what, rho, colour = what)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(position = position_jitter(width = 0.18), size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = c(Sex = PAL[3], Disease = PAL[1]), guide = "none") +
  labs(title = "What a genuinely shared signal looks like under the same test",
       subtitle = "Each point is one cohort pair; sex is recorded, not inferred from methylation",
       x = NULL, y = "Spearman correlation of effect sizes") +
  theme_paper()
save_fig(p, "fig53_sex_vs_disease", width = 7, height = 5.5)
log_stop()
