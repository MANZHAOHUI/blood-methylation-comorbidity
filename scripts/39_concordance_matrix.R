# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 39: concordance matrix of effect sizes (primary criterion of the project)
#
# Design: five disorders with two independent cohorts each, hence
#   within-disease positive controls  5 pairs  - two independent cohorts should agree; calibrates the power of the criterion
#   cross-disease comparisons        40 pairs  - 10 disease pairs x 4 cohort pairings each
#
# Each disease pair has four pairings, the most important point of the design: a claim of cross-disease sharing
# must hold in all four pairings, rather than reporting one favourable pair. This is the weakest link of the existing comorbidity
# literature - with a single pairing it cannot be judged whether a result is real.
#
# The criterion is the one fixed in the first round:
#   statistic  Spearman correlation of effect sizes over the common probes, and sign concordance
#   null       effect sizes recomputed after randomly permuting case-control labels within a cohort. Permutation retains all
#           technical and biological structure of the cohort and only breaks the link with the disease label, so residual
#           inflation unrelated to group enters the null and cancels out
#   algorithm  Frisch-Waugh-Lovell: M and the group vector are each residualised on the remaining covariates,
#           after which the group coefficient of each probe reduces to an inner product of two vectors
#           and each permutation is a single matrix multiplication. The coefficients are identical to those of limma,
#           which the script checks itself
#
# Computation: in each pair the smaller cohort is permuted, and pairs are processed grouped by the permuted cohort,
# so that the residual matrix of each data set is computed once and only one set is in memory at a time.
#
# This step produces the concordance matrix only, without direction, stratification or overlap analyses - those come in later steps;
# the order was learned through three iterations in the first round: direction and replication first, mechanism afterwards.
#
# Environment: Windows + RStudio. Can run unattended
# Location: <project root>/scripts/39_concordance_matrix.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "39_concordance_matrix"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "39_concordance")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260917)
N_PERM    <- 1000L
BETA_CLIP <- 0.001

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
LAB <- c(parkinson = "Parkinson's disease", schizophrenia = "schizophrenia",
         alzheimer = "Alzheimer's disease", depression = "depression", ibd = "inflammatory bowel disease")
LAB_EN <- c(parkinson = "Parkinson", schizophrenia = "Schizophrenia",
            alzheimer = "Alzheimer", depression = "Depression", ibd = "IBD")

ew <- list(); dz <- list(); nsmp <- integer(0)
for (d in names(DIS)) for (g in DIS[[d]]) {
  f1 <- file.path(RES_DIR, base::paste0(g, "_ewas_final.rds"))
  f2 <- file.path(RES_DIR, base::paste0(g, "_design.rds"))
  if (!file.exists(f1) || !file.exists(f2)) stop("Results or design matrix missing for ", g, "")
  ew[[g]] <- base::readRDS(f1); dz[[g]] <- base::readRDS(f2)
  nsmp[g] <- base::nrow(dz[[g]]$design)
}
say("Results and design matrices of the ten cohorts loaded; sample sizes ",
    base::paste(sprintf("%s=%d", base::names(nsmp), nsmp), collapse = ", "))

# --------------------------------------------------------------
# List of pairs
# --------------------------------------------------------------
pairs <- list()
for (d in names(DIS))
  pairs[[length(pairs) + 1]] <- list(
    kind = "within", dpair = d, a = DIS[[d]][1], b = DIS[[d]][2],
    id = base::paste0("PC-", LAB_EN[[d]]))
dn <- names(DIS)
for (i in 1:(length(dn) - 1)) for (j in (i + 1):length(dn)) {
  for (x in DIS[[dn[i]]]) for (y in DIS[[dn[j]]]) {
    pairs[[length(pairs) + 1]] <- list(
      kind = "cross", dpair = base::paste0(dn[i], "-", dn[j]), a = x, b = y,
      id = base::paste0(LAB_EN[[dn[i]]], "|", LAB_EN[[dn[j]]]))
  }
}
# In each pair the smaller cohort is permuted
for (k in base::seq_along(pairs)) {
  P <- pairs[[k]]
  pairs[[k]]$perm <- if (nsmp[[P$a]] <= nsmp[[P$b]]) P$a else P$b
  pairs[[k]]$fixed <- if (pairs[[k]]$perm == P$a) P$b else P$a
}
say("Pairs: within-disease ", base::sum(base::sapply(pairs, function(x) x$kind == "within")),
    ", cross-disease ", base::sum(base::sapply(pairs, function(x) x$kind == "cross")),
    ", total ", length(pairs))
grp <- base::split(base::seq_along(pairs), base::sapply(pairs, function(x) x$perm))
say("Grouped by permuted cohort: ", base::paste(sprintf("%s(%d pairs)", base::names(grp),
                                              base::sapply(grp, length)),
                                      collapse = ", "))
say("")

# --------------------------------------------------------------
make_state <- function(g) {
  D <- dz[[g]]
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
  des <- D$design
  jg <- base::which(base::colnames(des) == "groupcase")
  X0 <- des[, -jg, drop = FALSE]
  q0 <- qr(X0)
  R <- qr.resid(q0, base::t(M))
  rm(M); invisible(gc())
  list(R = R, q0 = q0, probes = D$probes, group = D$group)
}
coefs <- function(S, gvec) {
  gres <- qr.resid(S$q0, gvec)
  base::as.numeric(base::crossprod(gres, S$R)) / base::sum(gres^2)
}

rows <- list(); nulls <- list()
for (pg in base::names(grp)) {
  idx <- grp[[pg]]
  say("==== Residualising ", pg, " (in ", length(idx), " pairs) ====")
  t0 <- Sys.time()
  S <- make_state(pg)
  say("  Time: ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))),
      " min; residual matrix ", base::nrow(S$R), " x ", base::ncol(S$R))
  gnum <- base::as.numeric(S$group == "case")
  base_coef <- coefs(S, gnum)
  chk <- stats::cor(base_coef, ew[[pg]]$effect[base::match(S$probes, ew[[pg]]$probe)],
                    use = "complete.obs")
  say("  Correlation between coefficients recomputed by residualisation and limma results = ", sprintf("%.6f", chk),
      if (chk < 0.999) "  [WARNING] the two should be almost identical" else "")

  for (k in idx) {
    P <- pairs[[k]]
    A <- ew[[P$fixed]]
    cm <- base::intersect(A$probe, S$probes)
    ia <- base::match(cm, A$probe); ib <- base::match(cm, S$probes)
    ea <- A$effect[ia]; eb <- base_coef[ib]
    rho <- stats::cor(ea, eb, method = "spearman")
    dir <- base::mean(base::sign(ea) == base::sign(eb))
    ra <- base::rank(ea)
    pr <- numeric(N_PERM); pd <- numeric(N_PERM)
    for (i in base::seq_len(N_PERM)) {
      cf <- coefs(S, base::sample(gnum))[ib]
      pr[i] <- stats::cor(ra, base::rank(cf))
      pd[i] <- base::mean(base::sign(ea) == base::sign(cf))
    }
    p_rho <- (1 + base::sum(pr >= rho)) / (N_PERM + 1)
    p_dir <- (1 + base::sum(pd >= dir)) / (N_PERM + 1)
    z <- (rho - base::mean(pr)) / stats::sd(pr)
    say(sprintf("  %-26s %-9s vs %-9s probes %6d  rho=%+.4f  sign %.4f  Z=%+5.2f  P=%.4f",
                P$id, P$fixed, P$perm, length(cm), rho, dir, z, p_rho))
    rows[[length(rows) + 1]] <- data.frame(
      id = P$id, kind = P$kind, dpair = P$dpair, fixed = P$fixed, permuted = P$perm,
      n_probe = length(cm), rho = rho, dir = dir, null_mean = base::mean(pr),
      null_sd = stats::sd(pr), z = z, p_rho = p_rho, p_dir = p_dir,
      stringsAsFactors = FALSE)
    nulls[[length(nulls) + 1]] <- data.frame(
      id = P$id, kind = P$kind, dpair = P$dpair, pairkey = base::paste0(P$fixed, "_", P$perm),
      rho = pr, obs = rho, stringsAsFactors = FALSE)
  }
  rm(S, base_coef); invisible(gc())
  say("")
}

tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "concordance_all_pairs.csv"), row.names = FALSE)
write.csv(base::do.call(rbind, nulls), file.path(OUT_DIR, "permutation_nulls.csv"),
          row.names = FALSE)

say("=== Within-disease positive controls ===")
w <- tab[tab$kind == "within", ]
for (i in base::seq_len(base::nrow(w)))
  say(sprintf("  %-26s rho=%+.4f  sign %.4f  Z=%+5.2f  P=%.4f",
              w$id[i], w$rho[i], w$dir[i], w$z[i], w$p_rho[i]))
say("  Reaching P < 0.05: ", base::sum(w$p_rho < 0.05), " / ", base::nrow(w),
    ". This row decides how all cross-disease results below are read: for a disorder without significant within-disease replication, ")
say("  a negative cross-disease result cannot be read as 'not shared', only as 'not detectable'.")
say("")

say("=== Cross-disease comparisons: the four pairings summarised by disease pair ===")
cr <- tab[tab$kind == "cross", ]
sumr <- list()
for (dp in base::unique(cr$dpair)) {
  s <- cr[cr$dpair == dp, ]
  nsig <- base::sum(s$p_rho < 0.05)
  same <- base::all(base::sign(s$rho) == base::sign(s$rho[1]))
  say(sprintf("  %-28s rho of the four pairings: %s; significant %d/4; direction %s",
              dp, base::paste(sprintf("%+.3f", s$rho), collapse = " "),
              nsig, if (same) "consistent" else "inconsistent"))
  sumr[[length(sumr) + 1]] <- data.frame(
    dpair = dp, n_sig = nsig, sign_consistent = same,
    rho_min = base::min(s$rho), rho_median = stats::median(s$rho),
    rho_max = base::max(s$rho), stringsAsFactors = FALSE)
}
ss <- base::do.call(rbind, sumr)
write.csv(ss, file.path(OUT_DIR, "cross_disease_summary.csv"), row.names = FALSE)
say("")
say("  Disease pairs with all four pairings significant and consistent in direction: ",
    base::sum(ss$n_sig == 4 & ss$sign_consistent), " / ", base::nrow(ss))
say("  Interpretation: one or two significant pairings out of four are chance at the pairing level and do not support a claim of sharing.")

say("")
say("=== Plotting ===")
tab$panelkind <- factor(ifelse(tab$kind == "within", "Within-disease (positive control)",
                               "Cross-disease"),
                        levels = c("Within-disease (positive control)", "Cross-disease"))
p <- ggplot(tab, aes(stats::reorder(id, rho), rho, fill = panelkind)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Within-disease (positive control)" = PAL[3],
                               "Cross-disease" = PAL[1]), name = NULL) +
  labs(title = "Effect-size concordance across all cohort pairs",
       subtitle = "Each bar is one pair of cohorts; within-disease pairs calibrate what a real shared signal looks like",
       x = NULL, y = "Spearman correlation of effect sizes") +
  theme_paper() +
  theme(axis.text.y = element_text(size = 6))
save_fig(p, "fig41_concordance_all_pairs", width = 9, height = 9)

crs <- merge(cr, ss[, c("dpair", "n_sig")], by = "dpair")
p <- ggplot(crs, aes(dpair, rho)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = p_rho < 0.05), size = 2, alpha = 0.85) +
  scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = PAL[2]),
                      labels = c("ns", "P < 0.05"), name = NULL) +
  coord_flip() +
  labs(title = "Each disease pair tested in four independent cohort pairings",
       subtitle = "A claim of sharing should hold in all four, not in one",
       x = NULL, y = "Spearman correlation of effect sizes") +
  theme_paper()
save_fig(p, "fig42_cross_disease_fourway", width = 9, height = 5.5)
log_stop()
