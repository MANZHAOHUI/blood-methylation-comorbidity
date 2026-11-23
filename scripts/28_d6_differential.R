# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 28: carry out the pre-specified decision D6 - test of the difference between primary criterion and negative control
#
# D6 requires that the PD-IBD correlation be significantly higher than the two correlations involving the negative control.
# So far each pair was only compared with its own permutation null and the difference was never tested directly;
# this test was written into the plan but not executed. It must be done now, not because it
# might change the conclusion - at the global scale PD-SCZ is +0.026 and PD-IBD is -0.011,
# so the direction is clear - but because "not doing it because the answer is known" is exactly the practice
# that this paper criticises. A pre-specified test must be executed and reported as it is.
#
# Statistic: D = rho(shared cohort, target cohort) - rho(shared cohort, negative control)
# Null: permute the case-control labels of the shared cohort and recompute both correlations before taking the difference.
#   This preserves the dependence between the two correlations. Permuting the two sides separately would overestimate the variance of the difference
#   and make the test conservative, which would be wrong.
#
# A power check is run as well: with the target replaced by the within-disease replication cohort, the difference should be significantly positive.
# If that cell is not significant either, the test itself has no power and a negative cannot be interpreted.
#
# Environment: Windows + RStudio (peak memory about 16 GB, when residualising GSE145361)
# Location: <project root>/scripts/28_d6_differential.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "28_d6_differential"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "28_d6")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 1000L
BETA_CLIP <- 0.001

# Each row: permute the shared set and compare its correlations with target and control
TESTS <- list(
  list(id = "D6-1", shared = "GSE145361", target = "GSE87648",
       control = "GSE152027", kind = "primary criterion (PD as the shared side)"),
  list(id = "D6-2", shared = "GSE87648", target = "GSE145361",
       control = "GSE152027", kind = "primary criterion (IBD as the shared side)"),
  list(id = "PWR",  shared = "GSE87648", target = "GSE112611",
       control = "GSE152027", kind = "power check (within-disease replication vs negative control)"))

ew <- list()
for (g in c("GSE145361", "GSE87648", "GSE112611", "GSE152027")) {
  f <- file.path(RES_DIR, paste0(g, "_ewas_final.rds"))
  if (file.exists(f)) ew[[g]] <- readRDS(f)
}
dz <- list()
for (g in names(ew)) {
  f <- file.path(RES_DIR, paste0(g, "_design.rds"))
  if (file.exists(f)) dz[[g]] <- readRDS(f)
}
say("Loaded ", length(ew), " result sets")
say("")

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
  R <- qr.resid(q0, t(M))
  rm(M); invisible(gc())
  list(R = R, q0 = q0, probes = D$probes, group = D$group)
}
coefs <- function(S, gvec) {
  gres <- qr.resid(S$q0, gvec)
  as.numeric(crossprod(gres, S$R)) / sum(gres^2)
}

rows <- list(); nulls <- list()
for (T in TESTS) {
  if (is.null(ew[[T$shared]]) || is.null(ew[[T$target]]) || is.null(ew[[T$control]])) next
  say("==== ", T$id, "  ", T$kind, " ====")
  say("  Shared side ", T$shared, "; target ", T$target, "; negative control ", T$control)

  cm <- Reduce(intersect, list(ew[[T$shared]]$probe, ew[[T$target]]$probe,
                               ew[[T$control]]$probe))
  say("  Probes common to the three: ", length(cm))
  eT <- ew[[T$target]]$effect[match(cm, ew[[T$target]]$probe)]
  eC <- ew[[T$control]]$effect[match(cm, ew[[T$control]]$probe)]
  eS <- ew[[T$shared]]$effect[match(cm, ew[[T$shared]]$probe)]
  rT <- rank(eT); rC <- rank(eC)

  rho_t <- stats::cor(rank(eS), rT)
  rho_c <- stats::cor(rank(eS), rC)
  Dobs <- rho_t - rho_c
  say(sprintf("  Observed: rho with target = %+.4f, rho with negative control = %+.4f, difference D = %+.4f",
              rho_t, rho_c, Dobs))

  say("  Residualising ", T$shared, " ...")
  t0 <- Sys.time()
  S <- make_state(T$shared)
  say("  Time: ", sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  jj <- match(cm, S$probes)
  gnum <- as.numeric(S$group == "case")
  chk <- stats::cor(coefs(S, gnum)[jj], eS, use = "complete.obs")
  say("  Correlation between coefficients recomputed by residualisation and limma = ", sprintf("%.6f", chk),
      if (chk < 0.999) "  [WARNING] should be close to 1" else "")

  pd <- numeric(N_PERM)
  t0 <- Sys.time()
  for (i in seq_len(N_PERM)) {
    cf <- rank(coefs(S, sample(gnum))[jj])
    pd[i] <- stats::cor(cf, rT) - stats::cor(cf, rC)
  }
  say("  ", N_PERM, " permutations, time ",
      sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  pe <- (1 + sum(pd >= Dobs)) / (N_PERM + 1)
  z <- (Dobs - mean(pd)) / stats::sd(pd)
  say(sprintf("  Null: mean %+.4f, SD %.4f; one-sided empirical P = %.4f, Z = %+.2f",
              mean(pd), stats::sd(pd), pe, z))
  say("")

  rows[[length(rows) + 1]] <- data.frame(
    test = T$id, kind = T$kind, shared = T$shared, target = T$target,
    control = T$control, n_probe = length(cm), rho_target = rho_t,
    rho_control = rho_c, D = Dobs, null_mean = mean(pd),
    null_sd = stats::sd(pd), z = z, p = pe, stringsAsFactors = FALSE)
  nulls[[T$id]] <- data.frame(test = T$id, kind = T$kind, D = pd, obs = Dobs,
                              stringsAsFactors = FALSE)
  rm(S); invisible(gc())
}

tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "d6_differential.csv"), row.names = FALSE)
say("=== Summary ===")
print(tab[, c("test", "kind", "rho_target", "rho_control", "D", "z", "p")])
say("")
say("Interpretation: D6 requires the two primary-criterion rows to be significantly positive. If the power-check row is significantly positive, ")
say("      the test has power, and a non-significant primary criterion can be read as a true negative.")

say("")
say("=== Plotting ===")
nd <- do.call(rbind, nulls)
lab <- paste0(tab$test, "  ", tab$kind)
nd$panel <- factor(lab[match(nd$test, tab$test)], levels = lab)
p <- ggplot(nd, aes(D)) +
  geom_histogram(bins = 45, fill = "grey75", colour = "white", linewidth = 0.2) +
  geom_vline(aes(xintercept = obs), colour = PAL[2], linewidth = 0.7) +
  facet_wrap(~ panel, ncol = 3, scales = "free") +
  labs(title = "Pre-registered contrast between the target pair and the negative control",
       subtitle = "Labels are permuted in the shared cohort so both correlations are recomputed together, preserving their dependence",
       x = "Difference in Spearman correlation", y = "Permutations") +
  theme_paper()
save_fig(p, "fig38_d6_differential", width = 10, height = 4.2)

log_stop()
