# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 45: make the criterion of the five within-disease positive controls stable
#
# 44 exposed something more basic than quality confounding: the criterion itself is not stable enough.
# Results of 39 and 44 under the same model, side by side:
#   Parkinson's disease  rho is +0.0430 in both, while Z is 1.81 and 1.71 and P is 0.0310 and 0.0519
# Same model, same data: significance flipped merely because the random seed of the permutations differed.
# With 1000 permutations the Monte Carlo standard error at P of about 0.05 is 0.0069, and the gap between 0.031 and 0.052
# lies entirely within random fluctuation.
#
# This also shows that the threshold "rebuild if Z changes by 0.5" in 44 was badly chosen - it was set without knowing the Monte Carlo
# error, and Z itself fluctuates by that order.
#
# As for the quality covariate, three pieces of evidence together show that no rebuild is needed:
#   (1) the effect sizes of the ten sets correlate at 0.973 to 1.000 between the two models, almost unchanged;
#   (2) in the two sets truly related to quality, quality is already absorbed by the first methylation principal component
#       (correlation 0.941 in GSE87648, 0.979 in GSE144858);
#   (3) for the IBD pair Z falls from 4.71 to 2.98 not because the effect shrinks (rho falls by only 11.5%)
#       but because the SD of the null grows by 40% - near-collinearity destabilises the estimate; it is not an adjustment for confounding.
#
# This step does one thing: raise the number of permutations from 1000 to 10000, re-run the five within-disease positive controls
# and give the Monte Carlo standard error of P, so that "above or below the line" no longer depends on the random seed.
# The original model is used, without the quality covariate.
#
# One more item is reported: the lower limit of the effect size that can be excluded for each disorder. Without it readers cannot tell
# "not shared" from "not detectable" - and two of the five disorders are known to belong to the latter.
#
# Environment: Windows + RStudio. About three hours; can run overnight
# Location: <project root>/scripts/45_positive_controls_stable.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "45_positive_controls_stable"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "45_stable")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260920)
N_PERM    <- 10000L
BETA_CLIP <- 0.001
POWER_TARGET <- 0.80

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
LAB_EN <- c(parkinson = "Parkinson", schizophrenia = "Schizophrenia",
            alzheimer = "Alzheimer", depression = "Depression", ibd = "IBD")

state_of <- function(g) {
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
  jg <- base::which(base::colnames(D$design) == "groupcase")
  X0 <- D$design[, -jg, drop = FALSE]
  q0 <- qr(X0)
  R <- qr.resid(q0, base::t(M))
  rm(M); invisible(gc())
  list(q0 = q0, R = R, probes = D$probes, group = D$group)
}
coefs <- function(S, gvec) {
  gr <- qr.resid(S$q0, gvec)
  base::as.numeric(base::crossprod(gr, S$R)) / base::sum(gr^2)
}

rows <- list(); nulls <- list()
for (d in names(DIS)) {
  a <- DIS[[d]][1]; bcoh <- DIS[[d]][2]
  say("==== ", LAB_EN[[d]], ": fixed ", a, ", permuted ", bcoh, " ====")
  EA <- base::readRDS(file.path(RES_DIR, base::paste0(a, "_ewas_final.rds")))
  t0 <- Sys.time()
  S <- state_of(bcoh)
  say("  Residualisation time ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))),
      " min; matrix ", base::nrow(S$R), " x ", base::ncol(S$R))
  gnum <- base::as.numeric(S$group == "case")
  eb_all <- coefs(S, gnum)

  cm <- base::intersect(EA$probe, S$probes)
  ia <- base::match(cm, EA$probe); ib <- base::match(cm, S$probes)
  ea <- EA$effect[ia]; eb <- eb_all[ib]
  rho <- stats::cor(ea, eb, method = "spearman")
  dir <- base::mean(base::sign(ea) == base::sign(eb))
  ra <- base::rank(ea)
  pr <- numeric(N_PERM)
  t0 <- Sys.time()
  for (i in base::seq_len(N_PERM)) {
    pr[i] <- stats::cor(ra, base::rank(coefs(S, base::sample(gnum))[ib]))
    if (i %% 2000L == 0L)
      say("    permutation ", i, " / ", N_PERM, ", elapsed ",
          sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  }
  pv <- (1 + base::sum(pr >= rho)) / (N_PERM + 1)
  se <- base::sqrt(pv * (1 - pv) / N_PERM)
  nm <- base::mean(pr); ns <- stats::sd(pr)
  z <- (rho - nm) / ns
  say(sprintf("  probes %6d  rho=%+.4f  sign %.4f", length(cm), rho, dir))
  say(sprintf("  null mean %+.5f  SD %.5f  -> Z=%+5.2f", nm, ns, z))
  say(sprintf("  permutation P = %.4f, Monte Carlo standard error %.4f (95%% interval %.4f to %.4f)",
              pv, se, base::max(0, pv - 1.96 * se), pv + 1.96 * se))
  verdict <- if (pv + 1.96 * se < 0.05) "pass" else
             if (pv - 1.96 * se > 0.05) "fail" else "borderline (interval spans 0.05)"
  say("  Decision: ", verdict)
  # Detection limit: the SD of the null decides how large an effect this pair can exclude
  thr <- nm + stats::qnorm(0.95) * ns
  need <- nm + (stats::qnorm(0.95) + stats::qnorm(POWER_TARGET)) * ns
  say(sprintf("  Detection limit: one-sided 5%% critical rho = %.4f; %.0f%% power requires rho >= %.4f",
              thr, 100 * POWER_TARGET, need))
  rows[[length(rows) + 1]] <- data.frame(
    disease = LAB_EN[[d]], fixed = a, permuted = bcoh, n_probe = length(cm),
    rho = rho, dir = dir, null_mean = nm, null_sd = ns, z = z, p = pv, p_se = se,
    verdict = verdict, rho_crit = thr, rho_for_power = need,
    stringsAsFactors = FALSE)
  nulls[[length(nulls) + 1]] <- data.frame(disease = LAB_EN[[d]], rho = pr,
                                           obs = rho, stringsAsFactors = FALSE)
  rm(S); invisible(gc())
  say("")
}

tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "positive_controls_10k.csv"), row.names = FALSE)
write.csv(base::do.call(rbind, nulls), file.path(OUT_DIR, "positive_control_nulls.csv"),
          row.names = FALSE)
say("=== Five within-disease positive controls (10000 permutations) ===")
for (i in base::seq_len(base::nrow(tab)))
  say(sprintf("  %-14s rho=%+.4f  Z=%+5.2f  P=%.4f±%.4f  %s",
              tab$disease[i], tab$rho[i], tab$z[i], tab$p[i], tab$p_se[i], tab$verdict[i]))
say("")
say("  Passed: ", base::sum(tab$verdict == "pass"), ", borderline ",
    base::sum(base::grepl("^borderline", tab$verdict)), ", fail ",
    base::sum(tab$verdict == "fail"), " / 5")
say("")
say("=== Effect sizes that can be excluded for each disorder ===")
for (i in base::seq_len(base::nrow(tab)))
  say(sprintf("  %-14s critical rho %.4f; 80%% power requires rho >= %.4f; observed %+.4f",
              tab$disease[i], tab$rho_crit[i], tab$rho_for_power[i], tab$rho[i]))
say("  Note: for disorders that fail, one can only say 'no concordance above this limit is detectable in this design', ")
say("        not 'there is no concordance'. The two statements must be kept apart in the paper.")

say("")
say("=== Plotting ===")
nd <- base::do.call(rbind, nulls)
obs <- base::unique(nd[, c("disease", "obs")])
p <- ggplot(nd, aes(rho)) +
  geom_histogram(bins = 60, fill = "grey75", colour = NA) +
  geom_vline(data = obs, aes(xintercept = obs), colour = PAL[2], linewidth = 0.6) +
  facet_wrap(~ disease, scales = "free", ncol = 3) +
  labs(title = "Within-disease replication against its permutation null",
       subtitle = "Ten thousand label permutations per disease; the line marks the observed value",
       x = "Spearman correlation of effect sizes", y = "Permutations") +
  theme_paper()
save_fig(p, "fig48_positive_controls_10k", width = 10, height = 6)
log_stop()
