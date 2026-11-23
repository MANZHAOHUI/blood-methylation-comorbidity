# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 57: sensitivity analysis - the forty cross-disease comparisons under quality-adjusted models
#
# Background: step 43 found the sample-quality index (proportion of autosomal probes with beta between
# 0.3 and 0.7) associated with case-control status in three cohorts. Step 44 added the standardised index
# to the models and recomputed the five within-disease positive controls. This step does the same for the
# forty cross-disease comparisons, so that the whole concordance matrix is available under the
# quality-adjusted models.
#
# Model: for each of the ten cohorts, the design matrix of the primary model plus the standardised quality
# index, built exactly as in step 44 (missing index values, allowed for at most 2% of samples, are set to
# the cohort median; columns that become collinear are dropped). Both cohorts of every pair use the
# adjusted model.
#
# Statistic and null: as in step 39. Spearman correlation of the case-control effect vectors over the
# common probes; null by permuting case-control labels in the cohort that was permuted in the primary
# analysis (roles are read from the table of step 46); 1000 permutations per pair. The group coefficient of
# every probe is computed as M g~ / (g~' g~), where g~ is the group vector residualised on the other
# covariates. This equals the Frisch-Waugh-Lovell coefficient used in step 39 (residualising M as well does
# not change the inner product) and the limma coefficient; for the primary design the script checks this
# against the stored results and stops if they disagree.
#
# One difference in computation, none in the test: pairs that share a permuted cohort share its 1000 label
# permutations (each permuted effect vector is compared with the fixed vectors of all its pairs). Every
# pair still has its own valid permutation null of 1000 draws.
#
# The five within-disease pairs are not permuted again; their observed values under the adjusted models are
# recomputed and compared with the table of step 44, which holds their permutation results.
#
# Cohorts are processed in descending order of sample size. In every pair the smaller cohort is the
# permuted one, so the adjusted effect vector of the fixed cohort is always available when it is needed and
# every analysis set is read only once.
#
# This step changes no existing result. Output: output/57_quality_cross/
#   quality_adjusted_cross_pairs.csv     the forty pairs, primary values next to the adjusted values
#   quality_adjusted_cross_summary.csv   the four pairings summarised by disease pair
#   quality_adjusted_within_check.csv    observed within-disease values against step 44
#   quality_adjusted_model_info.csv      per cohort: parameters, quality term, agreement of effect vectors
#   quality_adjusted_nulls.csv           permutation nulls
#   quality_adjusted_key_numbers.csv     the numbers quoted in the manuscript
#   figure fig54_quality_adjusted_cross_pairs (PNG 600 dpi and PDF) in output/figures/
#
# Environment: Windows + RStudio. About two hours. Can run unattended
# Location: <project root>/scripts/57_quality_adjusted_cross_disease.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "57_quality_adjusted_cross_disease"

SET_DIR  <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR  <- file.path(DATA_ROOT, "ewas_final")
QC_FILE  <- file.path(PROJ_ROOT, "output", "43_qc", "bimodality_per_sample.csv")
Q44_FILE <- file.path(PROJ_ROOT, "output", "44_quality", "within_disease_quality_sensitivity.csv")
REF_FILE <- file.path(PROJ_ROOT, "output", "46_limits", "all_pairs_with_limits.csv")
OUT_DIR  <- file.path(PROJ_ROOT, "output", "57_quality_cross")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

N_PERM       <- 1000L     # as in step 39
BETA_CLIP    <- 0.001
POWER_TARGET <- 0.80
SEED         <- 20261124L
Z95 <- stats::qnorm(0.95); ZPW <- stats::qnorm(POWER_TARGET)

# --------------------------------------------------------------
# Functions
# --------------------------------------------------------------
read_rds <- function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  base::readRDS(f)
}
load_M <- function(g, D) {
  obj <- read_rds(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  mp <- base::setdiff(D$probes, base::rownames(obj$beta))
  ms <- base::setdiff(D$gsm, base::colnames(obj$beta))
  if (length(mp) || length(ms))
    stop(g, ": the analysis set lacks ", length(mp), " probes and ", length(ms),
         " samples of the design; the analysis set and the design do not belong together")
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  M <- base::log2(b / (1 - b)); rm(b); invisible(gc())
  if (base::anyNA(M)) stop(g, ": the M matrix contains missing values; the design was built on complete probes, so this must not happen")
  M
}
null_qr <- function(des) {
  jg <- base::which(base::colnames(des) == "groupcase")
  if (length(jg) != 1L) stop("The design matrix must contain exactly one column 'groupcase'; found ", length(jg))
  q0 <- base::qr(des[, -jg, drop = FALSE])
  if (q0$rank < base::ncol(des) - 1L) stop("The covariate part of the design matrix is rank deficient")
  q0
}
coef_group <- function(M, q0, gvec) {
  gres <- base::qr.resid(q0, gvec)
  base::as.numeric(M %*% gres) / base::sum(gres^2)
}
null_stats <- function(rho, pr) {
  nm <- base::mean(pr); ns <- stats::sd(pr)
  pv <- (1 + base::sum(pr >= rho)) / (length(pr) + 1)
  data.frame(null_mean = nm, null_sd = ns, z = (rho - nm) / ns, p = pv,
             rho_crit = nm + Z95 * ns, rho_power = nm + (Z95 + ZPW) * ns)
}
mins_since <- function(t0) sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins")))

# --------------------------------------------------------------
say("=== 1. Inputs ===")
for (f in c(QC_FILE, Q44_FILE, REF_FILE)) if (!file.exists(f)) stop("File not found: ", f)
REF <- read.csv(REF_FILE, stringsAsFactors = FALSE)
need <- c("id", "kind", "dpair", "fixed", "permuted", "n_probe", "rho", "p_rho", "rho_crit", "rho_power")
if (!base::all(need %in% base::names(REF)))
  stop("The table of step 46 lacks columns: ", base::paste(base::setdiff(need, base::names(REF)), collapse = ", "))
if (base::nrow(REF) != 45L || base::sum(REF$kind == "cross") != 40L || base::sum(REF$kind == "within") != 5L)
  stop("The table of step 46 should have 40 cross-disease and 5 within-disease pairs; found ",
       base::sum(REF$kind == "cross"), " and ", base::sum(REF$kind == "within"))
QS <- read.csv(QC_FILE, stringsAsFactors = FALSE)
if (!base::all(c("gse", "gsm", "frac_mid") %in% base::names(QS))) stop("The quality table of step 43 lacks the columns gse, gsm, frac_mid")
Q44 <- read.csv(Q44_FILE, stringsAsFactors = FALSE)
if (base::nrow(Q44) != 10L || !base::all(c("disease", "rho") %in% base::names(Q44)))
  stop("The table of step 44 should have two rows (primary model, quality-adjusted model) for each of five disorders")

COH <- base::unique(c(REF$fixed, REF$permuted))
if (length(COH) != 10L) stop("Expected ten cohorts in the table of step 46; found ", length(COH))
DZ <- list(); nsmp <- integer(0)
for (g in COH) {
  DZ[[g]] <- read_rds(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  if (base::nrow(DZ[[g]]$design) != length(DZ[[g]]$gsm) || length(DZ[[g]]$group) != length(DZ[[g]]$gsm))
    stop(g, ": design, sample list and group vector differ in length")
  nsmp[g] <- length(DZ[[g]]$gsm)
}
ORD <- base::names(nsmp)[base::order(-nsmp)]
say("Ten cohorts, in order of processing: ", base::paste(sprintf("%s=%d", ORD, nsmp[ORD]), collapse = ", "))
pos <- stats::setNames(base::seq_along(ORD), ORD)
bad <- pos[REF$fixed] >= pos[REF$permuted]
if (base::any(bad))
  stop("In ", base::sum(bad), " pairs the fixed cohort is not larger than the permuted cohort (first: ", REF$id[base::which(bad)[1]],
       "); the processing order assumed by this script does not hold")
say("Pairs: cross-disease ", base::sum(REF$kind == "cross"), " (", N_PERM, " permutations each), within-disease ",
    base::sum(REF$kind == "within"), " (observed value only, checked against step 44)")

# --------------------------------------------------------------
say("")
say("=== 2. Quality-adjusted models and permutations, cohort by cohort ===")
EFF <- list(); info <- list(); rows <- list(); nulls <- list(); wchk <- list()
t_all <- Sys.time()
for (g in ORD) {
  D <- DZ[[g]]
  say("==== ", g, " (", nsmp[[g]], " samples) ====")
  t0 <- Sys.time()
  M <- load_M(g, D)
  say("  M matrix ", base::nrow(M), " x ", base::ncol(M), ", loaded in ", mins_since(t0), " min")
  gnum <- base::as.numeric(D$group == "case")

  # Primary design: the recomputed coefficients must reproduce the stored limma results
  E <- read_rds(file.path(RES_DIR, base::paste0(g, "_ewas_final.rds")))
  e_old <- E$effect[base::match(D$probes, E$probe)]; rm(E)
  if (base::anyNA(e_old)) stop(g, ": probes of the design are missing from the stored results")
  e0 <- coef_group(M, null_qr(D$design), gnum)
  r_chk <- stats::cor(e0, e_old)
  say("  Recomputed coefficients against the stored limma results: r = ", sprintf("%.8f", r_chk),
      ", largest absolute difference ", sprintf("%.2e", base::max(base::abs(e0 - e_old))))
  if (!base::is.finite(r_chk) || r_chk < 0.999)
    stop(g, ": recomputed coefficients do not reproduce the stored results (r = ", r_chk, "); the design or the analysis set has changed, stopping")

  # Quality-adjusted design, built as in step 44
  q <- QS[QS$gse == g, , drop = FALSE]
  fq <- q$frac_mid[base::match(D$gsm, q$gsm)]
  n_miss <- base::sum(base::is.na(fq))
  if (n_miss / length(fq) > 0.02)
    stop(g, ": quality index missing for ", sprintf("%.1f%%", 100 * n_miss / length(fq)), " of the samples")
  fq[base::is.na(fq)] <- stats::median(fq, na.rm = TRUE)
  fqz <- (fq - base::mean(fq)) / stats::sd(fq)
  mpc <- base::grep("^mpc", base::colnames(D$design), value = TRUE)
  max_r <- if (length(mpc)) base::max(base::abs(base::sapply(mpc, function(k) stats::cor(D$design[, k], fqz)))) else NA_real_
  des1 <- base::cbind(D$design, quality = fqz)
  qq <- base::qr(des1)
  dropped <- character(0)
  if (qq$rank < base::ncol(des1)) {
    keepc <- base::sort(qq$pivot[base::seq_len(qq$rank)])
    dropped <- base::colnames(des1)[base::setdiff(base::seq_len(base::ncol(des1)), keepc)]
    des1 <- des1[, keepc, drop = FALSE]
  }
  if (!("groupcase" %in% base::colnames(des1))) stop(g, ": the group column became collinear after adding the quality index")
  has_q <- "quality" %in% base::colnames(des1)
  q1 <- null_qr(des1)
  e1 <- coef_group(M, q1, gnum)
  r01 <- stats::cor(e0, e1)
  say("  Samples without a quality index: ", n_miss, "; parameters ", base::ncol(D$design), " -> ", base::ncol(des1),
      if (length(dropped)) base::paste0("; columns dropped as collinear: ", base::paste(dropped, collapse = ", ")) else "",
      if (!has_q) "   [NOTE] the quality term itself is collinear with the existing covariates and was dropped" else "")
  say("  Correlation of the effect vectors of the two models: r = ", sprintf("%.4f", r01),
      if (length(mpc)) base::paste0("; largest absolute correlation of the index with the ", length(mpc),
                                    " methylation principal components ", sprintf("%.3f", max_r)) else "; no methylation principal components in this model")
  EFF[[g]] <- stats::setNames(e1, D$probes)
  info[[length(info) + 1]] <- data.frame(
    gse = g, n_samples = nsmp[[g]], n_probes = length(D$probes), n_quality_missing = n_miss,
    parameters_primary = base::ncol(D$design), parameters_adjusted = base::ncol(des1),
    quality_term_kept = has_q, columns_dropped = base::paste(dropped, collapse = ";"),
    max_abs_cor_quality_mpc = max_r, effect_cor_primary_vs_adjusted = r01,
    check_cor_with_stored = r_chk, stringsAsFactors = FALSE)

  # Pairs in which this cohort is the permuted one
  PG <- REF[REF$permuted == g, , drop = FALSE]
  if (base::nrow(PG)) {
    prep <- list()
    for (k in base::seq_len(base::nrow(PG))) {
      ef <- EFF[[PG$fixed[k]]]
      if (is.null(ef)) stop("The adjusted effect vector of ", PG$fixed[k], " is not available yet; processing order violated")
      cm <- base::intersect(base::names(ef), D$probes)
      if (length(cm) != PG$n_probe[k])
        stop(PG$id[k], " (", PG$fixed[k], " vs ", g, "): ", length(cm), " common probes, primary analysis ", PG$n_probe[k],
             "; the probe sets have changed, stopping")
      ia <- base::match(cm, base::names(ef)); ib <- base::match(cm, D$probes)
      ea <- base::unname(ef[ia]); eb <- e1[ib]
      prep[[k]] <- list(ib = ib, ra = base::rank(ea), n = length(cm),
                        rho = stats::cor(ea, eb, method = "spearman"),
                        dir = base::mean(base::sign(ea) == base::sign(eb)))
    }
    kw <- base::which(PG$kind == "within"); kc <- base::which(PG$kind == "cross")
    for (k in kw) {
      wchk[[length(wchk) + 1]] <- data.frame(
        id = PG$id[k], fixed = PG$fixed[k], permuted = g, n_probe = prep[[k]]$n,
        rho_primary = PG$rho[k], rho_adjusted = prep[[k]]$rho, stringsAsFactors = FALSE)
      say(sprintf("  %-24s %-9s vs %-9s probes %6d  primary rho=%+.4f | adjusted rho=%+.4f  (within-disease, observed value only)",
                  PG$id[k], PG$fixed[k], g, prep[[k]]$n, PG$rho[k], prep[[k]]$rho))
    }
    if (length(kc)) {
      set.seed(SEED + pos[[g]])
      pr <- base::matrix(NA_real_, N_PERM, length(kc))
      t0 <- Sys.time()
      for (i in base::seq_len(N_PERM)) {
        cf <- coef_group(M, q1, base::sample(gnum))
        for (jj in base::seq_along(kc))
          pr[i, jj] <- stats::cor(prep[[kc[jj]]]$ra, base::rank(cf[prep[[kc[jj]]]$ib]))
        if (i %% 250L == 0L) say("    permutation ", i, " / ", N_PERM, ", elapsed ", mins_since(t0), " min")
      }
      if (base::anyNA(pr)) stop("Missing values in the permutation null")
      for (jj in base::seq_along(kc)) {
        k <- kc[jj]
        st <- null_stats(prep[[k]]$rho, pr[, jj])
        say(sprintf("  %-24s %-9s vs %-9s probes %6d  primary rho=%+.4f P=%.4f | adjusted rho=%+.4f Z=%+5.2f P=%.4f",
                    PG$id[k], PG$fixed[k], g, prep[[k]]$n, PG$rho[k], PG$p_rho[k], prep[[k]]$rho, st$z, st$p))
        rows[[length(rows) + 1]] <- cbind(
          data.frame(id = PG$id[k], kind = PG$kind[k], dpair = PG$dpair[k], fixed = PG$fixed[k], permuted = g,
                     n_probe = prep[[k]]$n, n_perm = N_PERM, rho_primary = PG$rho[k], p_primary = PG$p_rho[k],
                     rho_crit_primary = PG$rho_crit[k], rho_power_primary = PG$rho_power[k],
                     rho = prep[[k]]$rho, dir = prep[[k]]$dir, stringsAsFactors = FALSE), st)
        nulls[[length(nulls) + 1]] <- data.frame(id = PG$id[k], fixed = PG$fixed[k], permuted = g,
                                                 rho = pr[, jj], obs = prep[[k]]$rho, stringsAsFactors = FALSE)
      }
      rm(pr)
    }
    rm(prep)
  } else say("  Not the permuted cohort of any pair; effect vector only")
  rm(M, D, e0, e1, e_old); invisible(gc())
  say("  Elapsed in total: ", mins_since(t_all), " min")
  say("")
}

# --------------------------------------------------------------
say("=== 3. Summary ===")
tab <- base::do.call(rbind, rows)
if (base::nrow(tab) != 40L) stop("Expected 40 cross-disease pairs; computed ", base::nrow(tab))
tab <- tab[base::match(base::paste(REF$fixed, REF$permuted)[REF$kind == "cross"], base::paste(tab$fixed, tab$permuted)), ]
if (base::anyNA(tab$rho)) stop("Computed pairs do not match the table of step 46")
tab$delta_rho <- tab$rho - tab$rho_primary
tab$sig_primary <- tab$p_primary < 0.05
tab$sig <- tab$p < 0.05
write.csv(tab, file.path(OUT_DIR, "quality_adjusted_cross_pairs.csv"), row.names = FALSE)
write.csv(base::do.call(rbind, nulls), file.path(OUT_DIR, "quality_adjusted_nulls.csv"), row.names = FALSE)
inf <- base::do.call(rbind, info)
write.csv(inf, file.path(OUT_DIR, "quality_adjusted_model_info.csv"), row.names = FALSE)

# Within-disease observed values against step 44 (second row of each disorder = quality-adjusted model)
wc <- base::do.call(rbind, wchk)
if (base::nrow(wc) != 5L) stop("Expected 5 within-disease pairs; found ", base::nrow(wc))
wc$disease <- base::sub("^PC-", "", wc$id)
wc$rho_step44 <- NA_real_; wc$p_step44 <- NA_real_
for (i in base::seq_len(base::nrow(wc))) {
  s <- Q44[Q44$disease == wc$disease[i], , drop = FALSE]
  if (base::nrow(s) != 2L) stop("Step 44 has ", base::nrow(s), " rows for ", wc$disease[i], "; expected 2")
  wc$rho_step44[i] <- s$rho[2]
  if ("p" %in% base::names(s)) wc$p_step44[i] <- s$p[2]
}
wc$agrees_with_step44 <- base::abs(wc$rho_adjusted - wc$rho_step44) < 1e-6
write.csv(wc, file.path(OUT_DIR, "quality_adjusted_within_check.csv"), row.names = FALSE)
say("  Within-disease pairs under the adjusted models (permutation results are those of step 44):")
for (i in base::seq_len(base::nrow(wc)))
  say(sprintf("    %-16s primary rho=%+.4f  adjusted rho=%+.4f  step 44 rho=%+.4f P=%.4f  %s",
              wc$disease[i], wc$rho_primary[i], wc$rho_adjusted[i], wc$rho_step44[i], wc$p_step44[i],
              if (wc$agrees_with_step44[i]) "agrees" else "[WARNING] differs from step 44"))
if (!base::all(wc$agrees_with_step44)) {
  say("  [WARNING] ", base::sum(!wc$agrees_with_step44), " within-disease value(s) differ from step 44. The quality table of step 43 or a model has changed since step 44 was run;")
  say("            re-run step 44 before the within-disease and cross-disease results are reported together.")
} else say("  All five observed values agree with step 44.")
ibd_adj <- wc$rho_adjusted[wc$disease == "IBD"]
ibd_pri <- wc$rho_primary[wc$disease == "IBD"]
if (length(ibd_adj) != 1L) stop("The within-disease pair of inflammatory bowel disease was not found")

say("")
say("  Cross-disease pairs by disease pair:")
sumr <- list()
for (dp in base::unique(tab$dpair)) {
  s <- tab[tab$dpair == dp, , drop = FALSE]
  if (base::nrow(s) != 4L) stop("Disease pair ", dp, " has ", base::nrow(s), " pairings; expected 4")
  say(sprintf("    %-28s adjusted rho %s; significant %d/4 (primary %d/4)",
              dp, base::paste(sprintf("%+.3f", s$rho), collapse = " "), base::sum(s$sig), base::sum(s$sig_primary)))
  sumr[[length(sumr) + 1]] <- data.frame(
    dpair = dp, n_sig_primary = base::sum(s$sig_primary), n_sig = base::sum(s$sig),
    sign_consistent = base::all(base::sign(s$rho) == base::sign(s$rho[1])),
    rho_min = base::min(s$rho), rho_median = stats::median(s$rho), rho_max = base::max(s$rho),
    power_min = base::min(s$rho_power), power_max = base::max(s$rho_power), stringsAsFactors = FALSE)
}
ss <- base::do.call(rbind, sumr)
write.csv(ss, file.path(OUT_DIR, "quality_adjusted_cross_summary.csv"), row.names = FALSE)

r_p <- stats::cor(tab$rho_primary, tab$rho)
r_s <- stats::cor(tab$rho_primary, tab$rho, method = "spearman")
say("")
say("  Forty cross-disease pairs, primary -> quality-adjusted:")
say("    nominally significant (P < 0.05): ", base::sum(tab$sig_primary), " -> ", base::sum(tab$sig))
say("    positive rho: ", base::sum(tab$rho_primary > 0), " -> ", base::sum(tab$rho > 0))
say(sprintf("    median rho %+.4f -> %+.4f; range %+.4f to %+.4f -> %+.4f to %+.4f",
            stats::median(tab$rho_primary), stats::median(tab$rho), base::min(tab$rho_primary), base::max(tab$rho_primary),
            base::min(tab$rho), base::max(tab$rho)))
say("    largest number of significant pairings per disease pair: ", base::max(ss$n_sig_primary), " -> ", base::max(ss$n_sig),
    "; disease pairs significant in all four pairings: ", base::sum(ss$n_sig_primary == 4), " -> ", base::sum(ss$n_sig == 4))
say(sprintf("    detection limits (%.0f%% power) %.4f to %.4f -> %.4f to %.4f",
            100 * POWER_TARGET, base::min(tab$rho_power_primary), base::max(tab$rho_power_primary),
            base::min(tab$rho_power), base::max(tab$rho_power)))
say(sprintf("    pairs with a detection limit below the within-disease rho of inflammatory bowel disease: %d (primary, rho %.4f) -> %d (adjusted, rho %.4f) of 40",
            base::sum(tab$rho_power_primary < ibd_pri), ibd_pri, base::sum(tab$rho_power < ibd_adj), ibd_adj))
say(sprintf("    agreement of the forty rho values between the two models: Pearson r = %.4f, Spearman rho = %.4f; largest absolute change %.4f",
            r_p, r_s, base::max(base::abs(tab$delta_rho))))
chg <- tab[tab$sig != tab$sig_primary, , drop = FALSE]
if (base::nrow(chg)) {
  say("    pairs whose status at P < 0.05 changed: ", base::nrow(chg))
  for (i in base::seq_len(base::nrow(chg)))
    say(sprintf("      %-24s %s vs %s: rho %+.4f -> %+.4f, P %.4f -> %.4f", chg$id[i], chg$fixed[i], chg$permuted[i],
                chg$rho_primary[i], chg$rho[i], chg$p_primary[i], chg$p[i]))
} else say("    no pair changed its status at P < 0.05")
say(sprintf("  Effect vectors of the two models, ten cohorts: r from %.4f to %.4f; quality term kept in %d of 10 models",
            base::min(inf$effect_cor_primary_vs_adjusted), base::max(inf$effect_cor_primary_vs_adjusted), base::sum(inf$quality_term_kept)))

key <- data.frame(
  item = c("n_sig_primary", "n_sig", "n_positive_primary", "n_positive",
           "rho_median_primary", "rho_median", "rho_min_primary", "rho_min", "rho_max_primary", "rho_max",
           "max_sig_per_dpair_primary", "max_sig_per_dpair", "dpairs_all_four_primary", "dpairs_all_four",
           "power_min_primary", "power_min", "power_max_primary", "power_max",
           "ibd_within_rho_primary", "ibd_within_rho_adjusted", "limits_below_ibd_primary", "limits_below_ibd",
           "rho_pearson_primary_vs_adjusted", "rho_spearman_primary_vs_adjusted", "max_abs_delta_rho",
           "n_status_changed", "effect_cor_min", "effect_cor_max", "n_quality_term_kept", "within_values_agree_with_step44"),
  value = c(base::sum(tab$sig_primary), base::sum(tab$sig), base::sum(tab$rho_primary > 0), base::sum(tab$rho > 0),
            stats::median(tab$rho_primary), stats::median(tab$rho), base::min(tab$rho_primary), base::min(tab$rho),
            base::max(tab$rho_primary), base::max(tab$rho),
            base::max(ss$n_sig_primary), base::max(ss$n_sig), base::sum(ss$n_sig_primary == 4), base::sum(ss$n_sig == 4),
            base::min(tab$rho_power_primary), base::min(tab$rho_power), base::max(tab$rho_power_primary), base::max(tab$rho_power),
            ibd_pri, ibd_adj, base::sum(tab$rho_power_primary < ibd_pri), base::sum(tab$rho_power < ibd_adj),
            r_p, r_s, base::max(base::abs(tab$delta_rho)),
            base::nrow(chg), base::min(inf$effect_cor_primary_vs_adjusted), base::max(inf$effect_cor_primary_vs_adjusted),
            base::sum(inf$quality_term_kept), base::as.numeric(base::all(wc$agrees_with_step44))),
  stringsAsFactors = FALSE)
write.csv(key, file.path(OUT_DIR, "quality_adjusted_key_numbers.csv"), row.names = FALSE)

# --------------------------------------------------------------
say("")
say("=== 4. Figure ===")
tab$status <- factor(base::ifelse(tab$sig & tab$sig_primary, "P < 0.05 in both",
                     base::ifelse(tab$sig, "P < 0.05, adjusted only",
                     base::ifelse(tab$sig_primary, "P < 0.05, primary only", "Not significant in either"))),
                     levels = c("Not significant in either", "P < 0.05 in both", "P < 0.05, primary only", "P < 0.05, adjusted only"))
lim <- base::range(c(tab$rho_primary, tab$rho)); lim <- lim + c(-1, 1) * 0.08 * base::diff(lim)
p <- ggplot(tab, aes(rho_primary, rho, colour = status)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(size = 2, alpha = 0.9) +
  scale_colour_manual(values = c("Not significant in either" = "grey60", "P < 0.05 in both" = PAL[2],
                                 "P < 0.05, primary only" = PAL[1], "P < 0.05, adjusted only" = PAL[3]),
                      drop = FALSE, name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Concordance, primary models (Spearman correlation)",
       y = "Concordance, quality-adjusted models (Spearman correlation)") +
  theme_paper() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 2))
save_fig(p, "fig54_quality_adjusted_cross_pairs", width = 5.5, height = 6)

say("")
say("Tables written to ", OUT_DIR)
log_stop()
