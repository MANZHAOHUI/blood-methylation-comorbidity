# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 43: bimodality QC of the ten cohorts
#
# Background: 42 found 15 samples in GSE87648 whose median proportion of autosomal probes between 0.3 and 0.7
# is 0.333, against 0.187 in normal samples, almost double (Wilcoxon p = 9.9e-08).
# Their misclassified sex is only a symptom; the substance is poor overall data quality.
#
# The implication goes beyond one cohort: our QC never used this index. The existing QC checked
# detection P values and missing rates, not the bimodality of beta. In GSE87648, 3.9% of samples of this kind
# passed QC, and the other nine sets probably contain some too.
#
# The direction of the effect is clear: poor-quality samples add noise and lower within-disease concordance. The main conclusion
# is precisely that "within-disease replication fails for three of five disorders" - if part of this is due to noise,
# the conclusion is overstated. It must therefore be checked, and before conclusions are drawn.
#
# Index: methylation beta values should be bimodal (near 0 or near 1) with few intermediate values. The proportion of probes between 0.3 and
# 0.7 is therefore a sample-level quality index; the higher, the more blurred the signal.
#
# The criterion is fixed beforehand and not changed after seeing results: within each cohort, samples above the median plus 3 median absolute
# deviations are flagged as suspect. A within-cohort criterion rather than a uniform threshold is used because platform,
# normalisation and submitter processing differ between cohorts and absolute levels are not comparable.
#
# Whether suspect samples are related to case-control status is tested as well - if concentrated in one group, removal would change the composition
# of the two groups, with consequences far worse than random loss.
#
# This step only diagnoses and writes lists; it changes no existing result.
#
# Environment: Windows + RStudio. The analysis sets are loaded one by one, about twenty minutes
# Location: <project root>/scripts/43_bimodality_qc.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "43_bimodality_qc"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "43_qc")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260919)
N_PROBE <- 50000L     # An estimate from a random sample is sufficient; the difference between full and sampled is far smaller than that between samples
MAD_K   <- 3          # Criterion: median plus 3 median absolute deviations

# The sample set is always the samples in the model (gsm of *_design.rds), in agreement with the EWAS, figures, tables and text.
# Additional samples of the analysis set (not modelled because of missing covariates, duplicates removed in step 54) do not take part.

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
LAB <- c(parkinson = "Parkinson's disease", schizophrenia = "schizophrenia",
         alzheimer = "Alzheimer's disease", depression = "depression", ibd = "inflammatory bowel disease")

rows <- list(); per_sample <- list()
for (d in names(DIS)) for (g in DIS[[d]]) {
  say("==== ", g, " (", LAB[[d]], ") ====")
  f <- file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds"))
  if (!file.exists(f)) { say("  Analysis set not found, skipped"); next }
  obj <- base::readRDS(f)
  B <- obj$beta; ph <- obj$pheno
  rm(obj); invisible(gc())
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  miss <- base::setdiff(D$gsm, base::colnames(B))
  if (length(miss)) stop(g, ": samples of the model not found in the analysis set: ", length(miss), ": ",
                         base::paste(utils::head(miss), collapse = ", "))
  dropc <- base::which(!base::colnames(B) %in% D$gsm)
  if (length(dropc))
    say("  Analysis set: ", base::ncol(B), " samples, of which ", length(dropc), " are not in the model and do not take part: ",
        base::paste(base::colnames(B)[dropc], collapse = ", "))
  keep <- base::match(D$gsm, base::colnames(B))
  B <- B[, keep, drop = FALSE]; ph <- ph[keep, , drop = FALSE]
  rm(D)
  ii <- base::sample(base::seq_len(base::nrow(B)), base::min(N_PROBE, base::nrow(B)))
  fr <- base::colMeans(B[ii, , drop = FALSE] > 0.3 & B[ii, , drop = FALSE] < 0.7,
                       na.rm = TRUE)
  rm(B); invisible(gc())
  md <- stats::median(fr); ma <- stats::mad(fr)
  cut <- md + MAD_K * ma
  flag <- fr > cut
  say("  Samples ", length(fr), "; proportion of intermediate probes: median ", sprintf("%.4f", md),
      ", median absolute deviation ", sprintf("%.4f", ma))
  say("  Criterion: median + ", MAD_K, " median absolute deviations = ", sprintf("%.4f", cut),
      "; above it: ", base::sum(flag), " (", sprintf("%.1f%%", 100 * base::mean(flag)), ")")
  if (base::sum(flag))
    say("  Value range of suspect samples: ", sprintf("%.4f to %.4f", base::min(fr[flag]), base::max(fr[flag])))

  # The phenotype column names differ between the two rounds and cannot be hard-coded. A failure to find them must be reported loudly -
  # last time NA was filled in silently, four cohorts skipped the test and the summary still said "unrelated in all".
  gcol <- base::grep("^(group|status|disease_status|case_control)$",
                     base::names(ph), ignore.case = TRUE, value = TRUE)[1]
  icol <- base::grep("^(gsm|geo_accession|sample|id)$",
                     base::names(ph), ignore.case = TRUE, value = TRUE)[1]
  if (base::is.na(gcol) || base::is.na(icol)) {
    say("  [FAIL] group column or ID column not found in the phenotype table; existing column names: ",
        base::paste(base::names(ph), collapse = ", "))
    grp <- base::rep(NA_character_, length(fr))
  } else {
    grp <- base::as.character(ph[[gcol]])[base::match(base::names(fr), ph[[icol]])]
    if (base::mean(base::is.na(grp)) > 0.05)
      say("  [FAIL] column '", icol, "' does not match the sample IDs; missing ",
          sprintf("%.1f%%", 100 * base::mean(base::is.na(grp))))
    else say("  Group taken from column '", gcol, "', ID taken from column '", icol, "'")
  }
  pv <- NA_real_
  tested <- FALSE
  if (base::sum(flag) > 0 && length(base::unique(grp[!base::is.na(grp)])) == 2L &&
      base::mean(base::is.na(grp)) <= 0.05) {
    tested <- TRUE
    ct <- base::table(base::ifelse(flag, "suspect", "normal"), grp)
    for (r in base::rownames(ct))
      say("    ", r, ": ", base::paste(sprintf("%s=%d", base::colnames(ct),
                                               base::as.integer(ct[r, ])), collapse = "; "))
    if (base::all(base::dim(ct) == 2)) {
      pv <- suppressWarnings(stats::fisher.test(ct)$p.value)
      say("    Fisher's exact test against group: p = ", format.pval(pv, digits = 3),
          if (!base::is.na(pv) && pv < 0.05)
            "  [NOTE] suspect samples are unevenly distributed between the groups; removal would change the composition of the two groups"
          else "  Balanced distribution; removal amounts to a random reduction in sample size")
    }
  }
  rows[[length(rows) + 1]] <- data.frame(
    gse = g, disease = LAB[[d]], n = length(fr), median_frac = md, mad_frac = ma,
    cutoff = cut, n_flag = base::sum(flag), pct_flag = 100 * base::mean(flag),
    fisher_p = pv, group_tested = tested,
    worst = base::max(fr), worst_vs_median = base::max(fr) / md,
    stringsAsFactors = FALSE)
  per_sample[[length(per_sample) + 1]] <- data.frame(
    gse = g, gsm = base::names(fr), frac_mid = base::as.numeric(fr),
    group = grp, flagged = flag, stringsAsFactors = FALSE)
  say("")
}

tab <- base::do.call(rbind, rows)
ps <- base::do.call(rbind, per_sample)
write.csv(tab, file.path(OUT_DIR, "bimodality_summary.csv"), row.names = FALSE)
write.csv(ps, file.path(OUT_DIR, "bimodality_per_sample.csv"), row.names = FALSE)
write.csv(ps[ps$flagged, c("gse", "gsm", "frac_mid", "group")],
          file.path(OUT_DIR, "flagged_samples.csv"), row.names = FALSE)

say("=== Summary ===")
print(tab)
say("")
say("  Suspect in the ten sets together: ", base::sum(tab$n_flag), " / ", base::sum(tab$n),
    " (samples in the model)",
    " (", sprintf("%.1f%%", 100 * base::sum(tab$n_flag) / base::sum(tab$n)), ")")
say("  Differences in median level between cohorts: ", sprintf("%.4f to %.4f", base::min(tab$median_frac),
                                        base::max(tab$median_frac)))
say("  Note: the median level differs between cohorts anyway (platform, normalisation and submitter processing differ), ")
say("        so the criterion uses the within-cohort median and median absolute deviation, without a uniform threshold.")
say("")
say("  Cohorts with completed group test: ", base::sum(tab$group_tested), " / ", base::nrow(tab))
if (base::any(!tab$group_tested))
  say("  [NOT TESTED] ", base::paste(tab$gse[!tab$group_tested], collapse = ", "),
      " - the conclusion does not cover these sets")
bad <- tab[!base::is.na(tab$fisher_p) & tab$fisher_p < 0.05, ]
if (base::nrow(bad)) {
  say("  [ACTION NEEDED] suspect samples are related to group: ",
      base::paste(sprintf("%s(p=%.3g)", bad$gse, bad$fisher_p), collapse = ", "),
      "; removal would change the composition of the two groups and must be judged separately.")
} else if (base::all(tab$group_tested)) {
  say("  In the cohorts tested, suspect samples are unrelated to group; removal amounts to a random reduction in sample size.")
} else {
  say("  No relation to group in the cohorts tested; the sets not tested are outside this conclusion.")
}

# The criterion is relative: it identifies the tail of each cohort's distribution, not absolute quality.
# Listing "the worst sample relative to the cohort median" as well shows which sets really have bad samples.
say("")
say("=== How much worse are the suspect samples: worst value relative to the cohort median ===")
o <- base::order(-tab$worst_vs_median)
for (i in o)
  say(sprintf("  %-10s median %.4f  worst %.4f  higher by %+.0f%%  suspect %d",
              tab$gse[i], tab$median_frac[i], tab$worst[i],
              100 * (tab$worst_vs_median[i] - 1), tab$n_flag[i]))
say("  Interpretation: values higher by 20-30%, i.e. the tail of the distribution, are not bad samples; only values more than twice as high")
say("        are qualitatively different. Uniform removal by the relative criterion would discard many good samples to remove a few truly bad ones.")

say("")
say("=== Plotting ===")
ps$panel <- factor(base::paste0(ps$gse, "\n", ps$disease[base::match(ps$gse, tab$gse)]),
                   levels = base::unique(base::paste0(tab$gse, "\n", tab$disease)))
p <- ggplot(ps, aes(frac_mid, panel, colour = flagged)) +
  geom_point(position = position_jitter(height = 0.18), size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = PAL[2]),
                      labels = c("passes", "flagged"), name = NULL) +
  labs(title = "Sample-level bimodality across the ten cohorts",
       subtitle = "Fraction of autosomal probes with beta between 0.3 and 0.7; higher means a blurrier signal",
       x = "Fraction of probes with intermediate beta", y = NULL) +
  theme_paper() + theme(axis.text.y = element_text(size = 6))
save_fig(p, "fig46_bimodality_qc", width = 9, height = 6.5)

say("")
say("Next, to be decided manually: whether suspect samples are removed. If so, the EWAS of the affected cohorts must be rebuilt ")
say("and the pairs involving them recomputed; the expected direction is that within-disease concordance may rise.")
log_stop()
