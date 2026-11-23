# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 27: external anchoring of the whole pipeline with published summary statistics
#
# This is the only external check of the project and was postponed for the longest time.
#
# The whole-blood methylation meta-analysis of Parkinson's disease published by Vallerga et al. in Nature Communications in 2020
# released the effect size and P value of every probe, and GSE145361 is
# one of its sample sets. The question is direct: do the effect sizes that we computed from
# the raw matrix agree with the published ones?
#
#   agree    - our QC, covariate specification, model and handling of the inflation factor have no
#             systematic bias, and the weak PD signal is a property of the data, not an artefact of our computation.
#   disagree - something in the pipeline is wrong and all earlier conclusions must be re-examined.
#
# Note that this is not independent replication: their meta-analysis includes this cohort, so what is compared is "can the
# published result be reproduced", not "does it hold in new data". The two must not be confused.
#
# They used the mixed linear models of OSCA (MOA / MOMENT), which are on a different scale from the linear
# model of limma, so rank correlation and sign concordance are primary and magnitudes are for reference only.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/27_external_anchor.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "27_external_anchor"

RES_DIR <- file.path(DATA_ROOT, "ewas_final")
VAL_DIR <- file.path(DATA_ROOT, "VALLERGA")
OUT_DIR <- file.path(PROJ_ROOT, "output", "27_anchor")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
suppressPackageStartupMessages(library(data.table))

FILES <- c(MOA = "Vallerga2020_NCOMMS_MWAS-meta-analysis-MOA.mlma",
           MOMENT = "Vallerga2020_NCOMMS_MWAS-meta-analysis-MOMENT.mlma")
OURS <- c(GSE145361 = "PD discovery", GSE111629 = "PD replication",
          GSE72774 = "PD sensitivity", GSE87648 = "IBD main (as negative control)")

ext <- list()
for (k in names(FILES)) {
  f <- file.path(VAL_DIR, FILES[[k]])
  if (!file.exists(f)) { say("File missing: ", f); next }
  d <- data.table::fread(f)
  say(k, ": ", nrow(d), " rows, column names ", paste(names(d), collapse = ", "))
  pc <- grep("^probe$|^Probe$", names(d), value = TRUE)[1]
  bc <- grep("^b$|^beta$|^Effect$", names(d), value = TRUE)[1]
  ppc <- grep("^p$|^P$|^pval", names(d), value = TRUE)[1]
  if (is.na(pc) || is.na(bc)) { say("  Probe or effect size column not found, skipped"); next }
  ext[[k]] <- data.frame(probe = as.character(d[[pc]]),
                         b = as.numeric(d[[bc]]),
                         p = if (is.na(ppc)) NA_real_ else as.numeric(d[[ppc]]),
                         stringsAsFactors = FALSE)
  say("  Effect size range ", sprintf("%.4g to %.4g", min(ext[[k]]$b, na.rm = TRUE),
                               max(ext[[k]]$b, na.rm = TRUE)))
}
say("")

rows <- list(); scat <- list()
for (k in names(ext)) {
  E <- ext[[k]]
  for (g in names(OURS)) {
    f <- file.path(RES_DIR, paste0(g, "_ewas_final.rds"))
    if (!file.exists(f)) next
    O <- readRDS(f)
    cm <- intersect(O$probe, E$probe)
    if (length(cm) < 1000L) { say("  ", k, " vs ", g, ": only ", length(cm), " common probes, skipped"); next }
    i1 <- match(cm, O$probe); i2 <- match(cm, E$probe)
    a <- O$effect[i1]; b <- E$b[i2]
    rho <- stats::cor(a, b, method = "spearman", use = "complete.obs")
    r <- stats::cor(a, b, use = "complete.obs")
    dir <- mean(sign(a) == sign(b), na.rm = TRUE)
    # Only their most significant sites, for a higher signal-to-noise ratio
    sel <- if (all(is.na(E$p))) integer(0) else which(E$p[i2] < 1e-4)
    rho_t <- if (length(sel) >= 30L)
      stats::cor(a[sel], b[sel], method = "spearman", use = "complete.obs") else NA_real_
    dir_t <- if (length(sel) >= 30L) mean(sign(a[sel]) == sign(b[sel]), na.rm = TRUE) else NA_real_
    say(sprintf("  %-7s vs %-10s common probes %6d  rho=%+.4f  r=%+.4f  sign %.3f | their %4d significant sites: rho=%s sign %s",
                k, g, length(cm), rho, r, dir, length(sel),
                if (is.na(rho_t)) "NA" else sprintf("%+.3f", rho_t),
                if (is.na(dir_t)) "NA" else sprintf("%.3f", dir_t)))
    rows[[length(rows) + 1]] <- data.frame(
      external = k, ours = g, role = OURS[[g]], n_common = length(cm),
      rho = rho, pearson = r, dir = dir, n_sig_ext = length(sel),
      rho_top = rho_t, dir_top = dir_t, stringsAsFactors = FALSE)
    if (k == "MOA") {
      s <- sample(seq_along(cm), min(40000L, length(cm)))
      scat[[g]] <- data.frame(ours = g, role = OURS[[g]],
                              x = b[s], y = a[s], stringsAsFactors = FALSE)
    }
  }
  say("")
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "external_anchor.csv"), row.names = FALSE)

say("=== Interpretation ===")
say("  GSE145361 and GSE111629 are component samples of that meta-analysis and should show a clear positive correlation; ")
say("  GSE87648 (inflammatory bowel disease) is unrelated and should be close to zero - it serves as the negative control here; ")
say("  if it correlates as well, what agrees is technical structure rather than the effect of Parkinson's disease.")
say("")

say("=== Plotting ===")
if (length(scat)) {
  sd2 <- do.call(rbind, scat)
  sd2$panel <- factor(paste0(sd2$ours, "\n", sd2$role),
                      levels = unique(paste0(sd2$ours, "\n", sd2$role)))
  p <- ggplot(sd2, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point_big(size = 0.3, alpha = 0.2, colour = PAL[1]) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = PAL[2], linewidth = 0.5) +
    facet_wrap(~ panel, nrow = 1, scales = "free") +
    labs(title = "Our effect estimates against the published meta-analysis",
         subtitle = "Published estimates come from a mixed linear model, so the scales differ; the question is whether the ordering agrees",
         x = "Published effect (Vallerga 2020, MOA)",
         y = "Effect estimated by this pipeline") +
    theme_paper()
  save_fig(p, "fig36_external_anchor", width = 10, height = 4.2)
}
if (!is.null(tab)) {
  tab$panel <- factor(paste0(tab$ours, "\n", tab$role),
                      levels = unique(paste0(tab$ours, "\n", tab$role)))
  p <- ggplot(tab, aes(panel, rho, fill = external)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_col(position = position_dodge(0.75), width = 0.68) +
    geom_text(aes(label = sprintf("%.3f", rho)),
              position = position_dodge(0.75), vjust = -0.4, size = 2.4) +
    scale_fill_manual(values = c(MOA = PAL[1], MOMENT = PAL[3]), name = NULL) +
    labs(title = "Agreement with the published Parkinson's meta-analysis",
         subtitle = "The inflammatory bowel disease cohort serves as a negative control and should sit at zero",
         x = NULL, y = "Spearman correlation of effect sizes") +
    theme_paper()
  save_fig(p, "fig37_anchor_summary", width = 8.5, height = 4.4)
}
log_stop()
