# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 10: graphical check of cohort composition and data quality
#
# This step yields no statistical conclusion. It has two purposes:
#   first, to see the composition, age and sex distributions, beta distributions, missingness and principal-component
#     structure of the six data sets at a glance; problems invisible in tables (bimodal distributions, batch strata,
#     a group crowded into a corner) only show up in figures;
#   second, several of these figures are themselves Figure 1 and supplementary material of the manuscript;
#     reviewers will certainly ask how the large age difference between the cohorts of the two diseases was handled,
#     and showing the figure is more convincing than arguing in words.
#
# Output: one 600 dpi PNG and one vector PDF per figure under output/figures/.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/10_qc_figures.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "10_qc_figures"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
OUT_DIR <- file.path(PROJ_ROOT, "output", "10_qc")
if (!dir.exists(SET_DIR)) stop("Directory does not exist: ", SET_DIR)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260910)
N_PROBE_DENSITY <- 20000L   # Number of probes sampled for the density curves
N_PROBE_PCA     <- 10000L   # Number of highly variable probes used for the principal components

# Filtered analysis objects are preferred; unfiltered ones are used if they do not exist
files <- list.files(SET_DIR, pattern = "_analysis_filtered\\.rds$", full.names = TRUE)
if (!length(files)) {
  files <- list.files(SET_DIR, pattern = "_analysis\\.rds$", full.names = TRUE)
  say("Note: filtered analysis objects not found; using the unfiltered version")
}
say("Reading ", length(files), " analysis objects")

# Role of the cohort in the paper; decides plotting order and colours
ROLE <- c(GSE145361 = "PD discovery", GSE111629 = "PD replication",
          GSE72774  = "PD sensitivity", GSE87648 = "IBD main",
          GSE112611 = "IBD paediatric", GSE152027 = "Negative control")
ORDER <- names(ROLE)

# --------------------------------------------------------------
# One pass collects the small data needed for plotting; large matrices are discarded right after reading
# --------------------------------------------------------------
pheno <- list(); dens <- list(); pcs <- list(); miss <- list()
for (f in files) {
  obj <- readRDS(f)
  g <- obj$gse
  say("Processing ", g, ": ", nrow(obj$beta), " probes x ", ncol(obj$beta), " samples")

  ph <- obj$pheno
  agecol <- grep("^age", names(ph), ignore.case = TRUE, value = TRUE)[1]
  sexcol <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  sx <- if (is.na(sexcol)) NA_character_ else toupper(substr(trimws(ph[[sexcol]]), 1, 1))
  pheno[[g]] <- data.frame(
    gse = g, disease = obj$disease,
    group = as.character(ph$group),
    age = if (is.na(agecol)) NA_real_ else suppressWarnings(as.numeric(ph[[agecol]])),
    sex = ifelse(sx %in% c("M", "F"), sx, NA_character_),
    stringsAsFactors = FALSE)

  idx <- sample(seq_len(nrow(obj$beta)), min(N_PROBE_DENSITY, nrow(obj$beta)))
  v <- as.numeric(obj$beta[idx, sample(ncol(obj$beta), min(60L, ncol(obj$beta)))])
  dens[[g]] <- data.frame(gse = g, beta = v[is.finite(v)], stringsAsFactors = FALSE)

  miss[[g]] <- data.frame(gse = g,
                          rate = colMeans(!is.finite(obj$beta)),
                          group = as.character(ph$group),
                          stringsAsFactors = FALSE)

  vr <- apply(obj$beta[sample(seq_len(nrow(obj$beta)),
                              min(60000L, nrow(obj$beta))), , drop = FALSE],
              1, stats::var, na.rm = TRUE)
  top <- names(sort(vr, decreasing = TRUE))[seq_len(min(N_PROBE_PCA, length(vr)))]
  m <- obj$beta[top, , drop = FALSE]
  m[!is.finite(m)] <- NA
  m <- m[stats::complete.cases(m), , drop = FALSE]
  pr <- stats::prcomp(t(m), center = TRUE, scale. = FALSE)
  ve <- round(100 * pr$sdev^2 / sum(pr$sdev^2), 1)
  pcs[[g]] <- data.frame(gse = g, PC1 = pr$x[, 1], PC2 = pr$x[, 2],
                         group = as.character(ph$group),
                         ve1 = ve[1], ve2 = ve[2], stringsAsFactors = FALSE)
  say("  PC1 explains ", ve[1], "%, PC2 explains ", ve[2], "%")
  rm(obj, m, pr); invisible(gc())
}

ph_all <- do.call(rbind, pheno)
ph_all$gse <- factor(ph_all$gse, levels = ORDER[ORDER %in% ph_all$gse])
ph_all$panel <- paste0(ph_all$gse, "\n", ROLE[as.character(ph_all$gse)])
ph_all$panel <- factor(ph_all$panel,
                       levels = unique(ph_all$panel[order(ph_all$gse)]))
ph_all$group <- factor(ph_all$group, levels = c("control", "case"))
write.csv(ph_all, file.path(OUT_DIR, "pheno_summary.csv"), row.names = FALSE)

GRP_COL <- c(control = PAL[1], case = PAL[2])

# --------------------------------------------------------------
say("")
say("=== Plotting ===")

# Figure 1: number of cases and controls per cohort
n_tab <- as.data.frame(table(ph_all$panel, ph_all$group))
names(n_tab) <- c("panel", "group", "n")
p1 <- ggplot(n_tab, aes(panel, n, fill = group)) +
  geom_col(width = 0.62, colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = 2.8, colour = "white") +
  scale_fill_manual(values = GRP_COL, name = NULL,
                    labels = c("Control", "Case")) +
  labs(title = "Cohort composition",
       subtitle = "Public blood 450K/EPIC datasets after QC and group assignment",
       x = NULL, y = "Samples") +
  theme_paper() + theme(axis.text.x = element_text(size = 8))
save_fig(p1, "fig01_cohort_composition", width = 8, height = 4.6)

# Figure 2: age distribution - the most important confounder in this study
p2 <- ggplot(subset(ph_all, !is.na(age)),
             aes(panel, age, fill = group)) +
  geom_violin(position = position_dodge(0.75), width = 0.85,
              linewidth = 0.25, colour = "grey30", alpha = 0.75) +
  geom_boxplot(position = position_dodge(0.75), width = 0.16,
               outlier.size = 0.3, linewidth = 0.25, fill = "white") +
  scale_fill_manual(values = GRP_COL, name = NULL,
                    labels = c("Control", "Case")) +
  labs(title = "Age distribution by cohort",
       subtitle = "PD cohorts are elderly; the IBD cohorts are young adults and children",
       caption = paste("Cohorts without an age field (GSE145361) or with mostly",
                       "missing age (GSE152027) are omitted here and handled",
                       "by DNAm age."),
       x = NULL, y = "Age (years)") +
  theme_paper() + theme(axis.text.x = element_text(size = 8))
save_fig(p2, "fig02_age_distribution", width = 8, height = 4.6)

# Figure 3: sex composition
sx <- as.data.frame(table(ph_all$panel, ph_all$group, ph_all$sex))
names(sx) <- c("panel", "group", "sex", "n")
sx <- subset(sx, n > 0)
p3 <- ggplot(sx, aes(group, n, fill = sex)) +
  geom_col(position = "fill", width = 0.68, colour = "grey25", linewidth = 0.2) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c(F = PAL[4], M = PAL[6]), name = NULL,
                    labels = c(F = "Female", M = "Male")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(labels = c(control = "Ctrl", case = "Case")) +
  labs(title = "Sex composition", x = NULL, y = "Proportion") +
  theme_paper()
save_fig(p3, "fig03_sex_composition", width = 9, height = 3.6)

# Figure 4: beta distribution - the cohorts were normalised differently, so comparability is checked first
d_all <- do.call(rbind, dens)
d_all$gse <- factor(d_all$gse, levels = ORDER[ORDER %in% d_all$gse])
p4 <- ggplot(d_all, aes(beta, colour = gse)) +
  geom_density(linewidth = 0.55, adjust = 1.2) +
  scale_colour_manual(values = PAL[seq_len(nlevels(d_all$gse))], name = NULL) +
  labs(title = "Beta-value distributions",
       subtitle = paste0("Random ", format(N_PROBE_DENSITY, big.mark = ","),
                         " probes x up to 60 samples per cohort"),
       x = "Beta value", y = "Density") +
  theme_paper()
save_fig(p4, "fig04_beta_density", width = 7, height = 4.2)

# Figure 5: sample missing rate
m_all <- do.call(rbind, miss)
m_all$gse <- factor(m_all$gse, levels = ORDER[ORDER %in% m_all$gse])
p5 <- ggplot(m_all, aes(gse, rate * 100, fill = gse)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.25, width = 0.6) +
  scale_fill_manual(values = PAL[seq_len(nlevels(m_all$gse))], guide = "none") +
  labs(title = "Per-sample missing rate",
       x = NULL, y = "Missing probes (%)") +
  theme_paper() + theme(axis.text.x = element_text(size = 8))
save_fig(p5, "fig05_sample_missingness", width = 7, height = 4)

# Figure 6: principal components per cohort - check for strata and whether case-control status is dominated by other structure
pc_all <- do.call(rbind, pcs)
pc_all$gse <- factor(pc_all$gse, levels = ORDER[ORDER %in% pc_all$gse])
pc_all$group <- factor(pc_all$group, levels = c("control", "case"))
lab <- unique(pc_all[, c("gse", "ve1", "ve2")])
pc_all$panel <- factor(paste0(pc_all$gse, "  (PC1 ", pc_all$ve1,
                              "%, PC2 ", pc_all$ve2, "%)"),
                       levels = paste0(lab$gse[order(lab$gse)], "  (PC1 ",
                                       lab$ve1[order(lab$gse)], "%, PC2 ",
                                       lab$ve2[order(lab$gse)], "%)"))
p6 <- ggplot(pc_all, aes(PC1, PC2, colour = group)) +
  geom_point_big(size = 0.7, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  scale_colour_manual(values = GRP_COL, name = NULL,
                      labels = c("Control", "Case")) +
  labs(title = "Principal components within each cohort",
       subtitle = paste0("Top ", format(N_PROBE_PCA, big.mark = ","),
                         " most variable probes; computed within cohort"),
       x = "PC1", y = "PC2") +
  theme_paper() +
  guides(colour = guide_legend(override.aes = list(size = 2)))
save_fig(p6, "fig06_pca_by_cohort", width = 9, height = 6)

say("")
say("Figures written: ", file.path(PROJ_ROOT, "output", "figures"))
log_stop()
