# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 23: apply the common practice to our data and see what it returns
#
# At this point the sharing between PD and IBD is negative under three independent tests:
#   site level (concordance of CpG effect sizes)   rho = -0.011, P = 0.62
#   composition level (twelve types, clr)          rho = +0.503, P = 0.065 (unstable)
#   composition level (lymphocyte subcomposition)  rho = -0.095, P = 0.60
# whereas the same criterion is consistently positive between the two independent IBD cohorts
#   (P = 0.020 / 0.0029 / 0.061 under the three treatments),
# so the criterion has power and the negatives are not due to a method that cannot detect.
#
# The question then is: what does the common practice of the comorbidity literature - take the differential sites of each disease,
# intersect them, and treat the genes in the intersection as a "shared molecular signature" - return on the same data?
# This step does exactly that, with our own results, following that workflow.
#
# If the common practice still produces hundreds of "shared sites" and a plausible gene list for a pair of diseases
# judged not to share by three tests, that is a positive result that can be shown directly:
# the criterion of such papers cannot fail and therefore carries no information.
#
# The same is done for the negative-control disease - if the intersections of schizophrenia with both
# are equally large, the practice does not even have specificity.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/23_template_falsepositive.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "23_template_falsepositive"

RES_DIR <- file.path(DATA_ROOT, "ewas_final")
REF_DIR <- file.path(DATA_ROOT, "probe_lists")
OUT_DIR <- file.path(PROJ_ROOT, "output", "23_template")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
TOPN <- c(500L, 1000L, 5000L, 10000L)     # The "take the top N" common in template papers
PCUT <- c(0.05, 0.01, 0.001)              # and "select by unadjusted P value"

ROLE <- c(GSE145361 = "PD", GSE111629 = "PD (replication)",
          GSE87648 = "IBD", GSE112611 = "IBD (paediatric)",
          GSE152027 = "Schizophrenia")
PAIRS <- list(
  list(id = "PD x IBD",            a = "GSE145361", b = "GSE87648",  kind = "comorbid pair"),
  list(id = "PD x SCZ",            a = "GSE145361", b = "GSE152027", kind = "negative-control pair"),
  list(id = "IBD x SCZ",           a = "GSE87648",  b = "GSE152027", kind = "negative-control pair"),
  list(id = "IBD x IBD(paed)",     a = "GSE87648",  b = "GSE112611", kind = "within-disease positive pair"))

ew <- list()
for (g in names(ROLE)) {
  f <- file.path(RES_DIR, paste0(g, "_ewas_final.rds"))
  if (file.exists(f)) ew[[g]] <- readRDS(f)
}
say("EWAS results loaded: ", length(ew), "")

anno <- NULL
af <- file.path(REF_DIR, "anno450k.rds")
if (file.exists(af)) {
  anno <- readRDS(af)
  say("Annotation loaded: ", nrow(anno), " rows, used to map sites to genes")
} else say("Annotation file not found; the gene-level part is skipped")
say("")

gene_of <- function(probes) {
  if (is.null(anno)) return(character(0))
  g <- anno[intersect(probes, rownames(anno)), "UCSC_RefGene_Name"]
  g <- unlist(strsplit(as.character(g), ";", fixed = TRUE))
  unique(g[nzchar(g)])
}

rows <- list()
for (P in PAIRS) {
  if (is.null(ew[[P$a]]) || is.null(ew[[P$b]])) next
  A <- ew[[P$a]]; B <- ew[[P$b]]
  common <- intersect(A$probe, B$probe)
  ia <- match(common, A$probe); ib <- match(common, B$probe)
  pa <- A$p[ia]; pb <- B$p[ib]
  ea <- A$effect[ia]; eb <- B$effect[ib]
  N <- length(common)
  say("==== ", P$id, " (", P$kind, ")  common probes ", N, " ====")

  for (n in TOPN) {
    sa <- common[order(pa)][seq_len(n)]
    sb <- common[order(pb)][seq_len(n)]
    ov <- intersect(sa, sb)
    exp_ov <- n * n / N
    idx <- match(ov, common)
    same_dir <- if (length(idx)) mean(sign(ea[idx]) == sign(eb[idx])) else NA_real_
    ng <- length(gene_of(ov))
    say(sprintf("  top %5d each: overlap %4d (random expectation %5.1f, fold %4.2f), same direction %s, genes involved %d",
                n, length(ov), exp_ov, length(ov) / exp_ov,
                if (is.na(same_dir)) "NA" else sprintf("%.2f", same_dir), ng))
    rows[[length(rows) + 1]] <- data.frame(
      pair = P$id, kind = P$kind, rule = paste0("top", n),
      n_sel = n, n_overlap = length(ov), expected = exp_ov,
      fold = length(ov) / exp_ov, same_dir = same_dir, n_gene = ng,
      stringsAsFactors = FALSE)
  }
  for (p0 in PCUT) {
    sa <- common[pa < p0]; sb <- common[pb < p0]
    ov <- intersect(sa, sb)
    exp_ov <- length(sa) * length(sb) / N
    idx <- match(ov, common)
    same_dir <- if (length(idx)) mean(sign(ea[idx]) == sign(eb[idx])) else NA_real_
    ng <- length(gene_of(ov))
    say(sprintf("  P < %-6g each: %6d and %6d sites, overlap %4d (expectation %5.1f, fold %4.2f), same direction %s, genes %d",
                p0, length(sa), length(sb), length(ov), exp_ov,
                length(ov) / exp_ov,
                if (is.na(same_dir)) "NA" else sprintf("%.2f", same_dir), ng))
    rows[[length(rows) + 1]] <- data.frame(
      pair = P$id, kind = P$kind, rule = paste0("P<", p0),
      n_sel = min(length(sa), length(sb)), n_overlap = length(ov),
      expected = exp_ov, fold = length(ov) / exp_ov,
      same_dir = same_dir, n_gene = ng, stringsAsFactors = FALSE)
  }
  say("")
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "template_overlaps.csv"), row.names = FALSE)

# --------------------------------------------------------------
say("=== A 'shared gene list' that could go straight into a paper ===")
P <- PAIRS[[1]]
A <- ew[[P$a]]; B <- ew[[P$b]]
common <- intersect(A$probe, B$probe)
ia <- match(common, A$probe); ib <- match(common, B$probe)
sa <- common[order(A$p[ia])][seq_len(5000L)]
sb <- common[order(B$p[ib])][seq_len(5000L)]
ov <- intersect(sa, sb)
say("  By the most common practice (intersecting the top 5000 sites of each disease): ")
say("  Shared sites: ", length(ov), ", genes involved: ", length(gene_of(ov)), "")
gl <- gene_of(ov)
if (length(gl)) {
  say("  First 30 genes: ", paste(utils::head(sort(gl), 30), collapse = ", "))
  writeLines(sort(gl), file.path(OUT_DIR, "pseudo_shared_genes.txt"), useBytes = TRUE)
}
say("  Reminder: this list comes from a pair of diseases judged not to share by three tests.")
say("        It looks entirely capable of supporting a paper, and that is exactly the problem.")
say("")

# --------------------------------------------------------------
say("=== Plotting ===")
tab$panel <- factor(tab$pair, levels = vapply(PAIRS, function(x) x$id, ""))
sub <- tab[grepl("^top", tab$rule), ]
sub$rule <- factor(sub$rule, levels = paste0("top", TOPN))
p <- ggplot(sub, aes(rule, n_overlap, fill = kind)) +
  geom_col(width = 0.68) +
  geom_point(aes(y = expected), shape = 4, size = 2, colour = "grey25") +
  geom_text(aes(label = n_overlap), vjust = -0.4, size = 2.5) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c("comorbid pair" = PAL[2], "negative-control pair" = PAL[8],
                               "within-disease positive pair" = PAL[3]), name = NULL) +
  labs(title = "What the conventional overlap approach returns on these data",
       subtitle = "Bars are the observed overlap of top-ranked probes; crosses mark the overlap expected by chance",
       x = "Number of top probes taken from each disease",
       y = "Overlapping probes") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
save_fig(p, "fig30_template_overlap", width = 10, height = 4.6)

p <- ggplot(sub, aes(rule, fold, fill = kind)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey45", linewidth = 0.4) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.2f", fold)), vjust = -0.4, size = 2.5) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = c("comorbid pair" = PAL[2], "negative-control pair" = PAL[8],
                               "within-disease positive pair" = PAL[3]), name = NULL) +
  labs(title = "Overlap relative to chance",
       subtitle = "A value of one means the overlap is exactly what independent selections would produce",
       x = "Number of top probes taken from each disease",
       y = "Observed overlap divided by expected") +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
save_fig(p, "fig31_template_fold", width = 10, height = 4.6)

log_stop()
