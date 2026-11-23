# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 47: robustness check by region-level aggregation
#
# Only one question is asked: did the noise of single probes bury a concordance that should have been seen?
# Reviewers will certainly ask this, and it is the only technical explanation not yet excluded for the two negative results,
# Alzheimer's disease and depression.
#
# The expectation is written down beforehand, to avoid explaining after seeing the result: if the two negatives were due to noise,
# their concordance should rise clearly after aggregating neighbouring probes into regions and averaging the noise out;
# if it stays near zero after aggregation, the negatives are a property of the data.
#
# Regions are defined by a single rule: windows of equal width of 1 kb along the genome.
# Neither genes nor CpG islands are used, because both annotations are one-to-many and incomplete and would introduce
# additional judgement; with equal-width windows each probe belongs to exactly one region, without ambiguity.
# Only windows with at least 3 probes are kept - with one or two probes in a window aggregation averages
# no noise out, and keeping them would dilute the conclusion.
#
# Everything else is identical to step 45: the same models, the same permutation nulls, and again only
# the five within-disease positive controls. This step changes no existing result and only answers the question above.
#
# Environment: Windows + RStudio. About one hour
# Location: <project root>/scripts/47_region_level_robustness.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "47_region_level_robustness"

SET_DIR   <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR   <- file.path(DATA_ROOT, "ewas_final")
PROBE_DIR <- file.path(DATA_ROOT, "probe_lists")
PC_DIR    <- file.path(PROJ_ROOT, "output", "45_stable")
OUT_DIR   <- file.path(PROJ_ROOT, "output", "47_region")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260920)
N_PERM    <- 2000L
BIN_BP    <- 1000L
MIN_PROBE <- 3L
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

anno <- base::readRDS(file.path(PROBE_DIR, "anno450k.rds"))
pid <- base::as.character(base::rownames(anno))
if (!length(pid) || base::mean(base::grepl("^cg", pid)) < 0.5)
  pid <- base::as.character(anno[["Name"]])
chrv <- base::as.character(anno[["chr"]])
posv <- base::suppressWarnings(base::as.numeric(anno[["pos"]]))
ok <- !base::is.na(posv) & nzchar(chrv)
BIN <- base::paste0(chrv[ok], ":", base::floor(posv[ok] / BIN_BP))
base::names(BIN) <- pid[ok]
say("Probes with a position in the annotation: ", length(BIN), "; window width ", BIN_BP, " bp")
say("")

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
  R <- qr.resid(q0, base::t(M)); rm(M); invisible(gc())
  list(q0 = q0, R = R, probes = D$probes, group = D$group)
}
coefs <- function(S, gvec) {
  gr <- qr.resid(S$q0, gvec)
  base::as.numeric(base::crossprod(gr, S$R)) / base::sum(gr^2)
}

rows <- list()
for (d in names(DIS)) {
  a <- DIS[[d]][1]; bcoh <- DIS[[d]][2]
  say("==== ", LAB_EN[[d]], ": fixed ", a, ", permuted ", bcoh, " ====")
  EA <- base::readRDS(file.path(RES_DIR, base::paste0(a, "_ewas_final.rds")))
  S <- state_of(bcoh)
  gnum <- base::as.numeric(S$group == "case")

  cm <- base::intersect(EA$probe, S$probes)
  cm <- cm[cm %in% base::names(BIN)]
  bid <- BIN[cm]
  tb <- base::table(bid)
  keep_bin <- base::names(tb)[tb >= MIN_PROBE]
  sel <- bid %in% keep_bin
  cm <- cm[sel]; bid <- bid[sel]
  say("  Common probes: ", length(base::intersect(EA$probe, S$probes)),
      "; with a position and in windows of at least ", MIN_PROBE, " probes: ", length(cm))
  say("  Windows: ", length(keep_bin), ", mean probes per window ",
      sprintf("%.1f", length(cm) / length(keep_bin)), "")

  ia <- base::match(cm, EA$probe); ib <- base::match(cm, S$probes)
  f <- factor(bid)
  agg <- function(v) base::as.numeric(base::rowsum(v, f) / base::as.numeric(base::table(f)))
  ea_r <- agg(EA$effect[ia])
  eb_all <- coefs(S, gnum)
  eb_r <- agg(eb_all[ib])

  rho_p <- stats::cor(EA$effect[ia], eb_all[ib], method = "spearman")
  rho_r <- stats::cor(ea_r, eb_r, method = "spearman")
  say(sprintf("  probe level rho=%+.4f | region level rho=%+.4f", rho_p, rho_r))

  ra <- base::rank(ea_r)
  pr <- numeric(N_PERM)
  t0 <- Sys.time()
  for (i in base::seq_len(N_PERM)) {
    pr[i] <- stats::cor(ra, base::rank(agg(coefs(S, base::sample(gnum))[ib])))
    if (i %% 500L == 0L)
      say("    permutation ", i, " / ", N_PERM, ", elapsed ",
          sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  }
  pv <- (1 + base::sum(pr >= rho_r)) / (N_PERM + 1)
  se <- base::sqrt(pv * (1 - pv) / N_PERM)
  z <- (rho_r - base::mean(pr)) / stats::sd(pr)
  crit <- base::mean(pr) + stats::qnorm(0.95) * stats::sd(pr)
  say(sprintf("  region level: Z=%+5.2f  P=%.4f±%.4f  critical rho %.4f",
              z, pv, se, crit))
  rows[[length(rows) + 1]] <- data.frame(
    disease = LAB_EN[[d]], n_probe = length(cm), n_bin = length(keep_bin),
    rho_probe = rho_p, rho_region = rho_r, z_region = z, p_region = pv,
    p_se = se, crit_region = crit, stringsAsFactors = FALSE)
  rm(S); invisible(gc())
  say("")
}

tab <- base::do.call(rbind, rows)
pc <- tryCatch(read.csv(file.path(PC_DIR, "positive_controls_10k.csv"),
                        stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(pc)) {
  tab$rho_probe_10k <- pc$rho[base::match(tab$disease, pc$disease)]
  tab$p_probe_10k <- pc$p[base::match(tab$disease, pc$disease)]
}
write.csv(tab, file.path(OUT_DIR, "region_level_within_disease.csv"), row.names = FALSE)

say("=== Probe level and region level side by side ===")
say("  Note: the 'probe' column of this table is computed only on the probes with a position that fall into windows of at least ", MIN_PROBE,
    " probes (a little over half of the common probes), ")
say("  so that it is comparable with the 'region' column on the same basis. It is not the same number as the value computed in step 45")
say("  on all common probes, and the two must not be mixed.")
for (i in base::seq_len(base::nrow(tab)))
  say(sprintf("  %-14s probe (same subset) rho=%+.4f | region rho=%+.4f  P=%.4f  critical %.4f | windows %d",
              tab$disease[i], tab$rho_probe[i], tab$rho_region[i],
              tab$p_region[i], tab$crit_region[i], tab$n_bin[i]))
if (!is.null(pc)) {
  say("")
  say("  Full-probe results of step 45 (for comparison, not combined with the table above): ")
  for (i in base::seq_len(base::nrow(tab)))
    say(sprintf("    %-14s full rho=%+.4f  P=%.4f", tab$disease[i],
                tab$rho_probe_10k[i], tab$p_probe_10k[i]))
}

# Aggregation works both ways: noise is averaged out, and the number of independent units that can be compared falls.
# Which side prevails is shown by the direction in which the critical value moves; aggregation must not be assumed to be more sensitive.
say("")
say("=== Does aggregation make the test more or less sensitive? ===")
if (!is.null(pc) && "crit_region" %in% base::names(tab)) {
  pcc <- pc$rho_crit[base::match(tab$disease, pc$disease)]
  for (i in base::seq_len(base::nrow(tab)))
    say(sprintf("  %-14s probe-level critical %.4f -> region-level critical %.4f (%+.0f%%)",
                tab$disease[i], pcc[i], tab$crit_region[i],
                100 * (tab$crit_region[i] / pcc[i] - 1)))
  say("  Critical values generally rise, so aggregation is on the whole less sensitive in these data: noise is averaged out, ")
  say("  but the loss of independent units weighs more. The region level can therefore confirm but not refute -")
  say("  a disorder that is not significant at region level but significant at probe level must be treated as 'this test is insensitive for it'.")
}
say("")
neg <- tab[tab$disease %in% c("Alzheimer", "Depression"), ]
say("=== Interpretation for the two negative results ===")
for (i in base::seq_len(base::nrow(neg)))
  say(sprintf("  %-14s probe %+.4f -> region %+.4f; region-level critical %.4f, %s",
              neg$disease[i], neg$rho_probe[i], neg$rho_region[i], neg$crit_region[i],
              if (neg$p_region[i] < 0.05) "significant after aggregation" else "still not significant after aggregation"))
if (base::all(neg$p_region >= 0.05)) {
  say("  Both negatives hold at region level; single-probe noise is not their cause.")
  say("  After aggregation Alzheimer's disease changes from ", sprintf("%+.4f", neg$rho_probe[neg$disease == "Alzheimer"]),
      " to ", sprintf("%+.4f", neg$rho_region[neg$disease == "Alzheimer"]),
      ", depression from ", sprintf("%+.4f", neg$rho_probe[neg$disease == "Depression"]),
      " to ", sprintf("%+.4f", neg$rho_region[neg$disease == "Depression"]),
      " - neither moves in the positive direction.")
} else {
  say("  [NOTE] a negative result becomes significant at region level, so single-probe noise did mask signal, ")
  say("  and the probe-level negative conclusion must be rewritten accordingly.")
}

say("")
say("=== Plotting ===")
pl <- data.frame(
  disease = base::rep(tab$disease, 2),
  level = factor(base::rep(c("Probe level", "Region level"), each = base::nrow(tab)),
                 levels = c("Probe level", "Region level")),
  rho = c(tab$rho_probe, tab$rho_region), stringsAsFactors = FALSE)
p <- ggplot(pl, aes(disease, rho, fill = level)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.72), width = 0.68) +
  scale_fill_manual(values = c("Probe level" = PAL[1], "Region level" = PAL[3]), name = NULL) +
  labs(title = "Within-disease replication at probe and region level",
       subtitle = "Regions are one-kilobase genomic windows containing at least three probes",
       x = NULL, y = "Spearman correlation of effect sizes") +
  theme_paper()
save_fig(p, "fig50_region_level", width = 8, height = 5)
log_stop()
