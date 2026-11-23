# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 42: samples with discordant sex in GSE87648 - borderline misclassification or true mismatch?
#
# Result of 41: the sex check is possible in six of the ten sets; five have a concordance of 1.0000 and GSE87648 has
# 0.9447 (21 of 380 samples discordant). The shuffled control is always between 0.51 and 0.55, so the test
# has ample discriminating power; five sets reach 1.0000, so 0.9447 is not a ceiling of the method.
# Moreover, 201 samples are inferred female against 184 recorded, 17 more in one direction - random borderline misclassification
# should roughly cancel between the two directions; a one-sided shift is systematic.
#
# This data set is currently the strongest link: it carries the within-disease positive control of inflammatory bowel disease (rho=0.116),
# and the only two significant cells between schizophrenia and inflammatory bowel disease both involve it. It
# must therefore be clarified; until then both conclusions are provisional.
#
# This step uses two independent lines of evidence:
#
#   (1) X chromosome: women have one inactivated X whose CpG islands are intermediately methylated, so the proportion of
#       probes with intermediate methylation is high; men have a single X and most probes are near 0 or 1.
#   (2) Y chromosome: men have a Y and the probes give real signal; women have none and the probes show only
#       cross-hybridisation background, with low values and small variance.
#
#   A true mismatch is called only when both lines point against the record; if only one does, or the value
#   lies between the two clusters, the sample is called suspect, not mismatched. One line of evidence is not enough to overturn the record.
#
# If true mismatches are confirmed, a list of samples for removal is written; whether to remove them and re-run is decided manually,
# and this step changes no existing result.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/42_sex_mismatch_gse87648.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "42_sex_mismatch_gse87648"

RDS_DIR   <- file.path(DATA_ROOT, "rds")
RES_DIR   <- file.path(DATA_ROOT, "ewas_final")
PROBE_DIR <- file.path(DATA_ROOT, "probe_lists")
MAP_DIR   <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR   <- file.path(PROJ_ROOT, "output", "42_sexmismatch")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

GSE <- "GSE87648"
AMBIG <- 0.25    # Samples between the two cluster centres and farther than this proportion from the nearest centre are called suspect

anno <- base::readRDS(file.path(PROBE_DIR, "anno450k.rds"))
pid <- base::as.character(base::rownames(anno))
if (!length(pid) || base::mean(base::grepl("^cg", pid)) < 0.5)
  pid <- base::as.character(anno[["Name"]])
chrv <- base::as.character(anno[["chr"]])
chrX <- pid[chrv == "chrX"]; chrY <- pid[chrv == "chrY"]
say("X probes in the annotation: ", length(chrX), ", Y probes ", length(chrY))

D <- base::readRDS(file.path(RES_DIR, base::paste0(GSE, "_design.rds")))
sc <- base::grep("^sex", base::colnames(D$design), value = TRUE)[1]
lab_pos <- base::sub("^sex", "", sc)
sex_rec <- base::ifelse(D$design[, sc] == 1, lab_pos, "REF")
say("Recorded sex: ", base::paste(sprintf("%s=%d", base::names(base::table(sex_rec)),
                                        base::as.integer(base::table(sex_rec))),
                                collapse = "; "))

b <- base::readRDS(file.path(RDS_DIR, base::paste0(GSE, "_beta.rds")))
m <- read.csv(file.path(MAP_DIR, base::paste0(GSE, "_sample_map.csv")),
              stringsAsFactors = FALSE, colClasses = "character")
gs <- m$gsm[base::match(base::colnames(b), m$col_name)]
base::colnames(b) <- gs
b <- b[, !base::is.na(gs), drop = FALSE]
b <- b[, base::intersect(D$gsm, base::colnames(b)), drop = FALSE]
say("Matrix ", base::nrow(b), " x ", base::ncol(b))

ix <- base::which(base::rownames(b) %in% chrX)
iy <- base::which(base::rownames(b) %in% chrY)
say("X probes in the matrix: ", length(ix), ", Y probes ", length(iy))

fx <- base::colMeans(b[ix, , drop = FALSE] > 0.3 & b[ix, , drop = FALSE] < 0.7, na.rm = TRUE)
have_y <- length(iy) >= 100L
if (have_y) {
  my <- base::colMeans(b[iy, , drop = FALSE], na.rm = TRUE)
  sy <- base::apply(b[iy, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  say("Y probes available; their mean and SD serve as the second line of evidence")
} else say("[NOTE] fewer than 100 Y probes in the matrix; only the X line is available and no two-line decision is possible")
rm(b); invisible(gc())

split2 <- function(v, nm) {
  cl <- stats::kmeans(v, centers = 2, nstart = 25)
  ctr <- base::sort(cl$centers[, 1])
  d1 <- base::abs(v - ctr[1]); d2 <- base::abs(v - ctr[2])
  near <- base::pmin(d1, d2) / base::abs(base::diff(ctr))
  say("  ", nm, ": cluster centres ", sprintf("%.3f and %.3f", ctr[1], ctr[2]),
      ", distance ", sprintf("%.3f", base::diff(ctr)),
      "; between the two clusters: ", base::sum(near > AMBIG), " samples")
  list(hi = v > base::mean(ctr), near = near, ctr = ctr)
}
SEP_MIN_Y <- 0.15   # Lower limit of the distance between the two Y clusters: men and women should differ greatly on Y;
                    # a small distance means a unimodal distribution is being cut in two
say("")
say("=== Clustering by each line of evidence ===")
X <- split2(fx, "proportion of X probes with intermediate methylation")
pred_f_x <- X$hi                       # high = female
if (have_y) {
  Y <- split2(my, "mean of Y probes")
  sepY <- base::abs(base::diff(Y$ctr))
  if (sepY < SEP_MIN_Y) {
    say("  [FAIL] the distance between the two Y clusters is only ", sprintf("%.3f", sepY), ", below ", SEP_MIN_Y,
        ". Women have no Y and only cross-hybridisation background remains, so men and women should differ greatly; ")
    say("            such a small distance means the Y probes of this matrix have no discriminating power, and this line of evidence is void.")
    have_y <- FALSE
  } else pred_f_y <- !Y$hi             # low = female
}
if (!have_y)
  say("  Only the X line is available; the decisions below rest on it and are not described as 'two lines of evidence'")

rec_f <- sex_rec[base::match(base::names(fx), D$gsm)] == "REF"
n_ref <- base::sum(rec_f)
say("")
say("Recorded 'REF': ", n_ref, "; called female by the X line: ", base::sum(pred_f_x), " samples",
    if (have_y) base::paste0("; called female by the Y line: ", base::sum(pred_f_y), " samples") else "")
say(" (first confirming that REF corresponds to female: concordance with the X line ",
    sprintf("%.4f", base::mean(pred_f_x == rec_f)), ")")

res <- data.frame(gsm = base::names(fx), recorded = base::ifelse(rec_f, "F", "M"),
                  x_frac_mid = base::as.numeric(fx),
                  x_pred = base::ifelse(pred_f_x, "F", "M"),
                  x_borderline = X$near > AMBIG, stringsAsFactors = FALSE)
if (have_y) {
  res$y_mean <- base::as.numeric(my); res$y_sd <- base::as.numeric(sy)
  res$y_pred <- base::ifelse(pred_f_y, "F", "M")
  res$y_borderline <- Y$near > AMBIG
}

if (have_y) {
  res$verdict <- base::ifelse(
    res$x_pred == res$recorded & res$y_pred == res$recorded, "consistent",
    base::ifelse(res$x_pred != res$recorded & res$y_pred != res$recorded &
                 !res$x_borderline & !res$y_borderline, "true mismatch",
    base::ifelse(res$x_borderline | res$y_borderline, "suspect (intermediate value)", "suspect (the two lines of evidence disagree)")))
} else {
  # X only: "X clearly discordant" is called only when the value lies clearly on the side of the opposite-sex cluster; intermediate values are called suspect.
  # A single line of evidence cannot establish a sample mismatch; such samples are only flagged for manual confirmation.
  res$verdict <- base::ifelse(res$x_pred == res$recorded, "consistent",
                 base::ifelse(res$x_borderline, "suspect (intermediate value)", "X clearly discordant"))
}
tb <- base::table(res$verdict)
say("")
say("=== Summary of decisions ===")
for (k in base::names(tb)) say(sprintf("  %-20s %d", k, tb[[k]]))

bad <- res[res$verdict %in% c("true mismatch", "X clearly discordant"), , drop = FALSE]
sus <- res[base::grepl("^suspect|X only", res$verdict), , drop = FALSE]
if (base::nrow(bad)) {
  say("")
  say(if (have_y) "True mismatches: " else "X clearly discordant: ", base::nrow(bad),
      if (have_y) " (both lines of evidence against the record, neither at the boundary): " else " (X value clearly on the side of the opposite-sex cluster): ")
  for (i in base::seq_len(base::min(30L, base::nrow(bad))))
    say(sprintf("  %-12s recorded %s  X proportion %.3f  Y mean %.3f", bad$gsm[i], bad$recorded[i],
                bad$x_frac_mid[i], if (have_y) bad$y_mean[i] else NA_real_))
}
if (base::nrow(sus)) say("Suspect: ", base::nrow(sus), "; the record is not overturned on this basis")

write.csv(res, file.path(OUT_DIR, base::paste0(GSE, "_sex_detail.csv")), row.names = FALSE)
if (base::nrow(bad))
  write.csv(bad[, c("gsm", "recorded", "verdict")],
            file.path(OUT_DIR, base::paste0(GSE, "_exclude_candidates.csv")), row.names = FALSE)

# Is the mismatch related to case-control status? If concentrated in one group, the effect is far worse than random loss
grp <- D$group[base::match(res$gsm, D$gsm)]
if (base::nrow(bad)) {
  ct <- base::table(base::ifelse(res$gsm %in% bad$gsm, "discordant", "normal"), grp)
  say("")
  say("=== Mismatch and group ===")
  for (r in base::rownames(ct))
    say("  ", r, ": ", base::paste(sprintf("%s=%d", base::colnames(ct),
                                           base::as.integer(ct[r, ])), collapse = "; "))
  if (base::all(base::dim(ct) == 2)) {
    pv <- suppressWarnings(stats::fisher.test(ct)$p.value)
    say("  Fisher's exact test: p = ", format.pval(pv, digits = 3),
        if (pv < 0.05) "  [NOTE] mismatch is related to group; removal would change the composition of the two groups"
        else "  Mismatch is unrelated to group; removal amounts to a random reduction in sample size")
  }
}

# --------------------------------------------------------------
# Samples with intermediate values: sex-chromosome anomaly or poor data quality?
# Poor-quality samples have high background, and more probes across the whole genome fall between 0.3 and 0.7;
# this is visible on the autosomes too, whereas a true sex-chromosome anomaly would not show it.
say("")
say("=== Intermediate values: sex-chromosome anomaly or data quality? ===")
b2 <- base::readRDS(file.path(RDS_DIR, base::paste0(GSE, "_beta.rds")))
base::colnames(b2) <- m$gsm[base::match(base::colnames(b2), m$col_name)]
b2 <- b2[, base::intersect(res$gsm, base::colnames(b2)), drop = FALSE]
auto <- base::which(!(base::rownames(b2) %in% c(chrX, chrY)))
set.seed(1); auto <- base::sample(auto, base::min(50000L, length(auto)))
fa <- base::colMeans(b2[auto, , drop = FALSE] > 0.3 & b2[auto, , drop = FALSE] < 0.7,
                     na.rm = TRUE)
rm(b2); invisible(gc())
res$auto_frac_mid <- base::as.numeric(fa[base::match(res$gsm, base::names(fa))])
mid <- res$verdict == "suspect (intermediate value)"
say("  Proportion of autosomal probes between 0.3 and 0.7 (quality index, higher is worse): ")
say("    consistent samples     median ", sprintf("%.4f", stats::median(res$auto_frac_mid[res$verdict == "consistent"], na.rm = TRUE)))
if (base::any(mid))
  say("    intermediate values    median ", sprintf("%.4f", stats::median(res$auto_frac_mid[mid], na.rm = TRUE)),
      " (n=", base::sum(mid), ")")
if (base::nrow(bad))
  say("    X clearly discordant   median ", sprintf("%.4f", stats::median(res$auto_frac_mid[res$gsm %in% bad$gsm], na.rm = TRUE)),
      " (n=", base::nrow(bad), ")")
if (base::any(mid) && base::sum(res$verdict == "consistent") > 30) {
  pv <- stats::wilcox.test(res$auto_frac_mid[mid],
                           res$auto_frac_mid[res$verdict == "consistent"])$p.value
  say("  Difference between the two, Wilcoxon p = ", format.pval(pv, digits = 3),
      if (pv < 0.05) "  The overall data quality of samples with intermediate values is indeed poorer; they should be treated as a quality problem"
      else "  No difference in quality; intermediate values cannot be explained by data quality")
}

say("")
say("=== Plotting ===")
res$grp <- base::as.character(grp)
if (have_y) {
  p <- ggplot(res, aes(x_frac_mid, y_mean, colour = recorded, shape = verdict)) +
    geom_point(size = 1.9, alpha = 0.85) +
    labs(title = "Sex check for the adult IBD cohort, two independent signals",
         subtitle = "X-chromosome intermediate methylation against mean Y-chromosome beta",
         x = "Fraction of chrX probes with beta between 0.3 and 0.7",
         y = "Mean beta over chrY probes", colour = "Recorded", shape = "Verdict") +
    theme_paper()
} else {
  p <- ggplot(res, aes(x_frac_mid, fill = recorded)) +
    geom_histogram(bins = 60, alpha = 0.75, position = "identity") +
    labs(title = "Sex check for the adult IBD cohort",
         subtitle = "Only the X-chromosome signal is available in this matrix",
         x = "Fraction of chrX probes with beta between 0.3 and 0.7", y = "Samples",
         fill = "Recorded") +
    theme_paper()
}
save_fig(p, "fig45_sex_check_gse87648", width = 8, height = 5.5)

say("")
say("Next, to be decided manually: whether samples confirmed as true mismatches are removed; if so, the EWAS of this cohort must be re-run, ")
say("together with the pairs involving it (one within-disease positive control and eight cross-disease pairs).")
log_stop()
