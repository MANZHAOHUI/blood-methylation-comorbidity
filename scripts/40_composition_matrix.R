# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 40: concordance matrix of leukocyte composition (second axis)
#
# Uses the same 45 pairs and the same permutation nulls as step 39; only the quantity compared changes
# from the effects at four hundred thousand sites to the effects on twelve leukocyte fractions. The two axes measure
# different things: the site level asks "does the same CpG change in the same direction", the composition level asks
# "do cell fractions shift in the same direction". Diseases may well share the latter and not the former.
#
# This axis has a clear basis in the literature. Hannon et al. (eLife. 2021;10:e58430) reported, in a meta-analysis of psychosis in seven
# cohorts and 4483 individuals, significant case-control differences in blood cell fractions and in smoking exposure
# derived from methylation. This step extends the observation from one disorder to five
# and adds within-disease positive controls - a layer absent from the original paper.
#
# Three methodological points:
#
#   (1) Compositional data require the centred log-ratio transform. Fractions sum to 1, an increase of one type
#       is necessarily accompanied by decreases of the others, and comparing raw fractions creates artificial sign concordance. The first round
#       fell into exactly this trap: under the raw logit all 12 cell types agreed in direction, and with clr
#       the conclusion reversed completely.
#
#   (2) Interpretation always rests on the permutation null, not on a theoretical baseline. The -1/(D-1) often cited
#       for compositional data (Aitchison's spurious correlation) concerns the average correlation between two different components across samples,
#       which is not the correlation between two effect vectors across components computed here:
#       Pearson centres first, the clr vector already sums to zero, centring is a no-op,
#       and the expected value is therefore 0. The observed null means lie near 0, in agreement.
#
#   (2b) With only 12 components the null distribution of the correlation is bounded and far from normal, so Z is a poor
#       summary (Z of only 2.4 with a permutation P of 1e-4 can occur). The permutation P is always
#       decisive and Z is for reference only.
#
#   (3) With only 12 components the resolution of a rank correlation is too low (Spearman on 12 points takes only
#       a few possible values), so Pearson is the primary criterion and Spearman is reported alongside.
#       This does not contradict the use of Spearman in 39: with four hundred thousand points rank correlation is more robust, with twelve
#       points it discards too much information. Both come with permutation nulls.
#
# An additional section: case-control difference in the methylation-derived smoking score, also with within-disease replication.
# Absolute scores are not comparable across cohorts (different numbers of usable sites shift them as a whole), but
# the standardised case-control difference within a cohort is, so only the latter is compared.
#
# Environment: Windows + RStudio. All computations are on small matrices and finish in a few minutes
# Location: <project root>/scripts/40_composition_matrix.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "40_composition_matrix"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "40_composition")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260917)
N_PERM <- 10000L          # Small matrices, so permutations can be plentiful

DIS <- list(
  parkinson     = c("GSE145361", "GSE111629"),
  schizophrenia = c("GSE84727",  "GSE80417"),
  alzheimer     = c("GSE338167", "GSE144858"),
  depression    = c("GSE198904", "GSE273970"),
  ibd           = c("GSE87648",  "GSE112611")
)
LAB_EN <- c(parkinson = "Parkinson", schizophrenia = "Schizophrenia",
            alzheimer = "Alzheimer", depression = "Depression", ibd = "IBD")

# --------------------------------------------------------------
# 1. Estimate the composition effects cohort by cohort
# --------------------------------------------------------------
clr <- function(P) {
  # Zeros are replaced multiplicatively: set to 0.65 times the smallest non-zero value of the row, then renormalised;
  # this is standard for compositional data and avoids minus infinity when taking logs
  n0 <- base::sum(P <= 0, na.rm = TRUE)
  if (n0 > 0) {
    for (i in base::seq_len(base::nrow(P))) {
      z <- P[i, ] <= 0
      if (base::any(z)) {
        mn <- base::min(P[i, !z], na.rm = TRUE)
        P[i, z] <- 0.65 * mn
        P[i, ] <- P[i, ] / base::sum(P[i, ])
      }
    }
  }
  L <- base::log(P)
  list(x = L - base::rowMeans(L), n_zero = n0)
}

eff <- list(); smoke <- list(); stateD <- list()
for (d in names(DIS)) for (g in DIS[[d]]) {
  say("==== ", g, " (", d, ") ====")
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  cov <- read.csv(file.path(COV_DIR, base::paste0(g, "_covariates.csv")),
                  stringsAsFactors = FALSE)
  cov <- cov[base::match(D$gsm, cov$gsm), , drop = FALSE]
  ctc <- base::grep("^ct_", base::names(cov), value = TRUE)
  P <- base::as.matrix(cov[, ctc, drop = FALSE])
  base::colnames(P) <- base::sub("^ct_", "", ctc)
  say("  Cell types ", base::ncol(P), " types: ", base::paste(base::colnames(P), collapse = ", "))
  say("  Row sums of fractions: median ", sprintf("%.4f", stats::median(base::rowSums(P), na.rm = TRUE)),
      ", range ", sprintf("%.4f to %.4f", base::min(base::rowSums(P), na.rm = TRUE),
                         base::max(base::rowSums(P), na.rm = TRUE)))
  cc <- clr(P)
  say("  Zeros: ", cc$n_zero, ", replaced with 0.65 times the smallest non-zero value")
  Y <- cc$x

  # The design matrix is that of 39, with the cell fractions themselves removed (they are now the outcome)
  X <- D$design
  keep_col <- !base::grepl("^ct_", base::colnames(X))
  X <- X[, keep_col, drop = FALSE]
  jg <- base::which(base::colnames(X) == "groupcase")
  if (!length(jg)) stop(g, ": groupcase not found in the design matrix")
  X0 <- X[, -jg, drop = FALSE]
  q0 <- qr(X0)
  gv <- base::as.numeric(D$group == "case")
  gres <- qr.resid(q0, gv)
  Yres <- qr.resid(q0, Y)
  b <- base::as.numeric(base::crossprod(gres, Yres)) / base::sum(gres^2)
  base::names(b) <- base::colnames(Y)
  say("  Covariates: ", base::ncol(X0), " columns (removed ", base::sum(!keep_col), " cell fraction columns: ",
      base::paste(base::sub("^ct_", "", base::colnames(D$design)[!keep_col]), collapse = ", "), ")")
  # The outcome must never appear in the design as well. Cell fractions were removed; this is confirmed here once.
  if (base::any(base::grepl("^ct_", base::colnames(X0))))
    stop(g, ": cell fraction columns remain in the design matrix")
  say("  Composition effects (clr, cases minus controls): ")
  o <- base::order(-base::abs(b))
  say("    ", base::paste(sprintf("%s %+.4f", base::names(b)[o], b[o]), collapse = "  "))

  # Standardised case-control difference in the smoking score
  #
  # Two points need separate handling here, otherwise the numbers are meaningless:
  #   (1) The design matrix itself contains a smoking column. Residualising smoking on a design that contains it
  #       gives residuals identically zero (about 1e-16 in floating point), and the ratio is pure numerical noise.
  #       When smoking is the outcome, that column must therefore be removed from the design too.
  #   (2) The smoking score is missing for some samples (GSE112611 does not include the term under the first-round settings);
  #       a vector with NA passed to qr.resid raises an error. The decomposition is redone on the complete-case subset instead.
  sm <- cov$dnam_smoking
  ok_sm <- base::is.finite(sm)
  if (base::sum(ok_sm) < 50L) {
    dsd <- NA_real_
    say("  Smoking score: available ", base::sum(ok_sm), ", fewer than 50, not done for this cohort")
  } else {
    Xs <- X[, !base::grepl("^smoking", base::colnames(X)), drop = FALSE]
    js <- base::which(base::colnames(Xs) == "groupcase")
    Xs0 <- Xs[ok_sm, -js, drop = FALSE]
    qs <- qr(Xs0)
    gs <- qr.resid(qs, gv[ok_sm])
    ss <- qr.resid(qs, sm[ok_sm])
    sd_raw <- stats::sd(sm[ok_sm]); sd_res <- stats::sd(ss)
    if (sd_res < 1e-6 * base::max(sd_raw, 1e-12)) {
      dsd <- NA_real_
      say("  [FAIL] the residual SD of the smoking score is almost zero (", sprintf("%.3e", sd_res),
          "); it is completely explained by the design matrix and cannot be interpreted here")
    } else {
      bs <- base::sum(gs * ss) / base::sum(gs^2)
      dsd <- bs / sd_res
      say("  Smoking score: available ", base::sum(ok_sm), " samples; residual SD ",
          sprintf("%.4f", sd_res), " (originally ", sprintf("%.4f", sd_raw), ")")
      say("    standardised case-control difference ", sprintf("%+.3f", dsd), " SD")
    }
  }

  eff[[g]] <- b
  smoke[[length(smoke) + 1]] <- data.frame(gse = g, disease = d, smoke_d = dsd,
                                           stringsAsFactors = FALSE)
  stateD[[g]] <- list(q0 = q0, gv = gv, Yres = Yres, cells = base::colnames(Y))
  say("")
}

cells_all <- base::Reduce(base::intersect, base::lapply(eff, names))
say("Cell types common to the ten sets: ", length(cells_all), " types: ",
    base::paste(cells_all, collapse = ", "))
say("")

# --------------------------------------------------------------
# 2. The 45 comparisons
# --------------------------------------------------------------
pairs <- list()
for (d in names(DIS))
  pairs[[length(pairs) + 1]] <- list(kind = "within", dpair = d,
                                     a = DIS[[d]][1], b = DIS[[d]][2],
                                     id = base::paste0("PC-", LAB_EN[[d]]))
dn <- names(DIS)
for (i in 1:(length(dn) - 1)) for (j in (i + 1):length(dn))
  for (x in DIS[[dn[i]]]) for (y in DIS[[dn[j]]])
    pairs[[length(pairs) + 1]] <- list(
      kind = "cross", dpair = base::paste0(dn[i], "-", dn[j]), a = x, b = y,
      id = base::paste0(LAB_EN[[dn[i]]], "|", LAB_EN[[dn[j]]]))

perm_eff <- function(S, cells) {
  gr <- qr.resid(S$q0, base::sample(S$gv))
  v <- base::as.numeric(base::crossprod(gr, S$Yres)) / base::sum(gr^2)
  base::names(v) <- S$cells
  v[cells]
}

rows <- list()
for (P in pairs) {
  ea <- eff[[P$a]][cells_all]; eb <- eff[[P$b]][cells_all]
  r_p <- stats::cor(ea, eb)
  r_s <- stats::cor(ea, eb, method = "spearman")
  np <- numeric(N_PERM); ns <- numeric(N_PERM)
  Sb <- stateD[[P$b]]
  for (i in base::seq_len(N_PERM)) {
    v <- perm_eff(Sb, cells_all)
    np[i] <- stats::cor(ea, v)
    ns[i] <- stats::cor(ea, v, method = "spearman")
  }
  p_p <- (1 + base::sum(np >= r_p)) / (N_PERM + 1)
  p_s <- (1 + base::sum(ns >= r_s)) / (N_PERM + 1)
  z_p <- (r_p - base::mean(np)) / stats::sd(np)
  say(sprintf("  %-26s %-9s vs %-9s  Pearson %+.3f (Z=%+5.2f, P=%.4f)  Spearman %+.3f (P=%.4f)  null mean %+.3f",
              P$id, P$a, P$b, r_p, z_p, p_p, r_s, p_s, base::mean(np)))
  rows[[length(rows) + 1]] <- data.frame(
    id = P$id, kind = P$kind, dpair = P$dpair, a = P$a, b = P$b,
    r_pearson = r_p, z = z_p, p_pearson = p_p,
    r_spearman = r_s, p_spearman = p_s,
    null_mean = base::mean(np), null_sd = stats::sd(np), stringsAsFactors = FALSE)
}
tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "composition_all_pairs.csv"), row.names = FALSE)
write.csv(base::do.call(rbind, smoke), file.path(OUT_DIR, "smoking_by_cohort.csv"),
          row.names = FALSE)

say("")
say("=== Within-disease positive controls (composition axis) ===")
w <- tab[tab$kind == "within", ]
for (i in base::seq_len(base::nrow(w)))
  say(sprintf("  %-26s Pearson %+.3f  Z=%+5.2f  P=%.4f", w$id[i], w$r_pearson[i],
              w$z[i], w$p_pearson[i]))
say("  Reaching P < 0.05: ", base::sum(w$p_pearson < 0.05), " / ", base::nrow(w))
say("  Null means: median ", sprintf("%+.4f", stats::median(tab$null_mean)),
    ", range ", sprintf("%+.4f to %+.4f", base::min(tab$null_mean), base::max(tab$null_mean)),
    " (expected value 0; a large deviation indicates a problem with the permutation or the clr)")
say("  Note: with 12 components the null distribution of the correlation is bounded and non-normal; Z is a poor summary and the permutation P is decisive")

say("")
say("=== Cross-disease comparisons: the four pairings summarised by disease pair ===")
cr <- tab[tab$kind == "cross", ]
sumr <- list()
for (dp in base::unique(cr$dpair)) {
  s <- cr[cr$dpair == dp, ]
  say(sprintf("  %-28s Pearson of the four pairings: %s; significant %d/4",
              dp, base::paste(sprintf("%+.3f", s$r_pearson), collapse = " "),
              base::sum(s$p_pearson < 0.05)))
  sumr[[length(sumr) + 1]] <- data.frame(dpair = dp,
    n_sig = base::sum(s$p_pearson < 0.05),
    sign_consistent = base::all(base::sign(s$r_pearson) == base::sign(s$r_pearson[1])),
    r_median = stats::median(s$r_pearson), stringsAsFactors = FALSE)
}
ss <- base::do.call(rbind, sumr)
write.csv(ss, file.path(OUT_DIR, "composition_cross_summary.csv"), row.names = FALSE)
say("  Disease pairs with all four pairings significant and consistent in direction: ",
    base::sum(ss$n_sig == 4 & ss$sign_consistent), " / ", base::nrow(ss))

say("")
say("=== Direction common to the cohorts in the composition effects ===")
say("  Most cross-disease correlations are positive; how much of this agreement comes from one direction common to the cohorts needs to be established.")
say("  This section is an exploratory analysis added after the cross-disease correlations were seen to be mostly positive.")
EM <- base::do.call(rbind, base::lapply(base::unlist(DIS), function(g) eff[[g]][cells_all]))
base::rownames(EM) <- base::unlist(DIS)
dis_of <- stats::setNames(base::rep(names(DIS), each = 2), base::unlist(DIS))

# (0) Direction of each cell type across the ten cohorts: the most direct evidence, independent of any decomposition
sgn <- data.frame(cell = cells_all,
                  n_higher_in_cases = base::as.integer(base::colSums(EM > 0)),
                  n_lower_in_cases  = base::as.integer(base::colSums(EM < 0)),
                  stringsAsFactors = FALSE)
sgn <- sgn[base::order(-sgn$n_higher_in_cases), ]
say("  Per cell type: number of cohorts with cases above controls / below controls (of ", base::nrow(EM), " sets)")
for (i in base::seq_len(base::nrow(sgn)))
  say(sprintf("    %-8s %2d / %2d", sgn$cell[i], sgn$n_higher_in_cases[i], sgn$n_lower_in_cases[i]))
write.csv(sgn, file.path(OUT_DIR, "composition_sign_counts.csv"), row.names = FALSE)

# (1) Common direction: first singular vector without centring.
# Centring must not be applied: centring across cohorts would subtract the common component first, and the first principal component would then be
# "the direction of largest deviation of the cohorts from the common direction", which answers a different question.
# Before the decomposition the effect vector of each cohort is scaled to unit length. Pairwise concordance uses the Pearson correlation,
# which does not depend on vector length; without scaling the first direction is dominated by the cohort with the largest effects and gives "the direction
# of the cohort with the largest amplitude" rather than "the direction common to the cohorts". The unscaled result is listed at the end of this section for reference.
first_dir <- function(M) {
  s <- base::svd(M)
  v <- s$v[, 1]; base::names(v) <- base::colnames(M)
  list(v = v, share = s$d[1]^2 / base::sum(s$d^2))
}
EU <- EM / base::sqrt(base::rowSums(EM^2))
fd <- first_dir(EU)
u1 <- fd$v
if (u1[["Neu"]] < 0) u1 <- -u1          # The sign of the direction is arbitrary; neutrophils are set positive
share <- fd$share
o <- base::order(-base::abs(u1))
say("  [Main analysis: effect vectors scaled to unit length]")
say("  First direction: share of the total sum of squares ", sprintf("%.1f%%", 100 * share))
say("  Loadings (sorted by absolute value): ")
say("    ", base::paste(sprintf("%s %+.3f", base::names(u1)[o], u1[o]), collapse = "  "))
proj <- base::as.numeric(EU %*% u1); base::names(proj) <- base::rownames(EU)
say("  Projection of each cohort on this direction (projection of a unit vector, i.e. the cosine of its angle with the direction): ")
for (d in names(DIS))
  say(sprintf("    %-14s %s", LAB_EN[[d]],
              base::paste(sprintf("%s %+.3f", DIS[[d]], proj[DIS[[d]]]), collapse = "  ")))
say("  Cohorts with positive projection: ", base::sum(proj > 0), " / ", length(proj))

# (2) How much pairwise correlation remains after removing this direction, computed in two ways:
#   in-sample:  the direction is estimated from all ten sets and removed from the same ten. Its estimate uses the two sets it is removed from,
#           so part of their noise is removed with it, which overestimates the role of the direction (upper bound).
#   leave-out:  for each pair of cohorts the direction is estimated only from the cohorts of the other disorders and removed from the two sets.
#           The estimate is independent of the pair but rests on fewer cohorts and is less precise, which underestimates the role (lower bound).
EUr <- EU - base::outer(proj, u1)
rem <- base::do.call(rbind, base::lapply(pairs, function(P) {
  keep <- base::names(dis_of)[!dis_of %in% c(dis_of[[P$a]], dis_of[[P$b]])]
  v <- first_dir(EU[keep, , drop = FALSE])$v
  ra <- EU[P$a, ] - base::sum(EU[P$a, ] * v) * v
  rb <- EU[P$b, ] - base::sum(EU[P$b, ] * v) * v
  data.frame(id = P$id, kind = P$kind, dpair = P$dpair, a = P$a, b = P$b,
             n_ref_cohorts = length(keep),
             r_before = stats::cor(EU[P$a, ], EU[P$b, ]),
             r_after_insample = stats::cor(EUr[P$a, ], EUr[P$b, ]),
             r_after_leaveout = stats::cor(ra, rb), stringsAsFactors = FALSE)
}))
# Scaling does not change the Pearson correlation: the values before removal must equal the results of the 45 pairs above one by one
dev_r <- base::max(base::abs(rem$r_before -
  tab$r_pearson[base::match(base::paste(rem$a, rem$b), base::paste(tab$a, tab$b))]))
if (!base::is.finite(dev_r) || dev_r > 1e-8)
  stop("Correlations before removal differ from the results of the 45 pairs; largest difference ", sprintf("%.2e", dev_r))
md <- function(kind, col) stats::median(rem[[col]][rem$kind == kind])
say("  Median correlation of the 40 cross-disease pairs: before removal ", sprintf("%+.3f", md("cross", "r_before")),
    "; after in-sample removal ", sprintf("%+.3f", md("cross", "r_after_insample")),
    "; after leave-out removal ", sprintf("%+.3f", md("cross", "r_after_leaveout")))
say("  Median correlation of the  5 within-disease pairs: before removal ", sprintf("%+.3f", md("within", "r_before")),
    "; after in-sample removal ", sprintf("%+.3f", md("within", "r_after_insample")),
    "; after leave-out removal ", sprintf("%+.3f", md("within", "r_after_leaveout")))
say("  Number of the 40 cross-disease pairs with positive correlation: before removal ", base::sum(rem$r_before[rem$kind == "cross"] > 0),
    "; after in-sample removal ", base::sum(rem$r_after_insample[rem$kind == "cross"] > 0),
    "; after leave-out removal ", base::sum(rem$r_after_leaveout[rem$kind == "cross"] > 0))
say("  Interpretation: the in-sample and leave-out computations give the upper and lower bound of the role of this direction.")

write.csv(data.frame(cell = base::names(u1), loading = base::as.numeric(u1),
                     stringsAsFactors = FALSE),
          file.path(OUT_DIR, "composition_shared_axis_loadings.csv"), row.names = FALSE)
write.csv(data.frame(gse = base::rownames(EU), projection = base::as.numeric(proj),
                     stringsAsFactors = FALSE),
          file.path(OUT_DIR, "composition_shared_axis_projection.csv"), row.names = FALSE)
write.csv(rem, file.path(OUT_DIR, "composition_axis_removal_pairs.csv"), row.names = FALSE)

# (3) For reference: results without scaling, and which cohorts contribute the sum of squares of the first direction
fr <- first_dir(EM); ur <- fr$v; if (ur[["Neu"]] < 0) ur <- -ur
pr <- base::as.numeric(EM %*% ur); base::names(pr) <- base::rownames(EM)
EMr <- EM - base::outer(pr, ur)
raw_after <- base::sapply(pairs, function(P) stats::cor(EMr[P$a, ], EMr[P$b, ]))
raw_kind <- base::sapply(pairs, function(P) P$kind)
say("  [For reference: unscaled] first direction: share of the total sum of squares ", sprintf("%.1f%%", 100 * fr$share),
    "; share contributed by each cohort: ")
say("    ", base::paste(sprintf("%s %.0f%%", base::names(pr), 100 * pr^2 / base::sum(pr^2)), collapse = "  "))
say("    median correlation after in-sample removal: cross-disease ", sprintf("%+.3f", stats::median(raw_after[raw_kind == "cross"])),
    "; within-disease ", sprintf("%+.3f", stats::median(raw_after[raw_kind == "within"])))
say("    cosine of the angle with the scaled direction ", sprintf("%.3f", base::sum(ur * u1)))

write.csv(data.frame(
  item = c("scaling", "share_first_direction", "n_positive_projection",
           "cross_before", "cross_after_insample", "cross_after_leaveout",
           "within_before", "within_after_insample", "within_after_leaveout",
           "cross_n_positive_before", "cross_n_positive_after_leaveout",
           "raw_share_first_direction", "raw_max_cohort_share_of_first_direction",
           "raw_cross_after_insample", "raw_within_after_insample"),
  value = c("unit length", share, base::sum(proj > 0),
            md("cross", "r_before"), md("cross", "r_after_insample"), md("cross", "r_after_leaveout"),
            md("within", "r_before"), md("within", "r_after_insample"), md("within", "r_after_leaveout"),
            base::sum(rem$r_before[rem$kind == "cross"] > 0),
            base::sum(rem$r_after_leaveout[rem$kind == "cross"] > 0),
            fr$share, base::max(pr^2 / base::sum(pr^2)),
            stats::median(raw_after[raw_kind == "cross"]), stats::median(raw_after[raw_kind == "within"])),
  stringsAsFactors = FALSE),
  file.path(OUT_DIR, "composition_axis_summary.csv"), row.names = FALSE)

say("")
say("=== Case-control difference in the smoking score (standardised) ===")
sm <- base::do.call(rbind, smoke)
for (d in names(DIS)) {
  s <- sm[sm$disease == d, ]
  say(sprintf("  %-14s %s", LAB_EN[[d]],
              base::paste(sprintf("%s %+.3f", s$gse, s$smoke_d), collapse = "  ")))
}
say("  Disorders with the same sign in both cohorts: ",
    base::sum(base::sapply(names(DIS), function(d) {
      s <- sm[sm$disease == d, ]
      base::all(base::is.finite(s$smoke_d)) && s$smoke_d[1] * s$smoke_d[2] > 0
    })), " / 5")
say("  Note 1: absolute scores are not comparable across cohorts (different numbers of usable sites shift them as a whole); ")
say("        only the standardised difference within each cohort is used here.")
say("  Note 2: the model of this section is not adjusted for cell fractions - smoking itself affects cell composition, and removing composition")
say("        as a covariate would remove part of the smoking effect.")
say("  Note 3: the smoking score of the paediatric cohort should not be read as true smoking exposure; the sites of the score also")
say("        respond to inflammation and other factors, and large differences must be understood accordingly.")

say("")
say("=== Plotting ===")
hm <- eff[base::unlist(DIS)]
M <- base::do.call(rbind, base::lapply(hm, function(v) v[cells_all]))
dfm <- data.frame(
  gse = base::rep(base::rownames(M), times = base::ncol(M)),
  cell = base::rep(cells_all, each = base::nrow(M)),
  eff = base::as.numeric(M), stringsAsFactors = FALSE)
dfm$disease <- base::unlist(base::lapply(dfm$gse, function(g)
  LAB_EN[[base::names(DIS)[base::sapply(DIS, function(z) g %in% z)]]]))
dfm$gse <- factor(dfm$gse, levels = base::unlist(DIS))
dfm$cell <- factor(dfm$cell, levels = cells_all)
lim <- base::max(base::abs(dfm$eff), na.rm = TRUE)
p <- ggplot(dfm, aes(cell, gse, fill = eff)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = PAL[1], mid = "white", high = PAL[2], midpoint = 0,
                       limits = c(-lim, lim), name = "clr effect") +
  facet_grid(disease ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Case-control shift in leukocyte composition, by cohort",
       subtitle = "Centred log-ratio effects; two independent cohorts per disease",
       x = NULL, y = NULL) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.placement = "outside")
save_fig(p, "fig43_composition_heatmap", width = 9, height = 7)

tab$panelkind <- factor(ifelse(tab$kind == "within",
                               "Within-disease (positive control)", "Cross-disease"),
                        levels = c("Within-disease (positive control)", "Cross-disease"))
p <- ggplot(tab, aes(stats::reorder(id, r_pearson), r_pearson, fill = panelkind)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Within-disease (positive control)" = PAL[3],
                               "Cross-disease" = PAL[1]), name = NULL) +
  labs(title = "Concordance of leukocyte composition shifts across cohort pairs",
       subtitle = "Correlation of centred log-ratio effects over shared cell types",
       x = NULL, y = "Pearson correlation of clr effects") +
  theme_paper() + theme(axis.text.y = element_text(size = 6))
save_fig(p, "fig44_composition_all_pairs", width = 9, height = 9)
log_stop()
