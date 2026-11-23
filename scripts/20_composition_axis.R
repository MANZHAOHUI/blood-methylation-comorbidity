# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 20: the second axis - sharing of leukocyte composition
#
# The first axis asks "which sites change" (concordance of CpG effect sizes, done in 17/18).
# This one asks "which cell types change". Both use the same data but measure
# different things: the former is the methylation change within cell types (cell fractions
# are removed in the model), the latter is exactly the part that was removed. The two axes do not contain each other,
# and agreement or disagreement is interpretable either way.
#
#   both axes        -> sharing at both the composition level and within cells, the strongest evidence
#   this axis only   -> what is shared is a change in peripheral immune cell composition, a common downstream
#               inflammatory process rather than a common intracellular epigenetic programme
#   first axis only  -> sharing is within cell types and unrelated to immune cell fractions
#   neither axis     -> negative, and a negative supported by two independent measurements
#
# The procedure parallels the first axis exactly: effects estimated within cohorts, effect sizes compared across cohorts,
# null distribution from permuted case-control labels, the same positive and negative controls.
# There are only twelve dimensions instead of four hundred thousand, so the correlation coefficient itself is noisy;
# the sign concordance (how many of twelve types agree in direction) is more stable, and both are reported.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/20_composition_axis.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "20_composition_axis"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "20_composition")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

set.seed(20260915)
N_PERM <- 10000L          # With only twelve dimensions permutations are cheap, so more are run
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

ROLE <- c(GSE145361 = "PD discovery", GSE111629 = "PD replication",
          GSE72774 = "PD sensitivity", GSE87648 = "IBD main",
          GSE112611 = "IBD paediatric", GSE152027 = "Negative control")
ORDER <- names(ROLE)
AGE_SOURCE <- c(GSE145361 = "dnam", GSE152027 = "dnam", GSE111629 = "chrono",
                GSE72774 = "chrono", GSE87648 = "chrono", GSE112611 = "chrono")
USE_SMOKING <- c(GSE145361 = TRUE, GSE152027 = TRUE, GSE111629 = TRUE,
                 GSE72774 = TRUE, GSE87648 = TRUE, GSE112611 = FALSE)

PAIRS <- list(
  list(id = "PC1", kind = "positive control", a = "GSE145361", b = "GSE111629"),
  list(id = "PC2", kind = "positive control", a = "GSE87648",  b = "GSE112611"),
  list(id = "M1",  kind = "primary criterion",   a = "GSE145361", b = "GSE87648"),
  list(id = "M2",  kind = "secondary",     a = "GSE111629", b = "GSE87648"),
  list(id = "NC1", kind = "negative control", a = "GSE145361", b = "GSE152027"),
  list(id = "NC2", kind = "negative control", a = "GSE87648",  b = "GSE152027")
)

slide_of <- function(g, gsm) {
  mf <- file.path(MAP_DIR, paste0(g, "_sample_map.csv"))
  if (!file.exists(mf)) return(rep(NA_character_, length(gsm)))
  m <- read.csv(mf, stringsAsFactors = FALSE, colClasses = "character")
  bc <- rep(NA_character_, nrow(m))
  hit <- regexpr(BARCODE, m$col_name) > 0
  bc[hit] <- regmatches(m$col_name, regexpr(BARCODE, m$col_name))
  sub("_R[0-9]{2}C[0-9]{2}$", "", bc)[match(gsm, m$gsm)]
}

# Fractions are modelled after a logit transform: fractions are bounded between 0 and 1, and direct linear regression
# has unequal variance near the bounds; after the logit they are closer to normal
lg <- function(p) {
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
  log(p / (1 - p))
}

# --------------------------------------------------------------
say("=== Composition effects of each cohort ===")
eff <- list(); dat <- list()
for (g in names(ROLE)) {
  f <- file.path(SET_DIR, paste0(g, "_analysis_filtered.rds"))
  cf <- file.path(COV_DIR, paste0(g, "_covariates.csv"))
  if (!file.exists(f) || !file.exists(cf)) { say("  ", g, ": file missing, skipped"); next }
  obj <- readRDS(f); ph <- obj$pheno
  cov <- read.csv(cf, stringsAsFactors = FALSE)
  cov <- cov[match(colnames(obj$beta), cov$gsm), , drop = FALSE]

  d <- data.frame(group = factor(as.character(ph$group),
                                 levels = c("control", "case")))
  sx <- grep("^(sex|gender)$", names(ph), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(sx)) d$sex <- factor(toupper(substr(trimws(ph[[sx]]), 1, 1)))
  if (AGE_SOURCE[g] == "dnam") d$age <- cov$dnam_age_horvath else {
    ag <- grep("^age$|ageat", names(ph), ignore.case = TRUE, value = TRUE)[1]
    d$age <- suppressWarnings(as.numeric(ph[[ag]]))
  }
  if (isTRUE(USE_SMOKING[g])) d$smoking <- cov$dnam_smoking
  if ("scan_date" %in% names(ph)) d$batch <- factor(as.character(ph$scan_date))
  else {
    sl <- slide_of(g, colnames(obj$beta))
    if (sum(!is.na(sl)) > 0.95 * length(sl) &&
        length(unique(sl)) <= 0.25 * length(sl)) d$batch <- factor(sl)
  }

  ct <- grep("^ct_", names(cov), value = TRUE)
  Y <- as.matrix(cov[, ct, drop = FALSE])
  colnames(Y) <- sub("^ct_", "", ct)
  ok <- stats::complete.cases(d) & stats::complete.cases(Y)
  d <- d[ok, , drop = FALSE]; Y <- Y[ok, , drop = FALSE]
  tb <- table(d$group)
  say("  ", g, " (", ROLE[g], "): controls ", tb[["control"]], ", cases ", tb[["case"]],
      ", cell types ", ncol(Y))

  X <- stats::model.matrix(~ ., data = d)
  jg <- which(colnames(X) == "groupcase")
  b <- numeric(ncol(Y)); se <- numeric(ncol(Y)); p <- numeric(ncol(Y))
  for (k in seq_len(ncol(Y))) {
    fit <- stats::lm.fit(X, lg(Y[, k]))
    r <- fit$residuals
    s2 <- sum(r^2) / fit$df.residual
    XtXi <- chol2inv(fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE])
    b[k] <- fit$coefficients[jg]
    se[k] <- sqrt(s2 * XtXi[jg, jg])
    p[k] <- 2 * stats::pt(abs(b[k] / se[k]), fit$df.residual, lower.tail = FALSE)
  }
  names(b) <- colnames(Y)
  eff[[g]] <- b
  dat[[g]] <- list(Y = Y, X = X, jg = jg, group = d$group)
  sig <- which(stats::p.adjust(p, "BH") < 0.05)
  say("     cell types with FDR<0.05: ", length(sig),
      if (length(sig)) paste0(": ", paste(sprintf("%s(%+.2f)", colnames(Y)[sig],
                                                  b[sig]), collapse = ", ")) else "")
}
say("")

E <- do.call(rbind, eff)
write.csv(data.frame(gse = rownames(E), E), 
          file.path(OUT_DIR, "composition_effects.csv"), row.names = FALSE)

# --------------------------------------------------------------
recompute <- function(g, perm = FALSE) {
  D <- dat[[g]]
  X <- D$X
  if (perm) X[, D$jg] <- as.numeric(sample(D$group) == "case")
  out <- numeric(ncol(D$Y))
  for (k in seq_len(ncol(D$Y))) {
    fit <- stats::lm.fit(X, lg(D$Y[, k]))
    out[k] <- fit$coefficients[D$jg]
  }
  names(out) <- colnames(D$Y)
  out
}

say("=== Pairwise concordance ===")
res <- list(); nulls <- list()
for (P in PAIRS) {
  if (is.null(eff[[P$a]]) || is.null(eff[[P$b]])) next
  ct <- intersect(names(eff[[P$a]]), names(eff[[P$b]]))
  a <- eff[[P$a]][ct]; b <- eff[[P$b]][ct]
  rho <- stats::cor(a, b, method = "spearman")
  r <- stats::cor(a, b)
  dir <- mean(sign(a) == sign(b))
  say("==== ", P$id, "  ", P$kind, ": ", P$a, " vs ", P$b, " ====")
  say("  Cell types ", length(ct), "; Spearman rho = ", sprintf("%.3f", rho),
      ", Pearson r = ", sprintf("%.3f", r),
      ", sign concordance ", sum(sign(a) == sign(b)), " / ", length(ct))

  pr <- numeric(N_PERM); pd <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    bp <- recompute(P$b, perm = TRUE)[ct]
    pr[i] <- stats::cor(a, bp, method = "spearman")
    pd[i] <- mean(sign(a) == sign(bp))
  }
  p_rho <- (1 + sum(pr >= rho)) / (N_PERM + 1)
  p_dir <- (1 + sum(pd >= dir)) / (N_PERM + 1)
  z <- (rho - mean(pr)) / stats::sd(pr)
  say("  Null rho: mean ", sprintf("%.3f", mean(pr)),
      ", SD ", sprintf("%.3f", stats::sd(pr)))
  say("  Empirical P: rho ", format.pval(p_rho, digits = 3),
      "; sign concordance ", format.pval(p_dir, digits = 3),
      "; Z = ", sprintf("%.2f", z))
  res[[length(res) + 1]] <- data.frame(
    pair = P$id, kind = P$kind, a = P$a, b = P$b, n_ct = length(ct),
    rho = rho, r = r, dir = dir, null_mean = mean(pr),
    null_sd = stats::sd(pr), z = z, p_rho = p_rho, p_dir = p_dir,
    stringsAsFactors = FALSE)
  nulls[[P$id]] <- data.frame(pair = P$id, kind = P$kind, rho = pr, obs = rho,
                              stringsAsFactors = FALSE)
  say("")
}
rr <- do.call(rbind, res)
write.csv(rr, file.path(OUT_DIR, "composition_concordance.csv"), row.names = FALSE)
say("=== Summary ===")
print(rr[, c("pair", "kind", "rho", "dir", "z", "p_rho")])

# --------------------------------------------------------------
say("")
say("=== Plotting ===")
long <- do.call(rbind, lapply(rownames(E), function(g)
  data.frame(gse = g, cell = colnames(E), effect = as.numeric(E[g, ]),
             stringsAsFactors = FALSE)))
long$panel <- factor(paste0(long$gse, "\n", ROLE[long$gse]),
                     levels = paste0(ORDER, "\n", ROLE[ORDER]))
long$cell <- factor(long$cell, levels = colnames(E))
lim <- max(abs(long$effect), na.rm = TRUE)
p <- ggplot(long, aes(cell, panel, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = PAL[1], mid = "grey96", high = PAL[2],
                       midpoint = 0, limits = c(-lim, lim),
                       name = "Case minus control\n(logit scale)") +
  labs(title = "Case-control shift in leukocyte composition",
       subtitle = "Adjusted for sex, age, smoking score and batch; twelve cell types per cohort",
       x = NULL, y = NULL) +
  theme_paper() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
save_fig(p, "fig25_composition_effects", width = 9, height = 4.4)

lab <- paste0(rr$pair, "  ", rr$kind, "\n", rr$a, " vs ", rr$b)
rr$panel <- factor(lab, levels = lab)
nd <- do.call(rbind, nulls)
nd$panel <- factor(rr$panel[match(nd$pair, rr$pair)], levels = levels(rr$panel))
p <- ggplot(nd, aes(rho)) +
  geom_histogram(bins = 45, fill = "grey75", colour = "white", linewidth = 0.15) +
  geom_vline(aes(xintercept = obs), colour = PAL[2], linewidth = 0.7) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  labs(title = "Composition concordance against the permutation null",
       subtitle = paste0(format(N_PERM, big.mark = ","),
                         " permutations of the case-control label within cohort"),
       x = "Spearman correlation of composition effects across twelve cell types",
       y = "Permutations") +
  theme_paper()
save_fig(p, "fig26_composition_null", width = 10, height = 6)

log_stop()
