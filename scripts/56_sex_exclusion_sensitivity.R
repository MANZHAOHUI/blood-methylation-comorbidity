# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 56: sensitivity analysis - exclusion of the sex-discordant samples of GSE87648
#
# Background: step 42 found six samples of GSE87648 whose X-chromosome methylation lies clearly in the
# cluster of the opposite sex (five controls, one case). The primary analysis retained them. This step
# excludes them and recomputes every cohort pair that involves GSE87648:
#   within-disease   1 pair   inflammatory bowel disease, GSE87648 vs GSE112611 (10000 permutations, as in step 45)
#   cross-disease    8 pairs  GSE87648 against the two cohorts of each of the other four disorders
#                             (1000 permutations, as in step 39)
#
# Model: the design matrix of the primary model of GSE87648 (step 54) with the rows of the excluded samples
# removed. All covariates, including the methylation principal components of the primary model, are kept as
# they are; columns that become collinear after the exclusion (for example a batch level left without
# samples) are dropped. The models of the other nine cohorts are unchanged.
#
# Statistic and null: as in steps 39 and 45. Spearman correlation of the case-control effect vectors over
# the common probes; null by permuting case-control labels in the cohort that was permuted in the primary
# analysis (roles are read from the table of step 46). The group coefficient of every probe is computed as
# M g~ / (g~' g~), where g~ is the group vector residualised on the other covariates. This equals the
# Frisch-Waugh-Lovell coefficient used in steps 39 and 45 (residualising M as well does not change the
# inner product) and the limma coefficient; the script checks this against the stored results and stops
# if they disagree.
#
# One difference in computation, none in the test: the six pairs in which GSE87648 is the permuted cohort
# share the same 1000 label permutations (each permuted effect vector is compared with all six fixed
# vectors). Every pair still has its own valid permutation null of 1000 draws.
#
# This step changes no existing result. Output: output/56_sex_exclusion/
#   sex_exclusion_pairs.csv           the nine pairs, primary values next to the values after exclusion
#   sex_exclusion_nulls.csv           permutation nulls
#   sex_exclusion_cross_summary.csv   significant pairings per disease pair, before and after
#   sex_exclusion_key_numbers.csv     the numbers quoted in the manuscript
#
# Environment: Windows + RStudio. About one hour. Can run unattended
# Location: <project root>/scripts/56_sex_exclusion_sensitivity.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "56_sex_exclusion_sensitivity"

SET_DIR  <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR  <- file.path(DATA_ROOT, "ewas_final")
EXC_FILE <- file.path(PROJ_ROOT, "output", "42_sexmismatch", "GSE87648_exclude_candidates.csv")
REF_FILE <- file.path(PROJ_ROOT, "output", "46_limits", "all_pairs_with_limits.csv")
OUT_DIR  <- file.path(PROJ_ROOT, "output", "56_sex_exclusion")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

GSE           <- "GSE87648"
N_EXPECTED    <- 6L        # number of sex-discordant samples reported in the manuscript
N_PERM_WITHIN <- 10000L    # as in step 45
N_PERM_CROSS  <- 1000L     # as in step 39
BETA_CLIP     <- 0.001
POWER_TARGET  <- 0.80
SEED          <- 20261123L
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
check_against_stored <- function(g, e_new, probes, E) {
  e_old <- E$effect[base::match(probes, E$probe)]
  if (base::anyNA(e_old)) stop(g, ": probes of the design are missing from the stored results")
  r <- stats::cor(e_new, e_old)
  say("  ", g, ": recomputed coefficients against the stored limma results: r = ", sprintf("%.8f", r),
      ", largest absolute difference ", sprintf("%.2e", base::max(base::abs(e_new - e_old))))
  if (!base::is.finite(r) || r < 0.999)
    stop(g, ": recomputed coefficients do not reproduce the stored results (r = ", r, "); the design or the analysis set has changed, stopping")
  invisible(r)
}
null_stats <- function(rho, pr) {
  nm <- base::mean(pr); ns <- stats::sd(pr)
  pv <- (1 + base::sum(pr >= rho)) / (length(pr) + 1)
  data.frame(null_mean = nm, null_sd = ns, z = (rho - nm) / ns, p = pv,
             p_se = base::sqrt(pv * (1 - pv) / length(pr)),
             rho_crit = nm + Z95 * ns, rho_power = nm + (Z95 + ZPW) * ns)
}
mins_since <- function(t0) sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins")))

# --------------------------------------------------------------
say("=== 1. Inputs ===")
for (f in c(EXC_FILE, REF_FILE)) if (!file.exists(f)) stop("File not found: ", f)
REF <- read.csv(REF_FILE, stringsAsFactors = FALSE)
need <- c("id", "kind", "dpair", "fixed", "permuted", "n_probe", "rho", "p_rho", "rho_crit", "rho_power")
if (!base::all(need %in% base::names(REF)))
  stop("The table of step 46 lacks columns: ", base::paste(base::setdiff(need, base::names(REF)), collapse = ", "))
if (base::nrow(REF) != 45L) stop("The table of step 46 should have 45 pairs; found ", base::nrow(REF))
PAIRS <- REF[REF$fixed == GSE | REF$permuted == GSE, , drop = FALSE]
if (base::nrow(PAIRS) != 9L || base::sum(PAIRS$kind == "within") != 1L)
  stop("Expected 1 within-disease and 8 cross-disease pairs involving ", GSE, "; found ",
       base::sum(PAIRS$kind == "within"), " and ", base::sum(PAIRS$kind == "cross"))
say("Pairs involving ", GSE, ": ", base::nrow(PAIRS), " (", GSE, " permuted in ",
    base::sum(PAIRS$permuted == GSE), ", fixed in ", base::sum(PAIRS$fixed == GSE), ")")

EX <- read.csv(EXC_FILE, stringsAsFactors = FALSE)
if (!("gsm" %in% base::names(EX))) stop("The exclusion list of step 42 has no column 'gsm'")
ex <- base::unique(base::as.character(EX$gsm))
D <- read_rds(file.path(RES_DIR, base::paste0(GSE, "_design.rds")))
E <- read_rds(file.path(RES_DIR, base::paste0(GSE, "_ewas_final.rds")))
if (base::nrow(D$design) != length(D$gsm) || length(D$group) != length(D$gsm))
  stop(GSE, ": design, sample list and group vector differ in length")
if (!base::all(ex %in% D$gsm))
  stop("Samples of the exclusion list that are not in the current model of ", GSE, ": ",
       base::paste(base::setdiff(ex, D$gsm), collapse = ", "))
if (length(ex) != N_EXPECTED)
  stop("The exclusion list holds ", length(ex), " samples; the manuscript reports ", N_EXPECTED,
       ". Step 42 has changed; check before running this step")
gex <- base::as.character(D$group[base::match(ex, D$gsm)])
say("Samples to exclude: ", length(ex), " (controls ", base::sum(gex == "control"), ", cases ",
    base::sum(gex == "case"), "): ", base::paste(ex, collapse = ", "))

# --------------------------------------------------------------
say("")
say("=== 2. ", GSE, ": primary model and model after exclusion ===")
t0 <- Sys.time()
M <- load_M(GSE, D)
say("  M matrix ", base::nrow(M), " x ", base::ncol(M), ", loaded in ", mins_since(t0), " min")
g_full <- base::as.numeric(D$group == "case")
e_full <- coef_group(M, null_qr(D$design), g_full)
check_against_stored(GSE, e_full, D$probes, E)

keep <- !(D$gsm %in% ex)
des_r <- D$design[keep, , drop = FALSE]
qq <- base::qr(des_r)
dropped <- character(0)
if (qq$rank < base::ncol(des_r)) {
  keepc <- base::sort(qq$pivot[base::seq_len(qq$rank)])
  dropped <- base::colnames(des_r)[base::setdiff(base::seq_len(base::ncol(des_r)), keepc)]
  des_r <- des_r[, keepc, drop = FALSE]
}
if (!("groupcase" %in% base::colnames(des_r))) stop("The group column became collinear after the exclusion; the model cannot be fitted")
M <- M[, keep, drop = FALSE]; invisible(gc())
g_r <- g_full[keep]
q_r <- null_qr(des_r)
e_red <- coef_group(M, q_r, g_r)
say("  Samples ", length(keep), " -> ", base::sum(keep), "; controls ", base::sum(g_full == 0), " -> ",
    base::sum(g_r == 0), ", cases ", base::sum(g_full == 1), " -> ", base::sum(g_r == 1))
say("  Parameters ", base::ncol(D$design), " -> ", base::ncol(des_r),
    if (length(dropped)) base::paste0("; columns dropped as collinear: ", base::paste(dropped, collapse = ", ")) else "; no column dropped")
r_eff <- stats::cor(e_full, e_red)
say("  Pearson correlation of the effect vector before and after exclusion: ", sprintf("%.5f", r_eff),
    " (", length(e_red), " probes)")

rows <- list(); nulls <- list()
add_result <- function(P, n_probe, n_perm, rho, dir, pr) {
  if (n_probe != P$n_probe)
    stop(P$id, " (", P$fixed, " vs ", P$permuted, "): ", n_probe, " common probes, primary analysis ", P$n_probe,
         "; the probe sets have changed, stopping")
  st <- null_stats(rho, pr)
  say(sprintf("  %-22s %-9s vs %-9s probes %6d  primary rho=%+.4f P=%.4f | after exclusion rho=%+.4f Z=%+5.2f P=%.4f",
              P$id, P$fixed, P$permuted, n_probe, P$rho, P$p_rho, rho, st$z, st$p))
  rows[[length(rows) + 1]] <<- cbind(
    data.frame(id = P$id, kind = P$kind, dpair = P$dpair, fixed = P$fixed, permuted = P$permuted,
               n_probe = n_probe, n_perm = n_perm, rho_primary = P$rho, p_primary = P$p_rho,
               rho_crit_primary = P$rho_crit, rho_power_primary = P$rho_power,
               rho = rho, dir = dir, stringsAsFactors = FALSE), st)
  nulls[[length(nulls) + 1]] <<- data.frame(id = P$id, fixed = P$fixed, permuted = P$permuted,
                                             rho = pr, obs = rho, stringsAsFactors = FALSE)
}

# --------------------------------------------------------------
say("")
say("=== 3. Pairs in which ", GSE, " is the permuted cohort ===")
PP <- PAIRS[PAIRS$permuted == GSE, , drop = FALSE]
if (base::any(PP$kind != "cross")) stop("A within-disease pair has ", GSE, " as the permuted cohort; this differs from the primary analysis")
prep <- list()
for (k in base::seq_len(base::nrow(PP))) {
  A <- read_rds(file.path(RES_DIR, base::paste0(PP$fixed[k], "_ewas_final.rds")))
  cm <- base::intersect(A$probe, D$probes)
  ia <- base::match(cm, A$probe); ib <- base::match(cm, D$probes)
  ea <- A$effect[ia]; eb <- e_red[ib]
  prep[[k]] <- list(ib = ib, ra = base::rank(ea), n = length(cm),
                    rho = stats::cor(ea, eb, method = "spearman"),
                    dir = base::mean(base::sign(ea) == base::sign(eb)))
  rm(A)
}
set.seed(SEED)
pr <- base::matrix(NA_real_, N_PERM_CROSS, base::nrow(PP))
t0 <- Sys.time()
for (i in base::seq_len(N_PERM_CROSS)) {
  cf <- coef_group(M, q_r, base::sample(g_r))
  for (k in base::seq_len(base::nrow(PP)))
    pr[i, k] <- stats::cor(prep[[k]]$ra, base::rank(cf[prep[[k]]$ib]))
  if (i %% 200L == 0L) say("    permutation ", i, " / ", N_PERM_CROSS, ", elapsed ", mins_since(t0), " min")
}
if (base::anyNA(pr)) stop("Missing values in the permutation null")
for (k in base::seq_len(base::nrow(PP)))
  add_result(PP[k, ], prep[[k]]$n, N_PERM_CROSS, prep[[k]]$rho, prep[[k]]$dir, pr[, k])
rm(M, pr, prep); invisible(gc())

# --------------------------------------------------------------
say("")
say("=== 4. Pairs in which ", GSE, " is the fixed cohort ===")
PF <- PAIRS[PAIRS$fixed == GSE, , drop = FALSE]
for (k in base::seq_len(base::nrow(PF))) {
  pg <- PF$permuted[k]
  n_perm <- if (PF$kind[k] == "within") N_PERM_WITHIN else N_PERM_CROSS
  say("  ---- ", pg, " (", PF$id[k], ", ", n_perm, " permutations) ----")
  Dp <- read_rds(file.path(RES_DIR, base::paste0(pg, "_design.rds")))
  Ep <- read_rds(file.path(RES_DIR, base::paste0(pg, "_ewas_final.rds")))
  t0 <- Sys.time()
  Mp <- load_M(pg, Dp)
  say("  M matrix ", base::nrow(Mp), " x ", base::ncol(Mp), ", loaded in ", mins_since(t0), " min")
  qp <- null_qr(Dp$design)
  gp <- base::as.numeric(Dp$group == "case")
  eb_all <- coef_group(Mp, qp, gp)
  check_against_stored(pg, eb_all, Dp$probes, Ep)
  cm <- base::intersect(D$probes, Dp$probes)
  ia <- base::match(cm, D$probes); ib <- base::match(cm, Dp$probes)
  ea <- e_red[ia]; eb <- eb_all[ib]
  rho <- stats::cor(ea, eb, method = "spearman")
  dir <- base::mean(base::sign(ea) == base::sign(eb))
  ra <- base::rank(ea)
  set.seed(SEED + k)
  prk <- numeric(n_perm)
  t0 <- Sys.time()
  for (i in base::seq_len(n_perm)) {
    prk[i] <- stats::cor(ra, base::rank(coef_group(Mp, qp, base::sample(gp))[ib]))
    if (i %% 1000L == 0L) say("    permutation ", i, " / ", n_perm, ", elapsed ", mins_since(t0), " min")
  }
  if (base::anyNA(prk)) stop("Missing values in the permutation null")
  add_result(PF[k, ], length(cm), n_perm, rho, dir, prk)
  rm(Mp, Dp, Ep, eb_all); invisible(gc())
}

# --------------------------------------------------------------
say("")
say("=== 5. Summary ===")
tab <- base::do.call(rbind, rows)
tab$delta_rho <- tab$rho - tab$rho_primary
tab$sig_primary <- tab$p_primary < 0.05
tab$sig <- tab$p < 0.05
write.csv(tab, file.path(OUT_DIR, "sex_exclusion_pairs.csv"), row.names = FALSE)
write.csv(base::do.call(rbind, nulls), file.path(OUT_DIR, "sex_exclusion_nulls.csv"), row.names = FALSE)

w <- tab[tab$kind == "within", ]
say(sprintf("  Within-disease pair (inflammatory bowel disease): primary rho=%+.4f P=%.4f; after exclusion rho=%+.4f Z=%+.2f P=%.4f (Monte Carlo SE %.4f)",
            w$rho_primary, w$p_primary, w$rho, w$z, w$p, w$p_se))
say(sprintf("    critical rho %.4f, %.0f%% power requires rho >= %.4f (primary %.4f and %.4f)",
            w$rho_crit, 100 * POWER_TARGET, w$rho_power, w$rho_crit_primary, w$rho_power_primary))
cr <- tab[tab$kind == "cross", ]
say("  Cross-disease pairs involving ", GSE, ": P < 0.05 in ", base::sum(cr$sig_primary), " of ", base::nrow(cr),
    " in the primary analysis and in ", base::sum(cr$sig), " of ", base::nrow(cr), " after exclusion")
say(sprintf("    largest absolute change in rho %.4f; correlation between primary and new rho over the %d pairs r = %.4f",
            base::max(base::abs(cr$delta_rho)), base::nrow(cr), stats::cor(cr$rho_primary, cr$rho)))
chg <- tab[tab$sig != tab$sig_primary, , drop = FALSE]
if (base::nrow(chg)) {
  say("  Pairs whose status at P < 0.05 changed: ", base::nrow(chg))
  for (i in base::seq_len(base::nrow(chg)))
    say(sprintf("    %-22s %s vs %s: P %.4f -> %.4f", chg$id[i], chg$fixed[i], chg$permuted[i], chg$p_primary[i], chg$p[i]))
} else say("  No pair changed its status at P < 0.05")

# The forty cross-disease pairs with the eight recomputed pairs substituted
CR <- REF[REF$kind == "cross", , drop = FALSE]
key_ref <- base::paste(CR$fixed, CR$permuted)
key_new <- base::paste(cr$fixed, cr$permuted)
if (!base::all(key_new %in% key_ref)) stop("Recomputed pairs do not match the table of step 46")
CR$p_new <- CR$p_rho; CR$rho_new <- CR$rho; CR$power_new <- CR$rho_power
j <- base::match(key_new, key_ref)
CR$p_new[j] <- cr$p; CR$rho_new[j] <- cr$rho; CR$power_new[j] <- cr$rho_power
sumr <- list()
for (dp in base::unique(CR$dpair)) {
  s <- CR[CR$dpair == dp, , drop = FALSE]
  sumr[[length(sumr) + 1]] <- data.frame(
    dpair = dp, n_pairings = base::nrow(s),
    n_sig_primary = base::sum(s$p_rho < 0.05), n_sig_after = base::sum(s$p_new < 0.05),
    involves_gse87648 = base::any(s$fixed == GSE | s$permuted == GSE), stringsAsFactors = FALSE)
}
ss <- base::do.call(rbind, sumr)
write.csv(ss, file.path(OUT_DIR, "sex_exclusion_cross_summary.csv"), row.names = FALSE)
say("")
say("  All forty cross-disease pairs, with the eight recomputed pairs substituted:")
for (i in base::which(ss$involves_gse87648))
  say(sprintf("    %-28s significant pairings %d/4 -> %d/4", ss$dpair[i], ss$n_sig_primary[i], ss$n_sig_after[i]))
say("    nominally significant: ", base::sum(CR$p_rho < 0.05), " -> ", base::sum(CR$p_new < 0.05), " of ", base::nrow(CR),
    "; largest number per disease pair: ", base::max(ss$n_sig_primary), " -> ", base::max(ss$n_sig_after),
    "; disease pairs significant in all four pairings: ", base::sum(ss$n_sig_primary == 4), " -> ", base::sum(ss$n_sig_after == 4))
say(sprintf("    cross-disease pairs with a detection limit below the within-disease rho of inflammatory bowel disease: %d -> %d of %d",
            base::sum(CR$rho_power < w$rho_primary), base::sum(CR$power_new < w$rho), base::nrow(CR)))

key <- data.frame(
  item = c("n_excluded", "excluded_controls", "excluded_cases", "n_samples_after", "n_parameters_after",
           "effect_vector_cor_before_vs_after",
           "ibd_within_rho_primary", "ibd_within_rho", "ibd_within_z", "ibd_within_p", "ibd_within_rho_crit", "ibd_within_rho_power",
           "cross8_n_sig_primary", "cross8_n_sig", "cross8_max_abs_delta_rho",
           "cross40_n_sig_primary", "cross40_n_sig", "cross40_max_per_dpair_primary", "cross40_max_per_dpair",
           "dpairs_all_four_primary", "dpairs_all_four",
           "cross40_limits_below_ibd_primary", "cross40_limits_below_ibd"),
  value = c(length(ex), base::sum(gex == "control"), base::sum(gex == "case"), base::sum(keep), base::ncol(des_r),
            r_eff,
            w$rho_primary, w$rho, w$z, w$p, w$rho_crit, w$rho_power,
            base::sum(cr$sig_primary), base::sum(cr$sig), base::max(base::abs(cr$delta_rho)),
            base::sum(CR$p_rho < 0.05), base::sum(CR$p_new < 0.05), base::max(ss$n_sig_primary), base::max(ss$n_sig_after),
            base::sum(ss$n_sig_primary == 4), base::sum(ss$n_sig_after == 4),
            base::sum(CR$rho_power < w$rho_primary), base::sum(CR$power_new < w$rho)),
  stringsAsFactors = FALSE)
write.csv(key, file.path(OUT_DIR, "sex_exclusion_key_numbers.csv"), row.names = FALSE)
say("")
say("Tables written to ", OUT_DIR)
log_stop()
