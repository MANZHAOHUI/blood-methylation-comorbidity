# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 32: phenotype cleaning and inclusion gate for the six additional cohorts
#
# All count gates of 30 passed, but they only checked group sizes, not age and missingness.
# Two points in the log must be handled first:
#
#   GSE80417 shows an age range of 18.0 to 891.0 - an age of 891 years cannot exist.
#   Unless removed, wrong values enter all models as outliers and bias the coefficient of the age term.
#
#   GSE84727 has only 665 / 847 ages convertible to numeric; 182 are missing (21%).
#   Complete-case modelling drops these individuals, and if missingness is related to group, the remaining
#   samples are no longer the original cases and controls, which creates selection bias. It is therefore not enough to look at
#   how many are missing; who is missing matters too.
#
# This step does four things, each of which can fail:
#   (1) flag implausible ages and set them to missing (18 to 110 for adult cohorts)
#   (2) test whether missing age is related to group (chi-square); if so, report it
#   (3) compute the final usable sample size of each cohort under complete-case analysis
#   (4) gate: cases >= 50 and controls >= 30, otherwise FAIL
#
# A false alarm of 30 is also recorded: GSE198904 reported "columns per sample 1.32" because of a wrong denominator -
# the matrix has 865 columns (1 probe name column + 864 data columns), 864 / 432 = 2, a clean paired layout,
# whereas the script used the 655 samples of the header as denominator. The data are fine; this step corrects the figure.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/32_pheno_qc_gate.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "32_pheno_qc_gate"

IN_DIR  <- file.path(PROJ_ROOT, "output", "30_ingest")
OUT_DIR <- file.path(PROJ_ROOT, "output", "32_pheno")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

AGE_MIN <- 18; AGE_MAX <- 110      # Plausible range for adult cohorts
MIN_CASE <- 50L; MIN_CTRL <- 30L   # Minimum group sizes, the standard used when the panel was defined

# The chronological age of GSE84727 is unusable: 154 of 414 cases lack it (37%) but only 28 of 433
# controls (6%), test of independence p < 2e-16. Complete-case modelling would drop a third of the
# cases, and the remainder would no longer represent the original case group. Decision D4 already allows "DNAm age when
# too many chronological ages are missing", as applied to GSE145361 and GSE152027, and the same is done here -
# DNAm age is computed from the methylation data and exists for every sample, which removes the missingness problem.
SPEC <- list(
  GSE84727  = list(disease = "schizophrenia", age_source = "dnam"),
  GSE80417  = list(disease = "schizophrenia", age_source = "chrono"),
  GSE198904 = list(disease = "depression",    age_source = "chrono"),
  GSE273970 = list(disease = "depression",    age_source = "dnam"),
  GSE338167 = list(disease = "alzheimer",     age_source = "chrono"),
  GSE144858 = list(disease = "alzheimer",     age_source = "chrono")
)

ok_all <- TRUE; summ <- list()
for (g in names(SPEC)) {
  S <- SPEC[[g]]
  f <- file.path(IN_DIR, paste0(g, "_pheno.csv"))
  say("==== ", g, " (", S$disease, ") ====")
  if (!file.exists(f)) { say("  [FAIL] not found: ", f); ok_all <- FALSE; next }
  ph <- read.csv(f, stringsAsFactors = FALSE)
  n0 <- nrow(ph)
  ph <- ph[!is.na(ph$group) & nzchar(ph$group), , drop = FALSE]
  # Respect the sub-study restriction written by 30
  if ("keep" %in% names(ph)) {
    n1 <- nrow(ph)
    ph <- ph[isTRUE(TRUE) & ph$keep %in% c(TRUE, "TRUE"), , drop = FALSE]
    if (nrow(ph) < n1) say("  Restricted by sub-study: ", n1, " reduced to ", nrow(ph), " samples")
  }
  say("  Samples with non-missing group and within the specified sub-study: ", nrow(ph), " (originally ", n0, ")")

  # ---- (1) age cleaning ----
  if (S$age_source == "dnam") {
    say("  Age source: DNAm age (this cohort only has ten-year groups); chronological age is not checked in this step")
    ph$age_clean <- NA_real_
    age_ok <- rep(TRUE, nrow(ph))
  } else {
    a <- suppressWarnings(as.numeric(ph$age))
    n_na0 <- sum(is.na(a))
    bad <- which(!is.na(a) & (a < AGE_MIN | a > AGE_MAX))
    if (length(bad)) {
      say("  [NOTE] ages outside ", AGE_MIN, " to ", AGE_MAX, ": ", length(bad),
          " samples, values: ", paste(sort(unique(a[bad])), collapse = ", "), ", set to missing")
      a[bad] <- NA_real_
    }
    ph$age_clean <- a
    age_ok <- !is.na(a)
    say("  Age: originally missing ", n_na0, ", missing after cleaning ", sum(is.na(a)),
        " ( ", sprintf("%.1f%%", 100 * mean(is.na(a))), ")",
        if (any(age_ok)) paste0(", valid range ",
          sprintf("%.1f to %.1f", min(a, na.rm = TRUE), max(a, na.rm = TRUE))) else "")
  }

  # ---- (2) is missingness related to group? ----
  if (S$age_source != "dnam" && any(!age_ok)) {
    tb <- table(missing_age = !age_ok, group = ph$group)
    say("  Cross-tabulation of missing age and group: ")
    for (r in rownames(tb))
      say("    ", if (r == "TRUE") "missing" else "non-missing", ": ",
          paste(sprintf("%s=%d", colnames(tb), as.integer(tb[r, ])), collapse = "; "))
    if (nrow(tb) == 2L && ncol(tb) == 2L && all(dim(tb) == 2)) {
      p <- suppressWarnings(stats::chisq.test(tb)$p.value)
      say("  Test of independence: p = ", format.pval(p, digits = 3),
          if (p < 0.05) "  [NOTE] missingness is related to group; complete-case modelling introduces selection bias, to be stated in the limitations"
          else "  Missingness is unrelated to group; complete-case modelling is acceptable")
    }
  }

  # ---- (3) complete cases ----
  sex_ok <- !is.na(ph$sex) & nzchar(as.character(ph$sex))
  comp <- age_ok & sex_ok
  tb2 <- table(factor(ph$group[comp], levels = c("control", "case")))
  say("  Sex missing: ", sum(!sex_ok), " samples")
  say("  Complete cases: controls ", tb2[["control"]], ", cases ", tb2[["case"]],
      ", total ", sum(comp), " (originally ", nrow(ph), ")")

  # ---- (4) gate ----
  pass <- tb2[["case"]] >= MIN_CASE && tb2[["control"]] >= MIN_CTRL
  say("  Gate (cases >= ", MIN_CASE, " and controls >= ", MIN_CTRL, "): ",
      if (pass) "pass" else "[FAIL]")
  if (!pass) ok_all <- FALSE

  ph$complete <- comp
  write.csv(ph, file.path(OUT_DIR, paste0(g, "_pheno_clean.csv")), row.names = FALSE)
  summ[[length(summ) + 1]] <- data.frame(
    gse = g, disease = S$disease, n_group = nrow(ph),
    n_age_missing = sum(!age_ok), n_sex_missing = sum(!sex_ok),
    n_complete = sum(comp), n_control = tb2[["control"]], n_case = tb2[["case"]],
    pass = pass, stringsAsFactors = FALSE)
  say("")
}

tab <- do.call(rbind, summ)
write.csv(tab, file.path(OUT_DIR, "pheno_qc_summary.csv"), row.names = FALSE)
say("=== Summary ===")
print(tab)

say("")
say("=== Final sample sizes of the five-disorder panel ===")
OLD <- data.frame(rbind(
  c("parkinson",     "GSE145361", 787, 851),
  c("parkinson",     "GSE111629", 237, 335),
  c("ibd",           "GSE87648",  176, 204),
  c("ibd",           "GSE112611",  74, 164)
), stringsAsFactors = FALSE)
names(OLD) <- c("disease", "gse", "n_control", "n_case")
OLD$n_control <- as.integer(OLD$n_control); OLD$n_case <- as.integer(OLD$n_case)
NEW <- tab[, c("disease", "gse", "n_control", "n_case")]
ALL <- rbind(OLD, NEW)
for (d in unique(ALL$disease)) {
  s <- ALL[ALL$disease == d, ]
  say(sprintf("  %-14s %s  total %d",
              d, paste(sprintf("%s %d/%d", s$gse, s$n_control, s$n_case), collapse = "  "),
              sum(s$n_control + s$n_case)))
}
say("  Total: ", sum(ALL$n_control + ALL$n_case), " individuals, ",
    length(unique(ALL$disease)), " disorders, ",
    choose(length(unique(ALL$disease)), 2), " cross-disease comparisons")
write.csv(ALL, file.path(OUT_DIR, "panel_final_n.csv"), row.names = FALSE)

say("")
say(if (ok_all) "All checks passed; the common probe set can be recomputed and the analysis sets built"
    else "[SOME CHECKS FAILED] resolve the problems flagged above first")
log_stop()
