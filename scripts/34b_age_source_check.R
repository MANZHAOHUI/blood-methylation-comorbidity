# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 34b: decide which age each cohort uses and check whether the proxy is credible
#
# 34 computed all three covariates, but two points about age must be settled first:
#
#   (1) In the two Alzheimer's disease cohorts DNAm age and chronological age correlate at only r = 0.65,
#       clearly below the 0.82 to 0.94 of the first round. Their chronological ages are complete, however,
#       and under D4 DNAm age only replaces chronological age when the latter is unusable - so these two sets
#       use chronological age anyway; the 0.65 does not affect the analysis but must be on record.
#
#   (2) GSE84727 is the only adult cohort forced to use DNAm age (21.5% of chronological ages
#       are missing, with missingness biased towards cases). The set that needs the proxy most
#       is the hardest one in which to verify it. It still has 665 samples with chronological age,
#       on which the two can be compared - the strongest verification available, and it must be done.
#
# This step only reads existing csv files, recomputes nothing and finishes within a minute.
# The output is a table of the age source of the ten cohorts, to be cited by the study plan and the Methods.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/34b_age_source_check.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "34b_age_source_check"

COV_DIR   <- file.path(DATA_ROOT, "covariates")
SET_DIR   <- file.path(DATA_ROOT, "analysis_sets")
ING_DIR   <- file.path(PROJ_ROOT, "output", "30_ingest")
PHENO_DIR <- file.path(PROJ_ROOT, "output", "32_pheno")
OUT_DIR   <- file.path(PROJ_ROOT, "output", "34_covariates")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)

AGE_MIN <- 18; AGE_MAX <- 110

# The four sets of the first round keep their established source; the six of this round are decided here
PANEL <- list(
  list(gse = "GSE145361", disease = "Parkinson's disease",     round = 1, source = "dnam",   why = "no chronological age field in this cohort"),
  list(gse = "GSE111629", disease = "Parkinson's disease",     round = 1, source = "chrono", why = "chronological age complete"),
  list(gse = "GSE87648",  disease = "inflammatory bowel disease",   round = 1, source = "chrono", why = "chronological age complete"),
  list(gse = "GSE112611", disease = "inflammatory bowel disease",   round = 1, source = "chrono", why = "chronological age complete"),
  list(gse = "GSE84727",  disease = "schizophrenia",   round = 2, source = "dnam",   why = "chronological age missing in 21.5%, biased towards cases"),
  list(gse = "GSE80417",  disease = "schizophrenia",   round = 2, source = "chrono", why = "chronological age missing in 5.8%, unrelated to group"),
  list(gse = "GSE198904", disease = "depression",       round = 2, source = "chrono", why = "chronological age complete"),
  list(gse = "GSE273970", disease = "depression",       round = 2, source = "dnam",   why = "age only in ten-year groups"),
  list(gse = "GSE338167", disease = "Alzheimer's disease", round = 2, source = "chrono", why = "chronological age complete"),
  list(gse = "GSE144858", disease = "Alzheimer's disease", round = 2, source = "chrono", why = "chronological age complete")
)

say("=== 1. GSE84727: validating DNAm age in the 665 samples with chronological age ===")
f_ing <- file.path(ING_DIR, "GSE84727_pheno.csv")
f_cov <- file.path(COV_DIR, "GSE84727_covariates.csv")
if (!file.exists(f_ing) || !file.exists(f_cov)) {
  say("  [FAIL] file missing; validation not possible")
} else {
  ing <- read.csv(f_ing, stringsAsFactors = FALSE)
  cov <- read.csv(f_cov, stringsAsFactors = FALSE)
  a <- base::suppressWarnings(base::as.numeric(ing$age))
  a[!base::is.na(a) & (a < AGE_MIN | a > AGE_MAX)] <- NA_real_
  m <- base::match(cov$gsm, ing$gsm)
  chrono <- a[m]; dnam <- cov$dnam_age_horvath
  ok <- base::is.finite(chrono) & base::is.finite(dnam)
  say("  Comparable samples: ", base::sum(ok), " / ", base::length(dnam))
  if (base::sum(ok) > 30) {
    r <- stats::cor(chrono[ok], dnam[ok])
    say("  Correlation r = ", sprintf("%.3f", r),
        "; mean difference ", sprintf("%+.1f", base::mean(dnam[ok] - chrono[ok])), " years",
        "; median absolute error ", sprintf("%.1f", stats::median(base::abs(dnam[ok] - chrono[ok]))), " years")
    grp <- ing$group[m]
    for (g in c("control", "case")) {
      s <- ok & grp == g
      if (base::sum(s) > 30)
        say("    ", g, ": n = ", base::sum(s), ", r = ",
            sprintf("%.3f", stats::cor(chrono[s], dnam[s])),
            ", mean difference ", sprintf("%+.1f", base::mean(dnam[s] - chrono[s])), " years")
    }
    say("  Interpretation: the proxy is considered usable at r of 0.85 or above; ",
        if (r >= 0.85) "this case passes" else "this case is low and must be stated in the limitations")
    # Whether the proxy is systematically biased between the two groups - more important than the overall correlation
    s1 <- ok & grp == "control"; s2 <- ok & grp == "case"
    if (base::sum(s1) > 30 && base::sum(s2) > 30) {
      d1 <- dnam[s1] - chrono[s1]; d2 <- dnam[s2] - chrono[s2]
      tt <- stats::t.test(d1, d2)
      say("  Difference in estimation bias between the two groups: ", sprintf("%+.2f", base::mean(d2) - base::mean(d1)),
          " years, p = ", format.pval(tt$p.value, digits = 3),
          if (tt$p.value < 0.05) "  [NOTE] the proxy is systematically biased between the groups and partly absorbs the group effect"
          else "  No difference in bias between the groups; the proxy introduces no group-related bias")
    }
  }
}
say("")

say("=== 2. Age source of the ten cohorts ===")
rows <- list()
for (P in PANEL) {
  fc <- file.path(COV_DIR, base::paste0(P$gse, "_covariates.csv"))
  n_dnam <- NA_integer_; med <- NA_real_
  if (file.exists(fc)) {
    cv <- read.csv(fc, stringsAsFactors = FALSE)
    if ("dnam_age_horvath" %in% base::names(cv)) {
      n_dnam <- base::sum(!base::is.na(cv$dnam_age_horvath))
      med <- stats::median(cv$dnam_age_horvath, na.rm = TRUE)
    }
  }
  say(sprintf("  %-10s %-12s round %d  uses %-6s  %s",
              P$gse, P$disease, P$round,
              if (P$source == "dnam") "DNAm age" else "chronological age", P$why))
  rows[[length(rows) + 1]] <- data.frame(
    gse = P$gse, disease = P$disease, round = P$round,
    age_source = P$source, reason = P$why,
    n_dnam_age = n_dnam, median_dnam_age = med, stringsAsFactors = FALSE)
}
tab <- base::do.call(rbind, rows)
write.csv(tab, file.path(OUT_DIR, "age_source_table.csv"), row.names = FALSE)
say("")
say("  Using DNAm age: ", base::sum(tab$age_source == "dnam"), " sets, ",
    "using chronological age: ", base::sum(tab$age_source == "chrono"), "")
say("  Age source table written, to be cited by the study plan and the Methods")

say("")
say("=== 3. Two notes for the Methods ===")
say("  (1) In the two Alzheimer's disease cohorts DNAm age and chronological age correlate at only about 0.65, ")
say("      below the 0.82 to 0.94 of the other cohorts. Both use chronological age, so the analysis is unaffected; ")
say("      the low value may reflect compression of the Horvath clock at old age, plus the imputation of")
say("      1015 and 586 reference sites during normalisation.")
say("  (2) The number of usable sites for the smoking score differs between cohorts (163 to 187), which shifts the scores as a whole; ")
say("      they are not comparable across cohorts. As a within-cohort covariate the shift is absorbed by the intercept and does not matter.")
log_stop()
