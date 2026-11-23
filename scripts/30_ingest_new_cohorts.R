# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 30: reading, verification and phenotype extraction of the five additional cohorts
#
# This step does one thing only: confirm that the downloaded matrices are usable and extract group, sex and age.
# No analysis is done. If a matrix has a problem, the earlier it is found the better.
#
# Each data set carries expected counts against which the script checks the parsed result. A mismatch
# is written to the log explicitly as FAIL and the script does not continue silently - in 19b a smoking field
# could pass as the group field precisely because the criterion could not fail.
#
# Known features of each cohort (from the dossier of 29, not guessed):
#   GSE84727   disease_status: 1 = control 433, 2 = case 414; the sentrixids field carries barcodes
#   GSE80417   disease status: 1 = control 322, 2 = case 353
#   GSE198904  diagnosis is mixed: MDD 361 / HC 100 / OBSERVEMDD0001 191 /
#              Excluded 2 / Not depressed 1. Only MDD and HC are taken, the rest removed
#   GSE273970  disease state: normal 147 / MDD 136. Age is given in ten-year groups
#              (such as 41-50), not as a number, so DNAm age is used throughout
#   GSE338167  disease state: AD 313 / control 308. EPIC v2 platform;
#              probe names may carry a suffix (cg00000029_TC21) that must be stripped before intersecting
#              with 450K. There are also ethnicity and cohort fields;
#              the cohort spans five population groups, which must be covariates later
#   GSE144858  disease state: AD 93 / MCI 111 / control 96. MCI is removed,
#              only AD and controls are kept
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/30_ingest_new_cohorts.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "30_ingest_new_cohorts"

RAW_DIR <- file.path(DATA_ROOT, "raw")          # series matrix files are placed here
RDS_DIR <- file.path(DATA_ROOT, "rds")
MAP_DIR <- file.path(PROJ_ROOT, "output", "30_ingest")
for (d in c(RDS_DIR, MAP_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
suppressPackageStartupMessages(library(data.table))

N_PEEK <- 5000L        # Read this many rows first for verification; the full file is read once everything is confirmed

SPEC <- list(
  GSE84727 = list(
    disease = "schizophrenia", platform = "450K",
    group_tag = "disease_status",
    map = c("1" = "control", "2" = "case"),
    expect = c(control = 433, case = 414),
    sex_tag = "Sex", age_tag = "age", age_source = "chrono"),
  GSE80417 = list(
    disease = "schizophrenia", platform = "450K",
    group_tag = "disease status",
    map = c("1" = "control", "2" = "case"),
    expect = c(control = 322, case = 353),
    sex_tag = "Sex", age_tag = "age", age_source = "chrono"),
  # Cohort 2 only: Cohort 1 consists of controls alone (HC 32, MDD 0) and is completely confounded with group;
  # it was identified and removed in 03c. The corresponding data matrix is DHRC (header matches 858/864).
  GSE198904 = list(
    disease = "depression", platform = "EPIC",
    group_tag = "diagnosis",
    map = c("HC" = "control", "MDD" = "case"),
    expect = c(control = 68, case = 361),
    restrict_tag = "cohort", restrict_keep = "Cohort 2",
    sex_tag = "gender", age_tag = "age", age_source = "chrono"),
  GSE273970 = list(
    disease = "depression", platform = "EPIC",
    group_tag = "disease state",
    map = c("normal" = "control", "major depressive disorder" = "case"),
    expect = c(control = 147, case = 136),
    sex_tag = "Sex", age_tag = NA_character_, age_source = "dnam"),
  GSE338167 = list(
    disease = "alzheimer", platform = "EPICv2",
    group_tag = "disease state",
    map = c("Healthy control" = "control", "Alzheimer's disease" = "case"),
    expect = c(control = 308, case = 313),
    sex_tag = "Sex", age_tag = "age", age_source = "chrono"),
  GSE144858 = list(
    disease = "alzheimer", platform = "450K",
    group_tag = "disease state",
    map = c("control" = "control", "Alzheimer's disease" = "case"),
    expect = c(control = 96, case = 93),
    sex_tag = "Sex", age_tag = "age", age_source = "chrono")
)

# Headers always come from the series matrix (phenotypes exist only there); the data table may be in
# the series matrix or in a supplementary matrix, and both are searched
find_header <- function(g) {
  list.files(RAW_DIR, pattern = paste0("^", g, "(-GPL[0-9]+)?_series_matrix\\.txt(\\.gz)?$"),
             full.names = TRUE, recursive = TRUE)
}
find_matrix <- function(g) {
  f <- list.files(RAW_DIR, pattern = paste0("^", g, "_.*\\.(txt|csv)\\.gz$"),
                  full.names = TRUE, recursive = TRUE)
  f <- f[file.info(f)$size > 5 * 1024^2]     # Files below 5 MB are bare headers, not data tables
  f
}
read_head <- function(f, n = 400L) {
  con <- if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}
split_line <- function(x) trimws(gsub('^"|"$', "", strsplit(x, "\t", fixed = TRUE)[[1]]))

ok_all <- TRUE
summ <- list()
for (g in names(SPEC)) {
  S <- SPEC[[g]]
  say("==== ", g, " (", S$disease, ", ", S$platform, ") ====")
  hs <- find_header(g)
  if (!length(hs)) { say("  [FAIL] series matrix not found; phenotypes cannot be extracted"); ok_all <- FALSE; next }
  fs <- find_matrix(g)
  if (!length(fs)) { say("  [FAIL] no matrix containing data found (larger than 5 MB)"); ok_all <- FALSE; next }
  f <- fs[which.max(file.info(fs)$size)]
  say("  Header ", basename(hs[1]), "; data ", basename(f),
      " (", sprintf("%.2f GB", file.info(f)$size / 1024^3), ")")

  hd <- read_head(hs[1], 400L)
  ib <- which(grepl("!series_matrix_table_begin", hd, fixed = TRUE))
  if (!length(ib)) { say("  [FAIL] data start marker not found in the header"); ok_all <- FALSE; next }
  gsm <- split_line(grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  say("  Samples ", length(gsm))

  # Group, sex, age
  pull <- function(tag) {
    for (one in grep("^!Sample_characteristics_ch1", hd, value = TRUE)) {
      v <- split_line(one)[-1]
      vv <- v[nzchar(v) & grepl(":", v, fixed = TRUE)]
      if (!length(vv)) next
      t0 <- trimws(sub(":.*$", "", vv[1]))
      if (tolower(t0) == tolower(tag)) return(trimws(sub("^[^:]*:", "", v)))
    }
    NULL
  }
  raw_grp <- pull(S$group_tag)
  if (is.null(raw_grp)) { say("  [FAIL] group field not found: '", S$group_tag, "'"); ok_all <- FALSE; next }
  grp <- unname(S$map[raw_grp])
  tb_raw <- sort(table(raw_grp), decreasing = TRUE)
  say("  Raw values of the group field: ",
      paste(sprintf("%s=%d", names(tb_raw), as.integer(tb_raw)), collapse = "; "))
  keep <- !is.na(grp)
  # Sub-study restriction: some data sets combine several sub-studies; only the specified one is kept
  if (!is.null(S$restrict_tag)) {
    rv <- pull(S$restrict_tag)
    if (is.null(rv)) {
      say("  [FAIL] restriction field not found: '", S$restrict_tag, "'"); ok_all <- FALSE; next
    }
    n_before <- sum(keep)
    keep <- keep & rv %in% S$restrict_keep
    say("  Restricted by '", S$restrict_tag, " = ", paste(S$restrict_keep, collapse = ", "),
        "': ", n_before, " reduced to ", sum(keep), " samples")
  }
  say("  Included: ", sum(keep), ", removed ", sum(!keep), " (neither case nor control, or not in the specified sub-study)")
  tb <- table(factor(grp[keep], levels = c("control", "case")))
  say("  Parsed: controls ", tb[["control"]], ", cases ", tb[["case"]])
  hit <- tb[["control"]] == S$expect[["control"]] && tb[["case"]] == S$expect[["case"]]
  say("  Expected: controls ", S$expect[["control"]], ", cases ", S$expect[["case"]],
      "  ", if (hit) "pass" else "[FAIL]")
  if (!hit) ok_all <- FALSE

  if (g == "GSE198904") {
    coh <- pull("cohort")
    if (!is.null(coh)) {
      ct <- table(coh[keep], grp[keep])
      say("  [KEY CHECK] cross-tabulation of sub-study and group (both sub-studies must contain cases and controls, ")
      say("             otherwise case-control status and sub-study are completely confounded and the data set is unusable): ")
      for (r in rownames(ct))
        say("    ", r, ": ", paste(sprintf("%s=%d", colnames(ct), as.integer(ct[r, ])),
                                   collapse = "; "))
      if (any(ct == 0)) {
        say("  [FAIL] a sub-study lacks one of the groups; group and sub-study are confounded")
        ok_all <- FALSE
      }
    } else say("  [NOTE] cohort field not found; confounding check not possible")
  }

  sx <- if (!is.na(S$sex_tag)) pull(S$sex_tag) else NULL
  if (is.null(sx)) say("  [NOTE] sex field not found: '", S$sex_tag, "'") else
    say("  Sex: ", paste(sprintf("%s=%d", names(table(sx)), as.integer(table(sx))),
                          collapse = "; "))
  ag <- if (!is.na(S$age_tag)) pull(S$age_tag) else NULL
  if (S$age_source == "dnam") {
    say("  Age: this cohort only has ten-year groups; DNAm age is used according to the plan, not extracted here")
  } else if (is.null(ag)) {
    say("  [NOTE] age field not found: '", S$age_tag, "'")
  } else {
    an <- suppressWarnings(as.numeric(ag))
    say("  Age: convertible to numeric ", sum(is.finite(an)), " / ", length(an),
        ", range ", sprintf("%.1f to %.1f", min(an, na.rm = TRUE), max(an, na.rm = TRUE)))
  }

  # Check of the data block: number of columns, whether values look like beta, whether columns are paired
  dt <- data.table::fread(f, skip = ib[1], nrows = N_PEEK, header = TRUE,
                          showProgress = FALSE)
  nm <- names(dt)
  say("  Data block: ", nrow(dt), " rows (trial read), ", ncol(dt), " columns")
  ncol_per <- (ncol(dt) - 1L) / length(gsm)
  say("  Columns per sample: ", sprintf("%.2f", ncol_per),
      if (abs(ncol_per - 2) < 0.01) " (paired, beta and detection P value side by side)"
      else if (abs(ncol_per - 1) < 0.01) " (single)" else "  [NOTE] not an integer; manual inspection needed")
  num <- suppressWarnings(as.matrix(dt[, 2:min(ncol(dt), 50L)]))
  num <- num[is.finite(num)]
  inrange <- mean(num >= 0 & num <= 1)
  say("  Proportion of values between 0 and 1: ", sprintf("%.3f", inrange),
      if (inrange > 0.95) " (looks like beta values)" else "  [FAIL] does not look like beta values; the supplementary matrix must be used instead")
  if (inrange <= 0.95) ok_all <- FALSE

  pid <- as.character(dt[[1]])
  npid <- sum(grepl("^cg", pid))
  suffixed <- sum(grepl("^cg[0-9]+_", pid))
  say("  Probe names: starting with cg ", npid, " / ", length(pid),
      "; with suffix ", suffixed,
      if (suffixed > 0) "  -> EPIC v2 naming; suffixes must be stripped before intersecting" else "")

  # The keep column must be written: the result of the sub-study restriction exists only in memory,
  # and without saving it downstream scripts would read all samples without restriction
  restrict_val <- if (!is.null(S$restrict_tag)) {
    rv <- pull(S$restrict_tag); if (is.null(rv)) NA else rv
  } else NA
  write.csv(data.frame(gsm = gsm, group_raw = raw_grp, group = grp,
                       sex = if (is.null(sx)) NA else sx,
                       age = if (is.null(ag)) NA else ag,
                       restrict = restrict_val, keep = keep,
                       stringsAsFactors = FALSE),
            file.path(MAP_DIR, paste0(g, "_pheno.csv")), row.names = FALSE)
  summ[[length(summ) + 1]] <- data.frame(
    gse = g, disease = S$disease, platform = S$platform, n_total = length(gsm),
    n_control = tb[["control"]], n_case = tb[["case"]],
    cols_per_sample = ncol_per, beta_like = inrange,
    probe_suffix = suffixed > 0, pass = hit, stringsAsFactors = FALSE)
  say("")
}

if (length(summ)) {
  tab <- do.call(rbind, summ)
  write.csv(tab, file.path(MAP_DIR, "ingest_summary.csv"), row.names = FALSE)
  say("=== Summary ===")
  print(tab)
}
say("")
say(if (ok_all) "All checks passed; analysis sets can be built and the common probes recomputed"
    else "[SOME CHECKS FAILED] resolve the problems flagged above before proceeding")
log_stop()
