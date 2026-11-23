# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 04: read the cohort matrices, check the format, take the common probes
#
# Design principle: only one data set in memory at a time. Each set is converted to a numeric matrix right after reading,
# saved as .rds and removed from memory. Peak memory is then set by the largest set,
# not by the sum of the six.
#
# The value range of every set is checked: beta values should lie between 0 and 1.
# If a file actually holds methylated/unmethylated signal intensities rather than beta values,
# this step reports it immediately instead of letting it slip through to modelling.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/04_load_and_common_probes.R
# ==============================================================

DATA_ROOT <- "D:/blood_methylation_data"
PROJ_ROOT <- "C:/blood_methylation_project"

if (!dir.exists(DATA_ROOT)) stop("Data directory does not exist: ", DATA_ROOT)
if (!dir.exists(PROJ_ROOT)) stop("Project directory does not exist: ", PROJ_ROOT)

RDS_DIR <- file.path(DATA_ROOT, "rds")
OUT_DIR <- file.path(PROJ_ROOT, "output", "04_load")
for (d in c(RDS_DIR, OUT_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
}

message("Data root: ", DATA_ROOT)
message("Current contents: ")
for (x in list.files(DATA_ROOT)) message("  ", x)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is missing; please run first: install.packages(\"data.table\")")
}
library(data.table)
setDTthreads(0L)

# --------------------------------------------------------------
# Files to read. role defines the role of the data set in the design.
# --------------------------------------------------------------
COHORTS <- list(
  list(gse = "GSE87648",  role = "IBD_main",    dis = "IBD",
       file = "GSE87648/GSE87648_processedMethWb.txt.gz"),
  list(gse = "GSE72774",  role = "PD_main",     dis = "PD",
       file = "GSE72774/GSE72774_datBetaNormalized.csv.gz"),
  list(gse = "GSE111629", role = "PD_super",    dis = "PD",
       file = "GSE111629/GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz"),
  list(gse = "GSE112611", role = "IBD_peds",    dis = "IBD",
       file = "GSE112611/GSE112611_beta_values.txt.gz"),
  list(gse = "GSE152027", role = "neg_control",  dis = "SCZ",
       file = "GSE152027/GSE152027_IOP_processed_signals.csv.gz"),
  list(gse = "GSE145361", role = "PD_large",    dis = "PD",
       file = "GSE145361/GSE145361_Vallerga2020_NCOMMS_AvgBeta_Matrix-file.txt.gz")
)

# With limited memory, processing can be restricted here to some of the sets
ONLY <- NULL   # e.g. c("GSE87648", "GSE72774")

report <- character(0)
add <- function(...) report <<- c(report, paste0(...))

probe_lists <- list()

# Before starting, confirm that source files exist at all. If none is found, stop at once and print what
# the directory actually contains, rather than reporting "source file missing" six times and ending with an empty report -
# that kind of output hides the real cause, a wrong path.
.present <- vapply(COHORTS, function(co) file.exists(file.path(DATA_ROOT, co$file)),
                   logical(1))
.cached  <- vapply(COHORTS, function(co)
  file.exists(file.path(RDS_DIR, paste0(co$gse, "_beta.rds"))), logical(1))
if (!any(.present | .cached)) {
  stop("No source file found under ", DATA_ROOT, ".\n",
       "Actual contents of this directory: ", paste(list.files(DATA_ROOT), collapse = ", "), "\n",
       "Please check that DATA_ROOT at the top of this script matches the location written by the download script.")
}
message("Source files found: ", sum(.present), "; already cached: ", sum(.cached), "")

for (co in COHORTS) {
  if (!is.null(ONLY) && !(co$gse %in% ONLY)) next

  src <- file.path(DATA_ROOT, co$file)
  rds <- file.path(RDS_DIR, paste0(co$gse, "_beta.rds"))
  pid <- file.path(RDS_DIR, paste0(co$gse, "_probes.rds"))

  message("\n==== ", co$gse, "  (", co$role, ") ====")

  if (file.exists(rds) && file.exists(pid)) {
    message(".rds exists, reading skipped")
    probe_lists[[co$gse]] <- readRDS(pid)
    add("### ", co$gse, " (", co$role, ")  - cached, probes ",
        length(probe_lists[[co$gse]]))
    next
  }

  if (!file.exists(src)) {
    message("  Source file does not exist, skipped: ", src)
    add("### ", co$gse, " (", co$role, ")  - source file missing, not processed")
    next
  }

  message("  Reading ", basename(src), "  (",
          round(file.size(src) / 1024^3, 2), " GB compressed)")
  dt <- fread(src, showProgress = TRUE)

  # First column = probe ID, remaining columns = samples
  id_col <- names(dt)[1]
  probes <- as.character(dt[[1]])
  dt[, (id_col) := NULL]

  samples <- names(dt)
  message("  Dimensions: ", length(probes), " probes x ", length(samples), " samples")
  message("  Name of first column \"", id_col, "\", probe ID examples: ",
          paste(head(probes, 3), collapse = ", "))

  # Convert all columns to numeric
  non_num <- names(dt)[!vapply(dt, is.numeric, logical(1))]
  if (length(non_num)) {
    message("  Non-numeric columns: ", length(non_num), "; coerced")
    dt[, (non_num) := lapply(.SD, function(x)
      suppressWarnings(as.numeric(as.character(x)))), .SDcols = non_num]
  }

  mat <- as.matrix(dt)
  rm(dt); invisible(gc())
  rownames(mat) <- probes
  colnames(mat) <- samples

  # Value range check
  smp_idx <- sample(seq_len(nrow(mat)), min(20000L, nrow(mat)))
  v <- as.numeric(mat[smp_idx, ])
  v <- v[is.finite(v)]
  rng <- range(v)
  na_frac <- mean(!is.finite(as.numeric(mat[smp_idx, ])))
  looks_beta <- rng[1] >= -0.001 && rng[2] <= 1.001

  message(sprintf("  value range %.4g to %.4g, missing %.3f%%",
                  rng[1], rng[2], 100 * na_frac))
  if (looks_beta) {
    message("  Decision: beta values")
  } else {
    message("  Decision: not beta values (most likely signal intensities); this set needs separate handling")
  }

  # Distribution of probe ID types
  n_cg <- sum(grepl("^cg", probes)); n_rs <- sum(grepl("^rs", probes))
  n_ch <- sum(grepl("^ch\\.", probes))
  message("  Probe types: cg ", n_cg, ", rs ", n_rs, ", ch ", n_ch,
          ", other ", length(probes) - n_cg - n_rs - n_ch)

  add("### ", co$gse, " (", co$role, ", ", co$dis, ")")
  add("  File: ", basename(src))
  add("  Dimensions: ", length(probes), " probes x ", length(samples), " samples")
  add("  Name of first column: ", id_col)
  add("  Sample ID examples: ", paste(head(samples, 3), collapse = " | "))
  add("  Probe ID examples: ", paste(head(probes, 3), collapse = ", "))
  add(sprintf("  Value range: %.4g to %.4g   missing %.3f%%   %s",
              rng[1], rng[2], 100 * na_frac,
              if (looks_beta) "classified as beta values" else "not beta values, needs separate handling"))
  add("  Probe types: cg ", n_cg, ", rs ", n_rs, ", ch ", n_ch)
  add("")

  saveRDS(mat, rds, compress = FALSE)
  saveRDS(probes, pid)
  message("  Saved ", basename(rds), "  (",
          round(file.size(rds) / 1024^3, 2), " GB)")

  probe_lists[[co$gse]] <- probes
  rm(mat); invisible(gc())
}

# --------------------------------------------------------------
# Common probes
# --------------------------------------------------------------
if (length(probe_lists) >= 2L) {
  message("\n==== Common probes ====")
  # Only probes starting with cg enter the intersection: rs are genotype probes, ch are non-CpG sites
  cg_lists <- lapply(probe_lists, function(p) p[grepl("^cg", p)])
  common <- Reduce(intersect, cg_lists)
  add("### Common probes")
  for (g in names(cg_lists)) {
    add("  ", g, " cg probes: ", length(cg_lists[[g]]))
  }
  add("  Data sets in the intersection: ", length(cg_lists), "")
  add("  Common cg probes: ", length(common))
  message("  Sets included: ", length(cg_lists), "; common cg probes: ", length(common))

  # Accumulate the intersection set by set to see which set reduces it
  acc <- cg_lists[[1]]; nm <- names(cg_lists)[1]
  add("  Cumulative intersection, set by set: ")
  add("    ", nm, " -> ", length(acc))
  for (i in 2:length(cg_lists)) {
    acc <- intersect(acc, cg_lists[[i]])
    add("    + ", names(cg_lists)[i], " -> ", length(acc))
    message("    + ", names(cg_lists)[i], " -> ", length(acc))
  }

  saveRDS(common, file.path(RDS_DIR, "common_cg_probes.rds"))
  writeLines(common, file.path(OUT_DIR, "common_cg_probes.txt"), useBytes = TRUE)
} else {
  add("### Common probes: fewer than 2 usable data sets, not computed")
}

rep_path <- file.path(OUT_DIR, "04_load_report.txt")
writeLines(report, rep_path, useBytes = TRUE)
message("\nReport: ", rep_path)
message(".rds cache: ", RDS_DIR)
message("\nNext: review 04_load_report.txt.")
