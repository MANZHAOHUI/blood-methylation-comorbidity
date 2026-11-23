# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 07: build the data objects used directly for modelling
#
# Four things are done for each data set:
#   (1) select the sample columns according to the mapping from 06, drop detection P value columns, rename columns to GSM
#   (2) restrict to the appropriate probe set according to decision D11 (primary analysis set for the four 450K sets,
#       replication set for the paediatric EPIC cohort)
#   (3) attach phenotypes and define a uniform group (case / control);
#       samples that cannot be assigned are removed and reported
#   (4) check missingness at sample and probe level
#
# Only one set in memory at a time. Output is one .rds per set, containing the beta matrix and the phenotype table.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/07_build_analysis_sets.R
# ==============================================================

DATA_ROOT <- "D:/blood_methylation_data"
PROJ_ROOT <- "C:/blood_methylation_project"

RDS_DIR  <- file.path(DATA_ROOT, "rds")
MAP_DIR  <- file.path(PROJ_ROOT, "output", "06_mapping")
META_DIR <- file.path(PROJ_ROOT, "output", "02_overlap")
OUT_DIR  <- file.path(PROJ_ROOT, "output", "07_analysis_sets")
SET_DIR  <- file.path(DATA_ROOT, "analysis_sets")
for (d in c(RDS_DIR, MAP_DIR, META_DIR)) if (!dir.exists(d)) stop("Directory does not exist: ", d)
for (d in c(OUT_DIR, SET_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
}

ANCESTRY_F <- file.path(DATA_ROOT, "VALLERGA", "Vallerga2020_NCOMMS_ancestry.txt")
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

probes_main <- readRDS(file.path(RDS_DIR, "probes_main_450k.rds"))
probes_repl <- readRDS(file.path(RDS_DIR, "probes_replication_epic.rds"))
message("Primary analysis set: ", length(probes_main), " probes; replication set: ", length(probes_repl), " probes")

# --------------------------------------------------------------
# Specification of each data set: group field, values for cases and controls, covariates kept
# --------------------------------------------------------------
SPEC <- list(
  GSE87648 = list(
    probes = "main", disease = "IBD",
    group_col = "simplified_diagnosis",
    case = c("CD", "UC"), control = c("HL", "HS"),
    keep = c("sex", "age", "smoking.status", "scan_date",
             "simplified_diagnosis", "patient_number")),
  GSE72774 = list(
    probes = "main", disease = "PD",
    group_col = "diseasestatus",
    case = "PD", control = "control",
    keep = c("sex", "age", "ethnicity", "dnamage",
             "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
             "smokingstatus", "packyearstotal",
             "levodopastatus", "levodopaperdayincludingagonists")),
  GSE111629 = list(
    probes = "main", disease = "PD",
    group_col = "disease.state",
    case = "Parkinson's disease (PD)", control = "PD-free control",
    keep = c("gender", "age", "ethnicity")),
  GSE152027 = list(
    probes = "main", disease = "SCZ",
    group_col = "status",
    case = "SCZ", control = "CON",
    keep = c("gender", "ageatbloodcollection")),
  GSE145361 = list(
    probes = "main", disease = "PD",
    group_col = "disease.state",
    case = "Parkinson's disease", control = "Control",
    keep = c("gender", "PC1", "PC2"),
    use_ancestry = TRUE),
  GSE112611 = list(
    probes = "repl", disease = "IBD",
    group_col = "diagnosis",
    case = "Crohn's disease", control = "non-IBD control",
    keep = c("gender", "age", "baseline.vs.follow.up"),
    baseline_only = TRUE)
)

report <- character(0)
add <- function(...) report <<- c(report, paste0(...))

for (gse in names(SPEC)) {
  sp <- SPEC[[gse]]
  message("\n==== ", gse, " (", sp$disease, ") ====")

  beta_f <- file.path(RDS_DIR, paste0(gse, "_beta.rds"))
  map_f  <- file.path(MAP_DIR, paste0(gse, "_sample_map.csv"))
  met_f  <- file.path(META_DIR, paste0(gse, "_metadata_bytag.csv"))
  miss_f <- c(beta_f, map_f, met_f)[!file.exists(c(beta_f, map_f, met_f))]
  if (length(miss_f)) {
    message("  File missing, set skipped: ", paste(basename(miss_f), collapse = ", "))
    add("### ", gse, "  Not processed, missing: ", paste(basename(miss_f), collapse = ", ")); add("")
    next
  }

  map <- read.csv(map_f, stringsAsFactors = FALSE, colClasses = "character")
  met <- read.csv(met_f, stringsAsFactors = FALSE, colClasses = "character")
  if (anyNA(map$gsm) || any(!nzchar(map$gsm))) {
    stop(gse, ": the mapping table still has unmatched columns; please resolve step 06 first")
  }

  beta <- readRDS(beta_f)
  # The matrix saved by 04 no longer has the probe ID column, so data column k is column k-1 here
  idx <- as.integer(map$col_pos_in_data) - 1L
  if (max(idx) > ncol(beta) || min(idx) < 1L) {
    stop(gse, ": column position out of range: matrix ", ncol(beta), " columns, largest position ", max(idx))
  }
  beta <- beta[, idx, drop = FALSE]
  colnames(beta) <- map$gsm

  probes <- if (sp$probes == "main") probes_main else probes_repl
  hit <- intersect(rownames(beta), probes)
  if (length(hit) != length(probes)) {
    stop(gse, " lacks probes of the probe set: ", length(probes) - length(hit), "")
  }
  beta <- beta[probes, , drop = FALSE]

  # ---- phenotypes ----
  ph <- met[match(colnames(beta), met$geo_accession), , drop = FALSE]
  if (anyNA(ph$geo_accession)) stop(gse, ": some columns have no matching GSM in the metadata")
  gcol <- sp$group_col
  if (!gcol %in% names(ph)) stop(gse, ": the metadata has no group field ", gcol)

  grp <- rep(NA_character_, nrow(ph))
  grp[ph[[gcol]] %in% sp$case] <- "case"
  grp[ph[[gcol]] %in% sp$control] <- "control"
  ph$group <- factor(grp, levels = c("control", "case"))

  drop_reason <- character(0)
  keep_i <- !is.na(ph$group)
  if (any(!keep_i)) {
    tb <- table(ph[[gcol]][!keep_i])
    drop_reason <- c(drop_reason,
      paste0("not assigned to a group: ", sum(!keep_i), " (",
             paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = "; "), ")"))
  }
  if (isTRUE(sp$use_ancestry)) {
    if (!file.exists(ANCESTRY_F)) stop("The ancestry file of Vallerga et al. is missing: ", ANCESTRY_F)
    anc <- read.delim(ANCESTRY_F, stringsAsFactors = FALSE)
    anc$key <- toupper(trimws(anc$Sample_ID))
    # The array barcode may be in the matrix column names or in the GEO sample titles; both are tried
    src <- ifelse(nzchar(map$col_name), map$col_name, "")
    bc <- toupper(regmatches(src, regexpr(BARCODE, src)))
    got <- regexpr(BARCODE, src) > 0
    if (sum(got) < length(src) && "sample_title" %in% names(ph)) {
      alt <- ph$sample_title
      hit2 <- regexpr(BARCODE, alt) > 0
      bc2 <- rep(NA_character_, length(alt))
      bc2[hit2] <- toupper(regmatches(alt, regexpr(BARCODE, alt)))
      bc_full <- rep(NA_character_, length(src))
      bc_full[got] <- bc
      bc_full[!got] <- bc2[!got]
    } else {
      bc_full <- rep(NA_character_, length(src)); bc_full[got] <- bc
    }
    n_bc <- sum(!is.na(bc_full))
    final_ok <- anc$key[anc$Final_dataset == "Y"]
    in_final <- !is.na(bc_full) & bc_full %in% final_ok
    ph$PC1 <- suppressWarnings(as.numeric(anc$PC1[match(bc_full, anc$key)]))
    ph$PC2 <- suppressWarnings(as.numeric(anc$PC2[match(bc_full, anc$key)]))
    drop_reason <- c(drop_reason,
      paste0("barcodes parsed: ", n_bc, " / ", length(bc_full),
             "; removed because not in the final unrelated set of Vallerga et al.: ", sum(keep_i & !in_final), " samples"))
    keep_i <- keep_i & in_final
  }

  if (isTRUE(sp$baseline_only)) {
    bl <- ph$baseline.vs.follow.up == "BL"
    drop_reason <- c(drop_reason, paste0("not baseline: ", sum(keep_i & !bl), " samples"))
    keep_i <- keep_i & bl
  }

  beta <- beta[, keep_i, drop = FALSE]
  ph <- ph[keep_i, , drop = FALSE]

  keep_cols <- intersect(c("geo_accession", sp$keep), names(ph))
  missing_cols <- setdiff(sp$keep, names(ph))
  ph <- ph[, c(keep_cols, "group"), drop = FALSE]

  # ---- plausibility of age ----
  # The age column of GSE87648 contains 0; an adult IBD cohort has no 0-year-olds, so this is a missing-value code.
  # Unless set to NA, these points would distort the linear model without any error message.
  for (cl in grep("age", names(ph), ignore.case = TRUE, value = TRUE)) {
    v <- suppressWarnings(as.numeric(ph[[cl]]))
    bad <- !is.na(v) & (v <= 0 | v > 110)
    if (any(bad)) {
      add("  [", cl, "] implausible values: ", sum(bad), " set to missing (values: ",
          paste(head(sort(unique(v[bad])), 5), collapse = ", "), ")")
      message("  ", cl, ": ", sum(bad), " implausible values, set to missing")
      v[bad] <- NA_real_
    }
    ph[[cl]] <- v
  }

  # ---- missingness check ----
  na_by_sample <- colMeans(!is.finite(beta))
  na_by_probe  <- rowMeans(!is.finite(beta))
  bad_samp <- sum(na_by_sample > 0.05)
  bad_prob <- sum(na_by_probe > 0.05)

  tb <- table(ph$group)
  message("  controls ", tb[["control"]], ", cases ", tb[["case"]],
          "; probes ", nrow(beta))
  message("  samples with missingness >5%: ", bad_samp, ", probes with missingness >5%: ", bad_prob)

  add("### ", gse, "  (", sp$disease, ")")
  add("  Probe set: ", if (sp$probes == "main") "primary analysis set" else "replication set",
      "  ", nrow(beta), "")
  add("  controls ", tb[["control"]], "; cases ", tb[["case"]],
      "; total ", ncol(beta))
  if (length(drop_reason)) add("  Removed: ", paste(drop_reason, collapse = "; "))
  if (length(missing_cols)) add("  Covariates absent from the metadata: ",
                               paste(missing_cols, collapse = ", "))
  add("  Sample missing rate: median ", sprintf("%.4f", median(na_by_sample)),
      ", maximum ", sprintf("%.4f", max(na_by_sample)),
      "; samples with missingness >5%: ", bad_samp)
  add("  Probe missing rate: maximum ", sprintf("%.4f", max(na_by_probe)),
      "; probes with missingness >5%: ", bad_prob)
  for (cl in setdiff(names(ph), c("geo_accession", "group"))) {
    v <- ph[[cl]]
    num <- suppressWarnings(as.numeric(v))
    if (sum(!is.na(num)) > 0.8 * length(v)) {
      add("  [", cl, "] numeric: median ", sprintf("%.3g", median(num, na.rm = TRUE)),
          ", missing ", sum(is.na(num)))
    } else {
      t2 <- sort(table(v, useNA = "ifany"), decreasing = TRUE)
      t2 <- t2[seq_len(min(8L, length(t2)))]
      add("  [", cl, "] ", paste(sprintf("%s=%d", names(t2), as.integer(t2)),
                                 collapse = "; "))
    }
  }
  add("")

  saveRDS(list(beta = beta, pheno = ph, gse = gse, disease = sp$disease,
               probe_set = sp$probes),
          file.path(SET_DIR, paste0(gse, "_analysis.rds")), compress = FALSE)
  message("  Saved ", gse, "_analysis.rds")
  rm(beta, ph); invisible(gc())
}

rep_path <- file.path(OUT_DIR, "07_analysis_sets_report.txt")
writeLines(report, rep_path, useBytes = TRUE)
message("\nReport: ", rep_path)
message("Objects for modelling: ", SET_DIR)
message("\nNext: review 07_analysis_sets_report.txt.")
