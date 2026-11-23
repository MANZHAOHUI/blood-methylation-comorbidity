# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 06: mapping of sample columns + two probe sets
#
# Part A: map the column names of each matrix to the GEO sample accessions (GSM).
#   The five matrices use five naming schemes: array barcodes, Sample1..N, hospital IDs,
#   whereas case-control status, age, sex and smoking are all attached to the GSM.
#   Without a mapping no model can be fitted; a wrong mapping raises no error and only assigns beta values to the wrong person.
#   This step therefore reports the match rate per set and lists unmatched columns and absent GSMs.
#
# Part B: determine the column layout of each matrix and flag detection P value columns and empty columns.
# Part C: build the two probe sets, primary analysis set and replication set, according to decision D11.
#
# This script does not read beta matrices; it reads only the first two lines of each file and the cached probe lists.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/06_map_samples_and_probes.R
# ==============================================================

DATA_ROOT <- "D:/blood_methylation_data"
PROJ_ROOT <- "C:/blood_methylation_project"

RDS_DIR <- file.path(DATA_ROOT, "rds")
HDR_DIR <- file.path(PROJ_ROOT, "output", "01_metadata")
OUT_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
for (d in c(RDS_DIR, HDR_DIR)) if (!dir.exists(d)) stop("Directory does not exist: ", d)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

FILES <- list(
  GSE87648  = list(path = "GSE87648/GSE87648_processedMethWb.txt.gz",  sep = "\t"),
  GSE72774  = list(path = "GSE72774/GSE72774_datBetaNormalized.csv.gz", sep = ","),
  GSE111629 = list(path = "GSE111629/GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz", sep = "\t"),
  GSE112611 = list(path = "GSE112611/GSE112611_beta_values.txt.gz",     sep = "\t"),
  GSE152027 = list(path = "GSE152027/GSE152027_IOP_processed_signals.csv.gz", sep = ","),
  GSE145361 = list(path = "GSE145361/GSE145361_Vallerga2020_NCOMMS_AvgBeta_Matrix-file.txt.gz", sep = "\t")
)

PLATFORM_450K <- c("GSE87648", "GSE72774", "GSE111629", "GSE152027", "GSE145361")
PLATFORM_EPIC <- "GSE112611"

BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"

report <- character(0)
add <- function(...) report <<- c(report, paste0(...))

# --------------------------------------------------------------
# Read the first two lines: header + first data row
# --------------------------------------------------------------
first_two <- function(gse) {
  p <- file.path(DATA_ROOT, FILES[[gse]]$path)
  if (!file.exists(p)) return(NULL)
  con <- gzfile(p, "rt"); on.exit(close(con), add = TRUE)
  l <- readLines(con, n = 2L, warn = FALSE)
  if (length(l) < 2L) return(NULL)
  sp <- function(x) gsub('^"|"$', "", strsplit(x, FILES[[gse]]$sep, fixed = TRUE)[[1]])
  list(header = sp(l[1]), row1 = sp(l[2]))
}

is_detection <- function(nm) grepl("detection|pval|p_?value", nm, ignore.case = TRUE)
ID_LABEL <- "^(id|id_ref|sampleid|targetid|probe.*|cpg.*|v1|x|)$"

# --------------------------------------------------------------
# Decide the column layout: is the first header field the name of the probe column or the first sample?
# --------------------------------------------------------------
layout_of <- function(h, r1) {
  n_h <- length(h); n_d <- length(r1)
  if (grepl(ID_LABEL, tolower(trimws(h[1])))) {
    # The first header field names the probe column; the others correspond to data columns 2..n
    data.frame(name = h[-1], pos = seq(2L, n_h), stringsAsFactors = FALSE)
  } else if (n_d == n_h + 1L) {
    # The header gives no name for the probe column; all header fields are samples, shifted right by one
    data.frame(name = h, pos = seq(2L, n_h + 1L), stringsAsFactors = FALSE)
  } else {
    data.frame(name = h[-1], pos = seq(2L, n_h), stringsAsFactors = FALSE)
  }
}

# --------------------------------------------------------------
# Metadata: title, description, supplementary file name
# --------------------------------------------------------------
geo_meta <- function(gse) {
  f <- file.path(HDR_DIR, paste0(gse, "_header.txt"))
  if (!file.exists(f)) return(NULL)
  l <- readLines(f, warn = FALSE, encoding = "UTF-8")
  fld <- function(tag, exact = TRUE) {
    pat <- if (exact) paste0("^", tag, "\t") else paste0("^", tag)
    x <- grep(pat, l, value = TRUE)
    if (!length(x)) return(NULL)
    lapply(x, function(one)
      trimws(gsub('^"|"$', "", strsplit(one, "\t", fixed = TRUE)[[1]])[-1]))
  }
  acc <- fld("!Sample_geo_accession")
  if (is.null(acc)) return(NULL)
  acc <- acc[[1]]
  n <- length(acc)
  grab <- function(tag) {
    v <- fld(tag, exact = FALSE)
    if (is.null(v)) return(rep("", n))
    v <- Filter(function(z) length(z) == n, v)
    if (!length(v)) return(rep("", n))
    do.call(paste, c(v, sep = " ; "))
  }
  data.frame(gsm = acc,
             title = grab("!Sample_title"),
             description = grab("!Sample_description"),
             suppl = grab("!Sample_supplementary_file"),
             stringsAsFactors = FALSE)
}

# --------------------------------------------------------------
# Several possible spellings of one string
# --------------------------------------------------------------
norm <- function(x) toupper(trimws(x))

variants_meta <- function(title, description, suppl) {
  s <- c(title, description, suppl)
  s <- s[nzchar(s)]
  out <- character(0)
  for (x in s) {
    out <- c(out,
             x,
             sub("\\s*\\[[^]]*\\]\\s*$", "", x),          # Remove the trailing [Whole blood]
             sub("^genomic DNA from\\s+", "", x),
             sub("^.*subject\\s+", "", x),
             sub("^([^ ;]+).*$", "\\1", x),               # First word
             sub("^.*[ ;]([^ ;]+)$", "\\1", x))           # Last word
    m <- regmatches(x, gregexpr(BARCODE, x))[[1]]         # Array barcode at any position
    if (length(m)) out <- c(out, m)
  }
  unique(norm(out[nzchar(out)]))
}

variants_col <- function(x) {
  out <- c(x,
           sub("^X", "", x),                              # Prefix added by R to make names syntactic
           sub("\\s*\\[[^]]*\\]\\s*$", "", x),
           sub("^([^ ]+).*$", "\\1", x))
  m <- regmatches(x, gregexpr(BARCODE, x))[[1]]
  if (length(m)) out <- c(out, m)
  unique(norm(out[nzchar(out)]))
}

# --------------------------------------------------------------
# Main loop
# --------------------------------------------------------------
for (gse in names(FILES)) {
  message("\n==== ", gse, " ====")
  ft <- first_two(gse)
  if (is.null(ft)) {
    message("  Source file does not exist or has fewer than two lines, skipped")
    add("### ", gse, "  Source file not usable"); add(""); next
  }
  meta <- geo_meta(gse)
  if (is.null(meta)) {
    message("  Cached header from 01 missing, skipped")
    add("### ", gse, "  Cached header from 01 missing"); add(""); next
  }

  cols <- layout_of(ft$header, ft$row1)
  n_det <- sum(is_detection(cols$name))
  n_blank <- sum(!nzchar(trimws(cols$name)))
  cols <- cols[nzchar(trimws(cols$name)) & !is_detection(cols$name), , drop = FALSE]

  message("  Header fields ", length(ft$header), ", data columns ", length(ft$row1),
          "; detection P value columns ", n_det, ", empty columns ", n_blank,
          "; sample columns ", nrow(cols))
  message("  GEO metadata: ", nrow(meta), " samples")

  # Build a "spelling -> GSM" lookup table; a spelling that points to several GSMs is ambiguous and removed
  key_list <- lapply(seq_len(nrow(meta)), function(i)
    variants_meta(meta$title[i], meta$description[i], meta$suppl[i]))
  lut <- data.frame(key = unlist(key_list, use.names = FALSE),
                    gsm = rep(meta$gsm, lengths(key_list)),
                    stringsAsFactors = FALSE)
  lut <- unique(lut)
  ambig <- names(which(table(lut$key) > 1L))
  if (length(ambig)) lut <- lut[!(lut$key %in% ambig), , drop = FALSE]
  lookup <- setNames(lut$gsm, lut$key)

  matched <- rep(NA_character_, nrow(cols))
  for (i in seq_len(nrow(cols))) {
    v <- variants_col(cols$name[i])
    hit <- unique(na.omit(lookup[v]))
    if (length(hit) == 1L) matched[i] <- hit
  }

  n_ok <- sum(!is.na(matched))
  n_dup <- sum(duplicated(matched[!is.na(matched)]))
  message("  Matched ", n_ok, " / ", nrow(cols), ", duplicate matches ", n_dup)

  add("### ", gse)
  add("  Header fields ", length(ft$header), "; data columns ", length(ft$row1),
      "; detection P value columns ", n_det, "; empty columns ", n_blank, "; sample columns ", nrow(cols))
  add("  GEO metadata: ", nrow(meta), " samples; ambiguous spellings ", length(ambig), "")
  add("  Matched ", n_ok, " / ", nrow(cols), "; duplicate matches ", n_dup)
  if (n_ok < nrow(cols)) {
    bad <- cols$name[is.na(matched)]
    add("  Unmatched matrix columns: ", length(bad), " (at most 8 shown): ",
        paste(head(bad, 8), collapse = " | "))
    message("  Examples of unmatched columns: ", paste(head(bad, 5), collapse = " | "))
  }
  miss <- setdiff(meta$gsm, matched)
  if (length(miss)) {
    add("  GSMs absent from the matrix: ", length(miss), " (at most 8 shown): ",
        paste(head(miss, 8), collapse = ", "))
    add("    their titles: ", paste(head(meta$title[meta$gsm %in% miss], 3), collapse = " | "))
    add("    their descriptions: ", paste(head(meta$description[meta$gsm %in% miss], 2), collapse = " | "))
  }
  add("")

  write.csv(data.frame(col_pos_in_data = cols$pos,
                       col_name = cols$name,
                       gsm = matched,
                       stringsAsFactors = FALSE),
            file.path(OUT_DIR, paste0(gse, "_sample_map.csv")), row.names = FALSE)
}

# --------------------------------------------------------------
# Part C: two probe sets
# --------------------------------------------------------------
message("\n==== Probe sets ====")
get_cg <- function(gse) {
  f <- file.path(RDS_DIR, paste0(gse, "_probes.rds"))
  if (!file.exists(f)) return(NULL)
  p <- readRDS(f); p[grepl("^cg", p)]
}
main_sets <- lapply(PLATFORM_450K, get_cg); names(main_sets) <- PLATFORM_450K
main_sets <- main_sets[!vapply(main_sets, is.null, logical(1))]
add("  Data sets in the primary analysis set: ", paste(names(main_sets), collapse = ", "))

add("### Probe sets")
if (length(main_sets) >= 2L) {
  main_common <- Reduce(intersect, main_sets)
  for (g in names(main_sets)) add("  ", g, " cg probes: ", length(main_sets[[g]]))
  add("  Primary analysis set (intersection of", length(main_sets), " 450K sets): ", length(main_common))
  message("  Primary analysis set: ", length(main_common))
  saveRDS(main_common, file.path(RDS_DIR, "probes_main_450k.rds"))
  epic <- get_cg(PLATFORM_EPIC)
  if (!is.null(epic)) {
    repl <- intersect(main_common, epic)
    add("  ", PLATFORM_EPIC, " cg probes: ", length(epic))
    add("  Replication set (after adding EPIC): ", length(repl),
        "; loss caused by the addition ", length(main_common) - length(repl))
    message("  Replication set: ", length(repl))
    saveRDS(repl, file.path(RDS_DIR, "probes_replication_epic.rds"))
  }
} else {
  add("  Not enough usable data sets, not computed")
}

rep_path <- file.path(OUT_DIR, "06_mapping_report.txt")
writeLines(report, rep_path, useBytes = TRUE)
message("\nReport: ", rep_path)
message("\nNext: review 06_mapping_report.txt.")
