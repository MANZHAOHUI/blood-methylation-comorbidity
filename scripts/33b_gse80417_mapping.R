# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 33b: sample correspondence of GSE80417
#
# Five of the six sets are built; only this one remains. Its column names are array barcodes
# (9704031032_R02C01) whereas the sample titles are study IDs (C0301), and the two do not match;
# the header has no sentrixids field as in GSE84727 either.
#
# The previous log printed only the six strategies with the highest success rates, and the supplementary-file fallback was last among the candidates,
# so it could not be seen whether it scored zero or was never evaluated - an omission in the output design;
# this step lists all of them.
#
# If no field of the header carries barcodes, the only remaining possibility is "the order of matrix columns equals the order of GEO
# samples". Identical order must not be assumed; it must be verified. The verification
# infers sex from the data themselves: the proportion of X-chromosome probes with intermediate methylation (0.3 to 0.7)
# is clearly higher in women than in men - women have one inactivated X whose CpG islands are intermediately
# methylated, whereas men have a single X and most probes are near 0 or 1. The inferred sex is compared sample by sample with the sex
# recorded in the header, and the order is accepted only at a concordance of 95% or more.
#
# This step changes no data; it only decides whether the correspondence is credible.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/33b_gse80417_mapping.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "33b_gse80417_mapping"

RAW_DIR   <- file.path(DATA_ROOT, "raw")
PROBE_DIR <- file.path(DATA_ROOT, "probe_lists")
OUT_DIR   <- file.path(PROJ_ROOT, "output", "33_build")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
suppressPackageStartupMessages(library(data.table))

GSE <- "GSE80417"
BARCODE <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
CONC_MIN <- 0.95        # Minimum required sex concordance

sl <- function(x) base::trimws(base::gsub('^"|"$', "",
  base::strsplit(x, "\t", fixed = TRUE)[[1]]))

hf <- file.path(RAW_DIR, GSE, base::paste0(GSE, "_series_matrix.txt.gz"))
con <- gzfile(hf, "rt"); hd <- readLines(con, n = 2000L, warn = FALSE); close(con)
gsm <- sl(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
ttl <- sl(base::grep("^!Sample_title", hd, value = TRUE)[1])[-1]
say("Samples: ", length(gsm))

# ==============================================================
say("")
say("=== 1. Sample-level fields in the header ===")
tags <- base::unique(base::sub("\t.*$", "", base::grep("^!Sample_", hd, value = TRUE)))
for (t in tags) {
  ln <- base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)
  for (one in ln) {
    v <- sl(one)[-1]
    nbc <- base::sum(base::regexpr(BARCODE, v) > 0)
    say(sprintf("  %-34s %d values; %d with barcodes; e.g. %s",
                t, length(v), nbc,
                base::paste(utils::head(v, 2), collapse = " | ")))
  }
}

# ==============================================================
say("")
say("=== 2. Checking the matrix column names strategy by strategy ===")
f <- list.files(file.path(RAW_DIR, GSE), pattern = "normalizedBetas", full.names = TRUE)
if (!length(f)) stop("beta matrix of GSE80417 not found")
dt <- data.table::fread(f[1], header = TRUE, check.names = FALSE,
                        nrows = 5L, showProgress = FALSE)
cn <- base::names(dt)[-1]
say("  Matrix columns: ", length(cn), "; e.g. ", base::paste(utils::head(cn, 3), collapse = " | "))
cands <- list(list(name = "GSM accession", key = gsm), list(name = "sample title", key = ttl))
for (t in tags) {
  for (one in base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)) {
    v <- sl(one)[-1]
    if (length(v) != length(gsm)) next
    cands[[length(cands) + 1]] <- list(name = base::paste0(t, " raw value"), key = v)
    k <- base::rep(NA_character_, length(v))
    hit <- base::regexpr(BARCODE, v) > 0
    if (base::any(hit)) {
      k[hit] <- base::regmatches(v, base::regexpr(BARCODE, v))
      cands[[length(cands) + 1]] <- list(name = base::paste0(t, " barcode within"), key = k)
    }
    vv <- base::trimws(base::sub("^[^:]*:", "", v))
    cands[[length(cands) + 1]] <- list(name = base::paste0(t, " after colon"), key = vv)
  }
}
best <- list(frac = -1)
for (one in cands) {
  fr <- base::mean(cn %in% one$key)
  say(sprintf("  %-44s %.3f", one$name, fr))
  if (fr > best$frac) best <- c(one, list(frac = fr))
}
say("  Highest success rate: ", sprintf("%.3f", best$frac), " (", best$name, ")")

found_field <- best$frac >= CONC_MIN
if (found_field) {
  say("")
  say("  A matching field exists in the header (", best$name, "); it can be used directly and order verification is not needed.")
  say("  The mapping strategy of 33 now scans all !Sample_ fields; simply re-run 33.")
}

# Order verification is used only if the header has no usable field
if (!found_field) {
# ==============================================================
say("")
say("=== 3. The assumption of identical order: verified by inferring sex from methylation ===")
say("  The header has no usable matching field; the only remaining possibility is 'order of matrix columns = order of GEO samples'.")
say("  It is tested below with the data themselves rather than assumed.")
sex_rec <- NULL
for (one in base::grep("^!Sample_characteristics_ch1", hd, value = TRUE)) {
  v <- sl(one)[-1]
  if (base::any(base::grepl("^Sex:", v, ignore.case = TRUE))) {
    sex_rec <- base::toupper(base::substr(base::trimws(base::sub("^[^:]*:", "", v)), 1, 1))
    break
  }
}
if (is.null(sex_rec)) {
  say("  [FAIL] no sex field in the header; verification not possible")
} else {
say("  Recorded sex: ", base::paste(sprintf("%s=%d", base::names(base::table(sex_rec)),
                                          base::as.integer(base::table(sex_rec))),
                                  collapse = "; "))

anno <- base::readRDS(file.path(PROBE_DIR, "anno450k.rds"))
pid <- base::as.character(base::rownames(anno))
if (!length(pid) || base::mean(base::grepl("^cg", pid)) < 0.5)
  pid <- base::as.character(anno[["Name"]])
chrX <- pid[base::as.character(anno[["chr"]]) == "chrX"]
say("  X-chromosome probes: ", length(chrX))

say("  Reading the full matrix ...")
dt <- data.table::fread(f[1], header = TRUE, check.names = FALSE, showProgress = FALSE)
probes <- base::as.character(dt[[1]])
ix <- base::which(probes %in% chrX)
say("  X probes in the matrix: ", length(ix))
MX <- base::as.matrix(dt[ix, -1, with = FALSE])
rm(dt); invisible(gc())

frac_mid <- base::colMeans(MX > 0.3 & MX < 0.7, na.rm = TRUE)
say("  Proportion of X probes with intermediate methylation: median ", sprintf("%.3f", stats::median(frac_mid)),
    ", range ", sprintf("%.3f to %.3f", base::min(frac_mid), base::max(frac_mid)))
cl <- stats::kmeans(frac_mid, centers = 2)
hi <- base::which.max(cl$centers)
pred <- base::ifelse(cl$cluster == hi, "F", "M")
say("  Centres of the two clusters: ", base::paste(sprintf("%.3f", base::sort(cl$centers)), collapse = " and "),
    "; split into ", base::sum(pred == "F"), " female / ", base::sum(pred == "M"), " male")

n <- base::min(length(pred), length(sex_rec))
conc <- base::mean(pred[1:n] == sex_rec[1:n])
say("  Compared under 'column order = sample order', sex concordance ", sprintf("%.3f", conc))
# Shuffling as a control: a similar concordance after shuffling would mean the test has no discriminating power
set.seed(1); conc_perm <- base::mean(base::sample(pred[1:n]) == sex_rec[1:n])
say("  Concordance after random shuffling: ", sprintf("%.3f", conc_perm), " (as a control)")

if (conc >= CONC_MIN) {
  say("  PASS: the order assumption holds; columns can be mapped to GEO samples by order")
  base::write.csv(data.frame(col_index = 1:n, col_name = cn[1:n], gsm = gsm[1:n],
                             sex_recorded = sex_rec[1:n], sex_predicted = pred[1:n],
                             stringsAsFactors = FALSE),
                  file.path(OUT_DIR, "GSE80417_order_map.csv"), row.names = FALSE)
  say("  Mapping table written, to be read directly by 33")
} else {
  say("  [FAIL] the order assumption does not hold. Samples of this data set cannot be matched reliably; consider using")
  say("  other schizophrenia cohorts instead (GSE152027 is already downloaded, GSE147221 can be obtained separately)")
}
}   # end of the sex field branch
}   # end of the branch "no usable field in the header"
log_stop()
