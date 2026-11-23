# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 33: build the analysis sets of the six additional cohorts (written to the actual file layouts)
#
# 33a established that the six data sets come in four different file layouts; this step handles them as they are, without assumptions:
#
#   GSE84727 / GSE80417   csv, first header field empty (the probe column is unnamed);
#                         column names are array barcodes whereas sample titles are study IDs,
#                         so the barcodes must be found in the characteristic fields
#   GSE198904             header with 864 fields, data rows with 865 - the header lacks the first field,
#                         and fread with header=TRUE would take the first data row as the header
#   GSE273970             column names of the form sample_N.Noob_Beta and
#                         sample_N.Detection_Pval, arranged in blocks, with gaps in the numbering
#   GSE338167             column names are barcodes, hidden at the end of the sample titles
#   GSE144858             a series matrix; skip to the line after the data start marker
#
# The sample mapping is not hard-coded: five strategies are tried in turn and the successful one is reported:
#   GSM directly, title directly, barcode from the title by regular expression, scan of all characteristic fields,
#   sample_N numbering. A success rate below 95% is a FAIL.
#
# An error of 33 is also fixed: the probe universe came out as 0 because grepl was called unqualified while
# BiocGenerics masks grep and grepl - the masking notice of 00_logging.R lists
# both functions. This step uses the base:: prefix throughout.
#
# Sets are processed and saved one by one; a failure midway does not affect cohorts already finished, and a re-run skips them.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/33_build_new_analysis_sets.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "33_build_new_analysis_sets"

RAW_DIR   <- file.path(DATA_ROOT, "raw")
PROBE_DIR <- file.path(DATA_ROOT, "probe_lists")
SET_DIR   <- file.path(DATA_ROOT, "analysis_sets")
PHENO_DIR <- file.path(PROJ_ROOT, "output", "32_pheno")
OUT_DIR   <- file.path(PROJ_ROOT, "output", "33_build")
for (d in c(SET_DIR, OUT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
suppressPackageStartupMessages(library(data.table))

FORCE     <- FALSE
MAF_CUT   <- 0.01
MISS_CUT  <- 0.05
MAP_MIN   <- 0.95
BARCODE   <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
BETA_SUF  <- "\\.(Noob_)?Beta$|\\.AVG_Beta$|_Beta$"
PVAL_SUF  <- "\\.Detection_Pval$|\\.Detection\\.Pval$|\\.pval$|_Pval$"

COHORTS <- list(
  GSE84727  = list(disease = "schizophrenia", platform = "450K"),
  GSE80417  = list(disease = "schizophrenia", platform = "450K"),
  GSE198904 = list(disease = "depression",    platform = "EPIC"),
  GSE273970 = list(disease = "depression",    platform = "EPIC"),
  GSE338167 = list(disease = "alzheimer",     platform = "EPICv2"),
  GSE144858 = list(disease = "alzheimer",     platform = "450K")
)

# ==============================================================
# 1. Unified probe universe
# ==============================================================
say("=== 1. Unified probe universe filtered according to D2 (450K annotation as the full set) ===")
anno <- readRDS(file.path(PROBE_DIR, "anno450k.rds"))
say("  Annotation object: ", base::paste(class(anno), collapse = "/"), ", ",
    base::paste(dim(anno), collapse = " x "))
u <- base::as.character(base::rownames(anno))
if (!length(u) || base::mean(base::grepl("^cg", u)) < 0.5) {
  u <- base::as.character(anno[["Name"]])
  say("  Row names not usable; probe IDs taken from the Name column")
}
u <- u[base::grepl("^cg", u)]
say("  cg probes: ", length(u))
if (!length(u)) stop("Probe IDs cannot be obtained; please inspect the annotation object manually")

idx <- base::match(u, base::as.character(base::rownames(anno)))
if (base::anyNA(idx)) idx <- base::match(u, base::as.character(anno[["Name"]]))
chr <- base::as.character(anno[["chr"]])[idx]
drop_sex <- u[chr %in% c("chrX", "chrY")]
say("  Sex-chromosome probes: ", length(drop_sex))

snp_hit <- base::rep(FALSE, length(u))
for (cc in c("CpG_maf", "SBE_maf")) {
  if (!(cc %in% base::colnames(anno))) next
  v <- base::suppressWarnings(base::as.numeric(anno[[cc]])[idx])
  snp_hit <- snp_hit | (!base::is.na(v) & v > MAF_CUT)
}
drop_snp <- u[snp_hit]
say("  MAF at CpG or SBE > ", MAF_CUT, ": ", length(drop_snp))

xr <- character(0)
for (fn in list.files(PROBE_DIR, pattern = "(non-specific|multimap)",
                      full.names = TRUE, ignore.case = TRUE)) {
  x <- tryCatch(readLines(fn, warn = FALSE), error = function(e) character(0))
  ids <- base::unique(base::unlist(base::regmatches(x, base::gregexpr("cg[0-9]{8}", x))))
  say("  ", base::basename(fn), ": ", length(ids), "")
  xr <- base::union(xr, ids)
}
drop_xr <- base::intersect(u, xr)
say("  Cross-reactive and multi-mapping in total: ", length(xr), ", of which on 450K: ", length(drop_xr))

UNIVERSE <- base::setdiff(u, base::union(base::union(drop_sex, drop_snp), drop_xr))
say("  Probe universe after filtering: ", length(UNIVERSE), " (removed ", length(u) - length(UNIVERSE),
    ", ", sprintf("%.1f%%", 100 * (1 - length(UNIVERSE) / length(u))), ")")
if (length(UNIVERSE) < 3e5) stop("The probe universe is abnormally small; please inspect the filtering steps manually")
base::saveRDS(UNIVERSE, file.path(PROBE_DIR, "universe_450k_filtered.rds"))
say("")

# ==============================================================
# Reading utilities
# ==============================================================
sl <- function(x, sep = "\t")
  base::trimws(base::gsub('^"|"$', "", base::strsplit(x, sep, fixed = TRUE)[[1]]))

read_hdr <- function(g) {
  f <- file.path(RAW_DIR, g, paste0(g, "_series_matrix.txt.gz"))
  con <- gzfile(f, "rt"); on.exit(close(con), add = TRUE)
  readLines(con, n = 2000L, warn = FALSE)
}
first_lines <- function(f, n = 2L) {
  con <- gzfile(f, "rt"); on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

# General reader: look at the first two lines to decide the delimiter and whether the header lacks its first field
read_matrix <- function(f) {
  is_series <- base::grepl("series_matrix", base::basename(f), fixed = TRUE)
  n_skip <- 0L
  if (is_series) {
    con <- gzfile(f, "rt"); on.exit(close(con), add = TRUE)
    k <- 0L
    repeat {
      ln <- readLines(con, n = 1L, warn = FALSE)
      if (!length(ln)) stop("Data start marker not found: ", base::basename(f))
      k <- k + 1L
      if (base::grepl("!series_matrix_table_begin", ln, fixed = TRUE)) break
    }
    close(con); on.exit()
    n_skip <- k
  }
  two <- if (n_skip == 0L) first_lines(f, 2L) else {
    con <- gzfile(f, "rt"); on.exit(close(con), add = TRUE)
    x <- readLines(con, n = n_skip + 2L, warn = FALSE); close(con); on.exit()
    x[(n_skip + 1L):(n_skip + 2L)]
  }
  sep <- if (base::grepl("\t", two[1], fixed = TRUE)) "\t" else ","
  h <- sl(two[1], sep); d <- sl(two[2], sep)
  say("    delimiter ", if (sep == "\t") "tab" else "comma",
      "; header ", length(h), " fields, data rows ", length(d), " fields")
  shifted <- length(h) == length(d) - 1L
  if (shifted) say("    header lacks the first field, treated as shifted")
  dt <- data.table::fread(f, skip = n_skip + 1L, header = FALSE,
                          sep = sep, check.names = FALSE, showProgress = FALSE,
                          fill = TRUE)
  # The last line of a series matrix is the end marker
  if (is_series) {
    last <- base::as.character(dt[[1]][base::nrow(dt)])
    if (!base::is.na(last) && base::startsWith(last, "!")) dt <- dt[-base::nrow(dt), ]
  }
  nm <- if (shifted) c("probe_id", h) else h
  if (length(nm) != base::ncol(dt)) stop("Number of column names ", length(nm), " does not match number of columns ", base::ncol(dt), "")
  data.table::setnames(dt, nm)
  dt
}

# Sample mapping: five strategies tried in turn
map_samples <- function(cols, gsm, ttl, hd) {
  cand <- list(
    list(name = "GSM directly", key = gsm),
    list(name = "sample title directly", key = ttl),
    list(name = "array barcode in title",
         key = { m <- base::regmatches(ttl, base::regexpr(BARCODE, ttl))
                 k <- base::rep(NA_character_, length(ttl))
                 k[base::regexpr(BARCODE, ttl) > 0] <- m; k }),
    list(name = "sample_N number in title",
         key = base::sub("^\\s*sample[ _]?([0-9]+).*$", "sample_\\1", ttl,
                         ignore.case = TRUE))
  )
  # Scan all !Sample_ fields: the raw value, the part after the colon, and the array barcode within it.
  # No field name is assumed - the barcodes of GSE84727 are in sentrixids, those of GSE80417 in
  # !Sample_description, which cannot be guessed beforehand; a full scan is the safest.
  for (t in base::unique(base::sub("\t.*$", "", base::grep("^!Sample_", hd, value = TRUE)))) {
    for (one in base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)) {
      v <- sl(one)[-1]
      if (length(v) != length(gsm)) next
      nm <- base::sub("^!Sample_", "", t)
      cand[[length(cand) + 1]] <- list(name = paste0(nm, " raw value"), key = v)
      vv <- base::trimws(base::sub("^[^:]*:", "", v))
      if (!identical(vv, v))
        cand[[length(cand) + 1]] <- list(name = paste0(nm, " after colon"), key = vv)
      hit <- base::regexpr(BARCODE, v) > 0
      if (base::any(hit)) {
        k <- base::rep(NA_character_, length(v))
        k[hit] <- base::regmatches(v, base::regexpr(BARCODE, v))
        cand[[length(cand) + 1]] <- list(name = paste0(nm, " barcode within"), key = k)
      }
    }
  }
  best <- list(frac = -1); score <- list()
  for (one in cand) {
    k <- one$key
    if (length(k) != length(gsm)) next
    f <- base::mean(cols %in% k)
    score[[length(score) + 1]] <- data.frame(name = one$name, frac = f,
                                             stringsAsFactors = FALSE)
    if (f > best$frac) best <- list(frac = f, name = one$name, key = k)
  }
  best$score <- if (length(score)) base::do.call(rbind, score) else NULL
  best
}

# ==============================================================
# 2. Build the sets one by one
# ==============================================================
summ <- list()
for (g in names(COHORTS)) {
  C <- COHORTS[[g]]
  out_rds <- file.path(SET_DIR, paste0(g, "_analysis_filtered.rds"))
  say("==== ", g, " (", C$disease, ", ", C$platform, ") ====")
  if (file.exists(out_rds) && !FORCE) { say("  Exists, skipped (set FORCE to rebuild)"); next }

  ph <- read.csv(file.path(PHENO_DIR, paste0(g, "_pheno_clean.csv")), stringsAsFactors = FALSE)
  ph <- ph[ph$complete %in% c(TRUE, "TRUE"), , drop = FALSE]
  say("  Complete cases ", base::nrow(ph), ": controls ", base::sum(ph$group == "control"),
      ", cases ", base::sum(ph$group == "case"))

  hd <- read_hdr(g)
  gsm <- sl(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  ttl <- sl(base::grep("^!Sample_title", hd, value = TRUE)[1])[-1]

  fs <- list.files(file.path(RAW_DIR, g), pattern = "\\.(txt|csv)\\.gz$", full.names = TRUE)
  fs <- fs[file.info(fs)$size > 5 * 1024^2]
  if (!length(fs)) { say("  [FAIL] no matrix containing data found"); next }
  f <- fs[base::which.max(file.info(fs)$size)]
  say("  Reading ", base::basename(f), " (", sprintf("%.2f GB", file.info(f)$size / 1024^3), ")")
  t0 <- Sys.time()
  dt <- read_matrix(f)
  say("    time ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))),
      " min; ", base::nrow(dt), " rows x ", base::ncol(dt), " columns")

  probes <- base::as.character(dt[[1]])
  cn <- base::names(dt)[-1]
  M <- base::as.matrix(dt[, -1, with = FALSE])
  rm(dt); invisible(gc())

  # ---- paired columns: by suffix first, then by repeated names ----
  has_suf <- base::mean(base::grepl(BETA_SUF, cn) | base::grepl(PVAL_SUF, cn))
  if (has_suf > 0.8) {
    keep <- base::grepl(BETA_SUF, cn)
    say("  Column names carry beta / detection P value suffixes; beta columns kept: ", base::sum(keep), "")
    M <- M[, keep, drop = FALSE]
    cn <- base::sub(BETA_SUF, "", cn[keep])
  } else if (base::mean(base::duplicated(cn)) > 0.4) {
    n2 <- length(cn) / 2
    inter <- base::all(cn[base::seq(1, length(cn), 2)] == cn[base::seq(2, length(cn), 2)])
    block <- length(cn) %% 2 == 0 && base::all(cn[1:n2] == cn[(n2 + 1):length(cn)])
    say("  Column names repeat in pairs: alternating ", inter, ", in blocks ", block)
    if (!inter && !block) { say("  [FAIL] layout cannot be determined"); next }
    i1 <- if (inter) base::seq(1, length(cn), 2) else 1:n2
    i2 <- if (inter) base::seq(2, length(cn), 2) else (n2 + 1):length(cn)
    k <- base::min(2000L, base::nrow(M))
    m1 <- base::colMeans(M[1:k, i1, drop = FALSE], na.rm = TRUE)
    m2 <- base::colMeans(M[1:k, i2, drop = FALSE], na.rm = TRUE)
    first_is_beta <- base::mean(m1 > m2, na.rm = TRUE) > 0.5
    say("  The group with the larger mean is taken as beta: ", if (first_is_beta) "first group" else "second group",
        " (median of means in the two groups ", sprintf("%.3f / %.3f", stats::median(m1, na.rm = TRUE),
                                   stats::median(m2, na.rm = TRUE)), ")")
    sel <- if (first_is_beta) i1 else i2
    M <- M[, sel, drop = FALSE]; cn <- cn[sel]
  } else say("  Single-column matrix")

  # ---- sample mapping ----
  mp <- map_samples(cn, gsm, ttl, hd)
  say("  Mapping strategy: ", mp$name, ", success rate ", sprintf("%.3f", mp$frac))
  if (mp$frac < MAP_MIN) {
    say("  [FAIL] mapping success rate below ", MAP_MIN, ". Column name examples: ",
        base::paste(utils::head(cn, 3), collapse = " | "))
    if (!is.null(mp$score)) {
      sc <- mp$score[base::order(-mp$score$frac), , drop = FALSE]
      say("  Success rate of each strategy (top 6): ")
      for (i in seq_len(min(6L, base::nrow(sc))))
        say(sprintf("    %-28s %.3f", sc$name[i], sc$frac[i]))
    }
    next
  }
  gsm_of_col <- gsm[base::match(cn, mp$key)]
  base::colnames(M) <- gsm_of_col
  keep_col <- base::which(!base::is.na(gsm_of_col) & gsm_of_col %in% ph$gsm)
  M <- M[, keep_col, drop = FALSE]
  ph <- ph[base::match(base::colnames(M), ph$gsm), , drop = FALSE]
  say("  Analysis samples obtained: ", base::ncol(M), " (controls ", base::sum(ph$group == "control"),
      ", cases ", base::sum(ph$group == "case"), ")")

  # ---- EPIC v2 probe names ----
  if (C$platform == "EPICv2") {
    bs <- base::sub("_.*$", "", probes)
    nd <- base::sum(base::duplicated(bs))
    say("  Suffixes stripped: duplicated probes ", nd, " rows, ",
        length(base::unique(bs[base::duplicated(bs)])), " sites, merged by the mean")
    if (nd > 0) {
      s <- base::rowsum(M, group = bs, na.rm = TRUE)
      n <- base::rowsum((!base::is.na(M)) * 1, group = bs, na.rm = TRUE)
      M <- s / n; M[!base::is.finite(M)] <- NA
    } else base::rownames(M) <- bs
  } else base::rownames(M) <- probes

  # ---- probe filtering ----
  n0 <- base::nrow(M)
  M <- M[base::rownames(M) %in% UNIVERSE, , drop = FALSE]
  say("  Intersection with the probe universe: ", n0, " -> ", base::nrow(M))
  miss <- base::rowMeans(base::is.na(M))
  M <- M[miss <= MISS_CUT, , drop = FALSE]
  say("  Missing rate above ", MISS_CUT, ": removed ", base::sum(miss > MISS_CUT), ", remaining ", base::nrow(M))

  q <- stats::quantile(M, c(0.01, 0.5, 0.99), na.rm = TRUE)
  say("  beta quantiles 1%/50%/99%: ", base::paste(sprintf("%.3f", q), collapse = " / "))
  if (q[[2]] < 0.2 || q[[2]] > 0.9) say("  [NOTE] the median is outside the usual range; please confirm manually")

  base::saveRDS(list(gse = g, disease = C$disease, platform = C$platform,
               beta = M, pheno = ph), out_rds)
  say("  Saved (", sprintf("%.2f GB", file.size(out_rds) / 1024^3), ")")
  summ[[length(summ) + 1]] <- data.frame(
    gse = g, disease = C$disease, platform = C$platform,
    n_sample = base::ncol(M), n_probe = base::nrow(M),
    n_control = base::sum(ph$group == "control"),
    n_case = base::sum(ph$group == "case"), stringsAsFactors = FALSE)
  rm(M); invisible(gc())
  say("")
}

if (length(summ)) {
  tab <- base::do.call(rbind, summ)
  write.csv(tab, file.path(OUT_DIR, "build_summary.csv"), row.names = FALSE)
  say("=== Summary ==="); print(tab)
}
log_stop()
