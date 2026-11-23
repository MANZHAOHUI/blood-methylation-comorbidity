# ==============================================================
# Shared component: reading GEO matrices and matching samples
#
# Sourced by script 34 and later. The reading logic that already works in 33 is factored out here
# so that each script does not carry its own copy - six data sets come in four file layouts, and repeated implementations would drift.
#
# Cases handled (all observed in practice, none assumed):
#   delimiter      tab or comma, decided from the first line
#   shifted header the header has one field fewer than the data rows (the probe column is unnamed)
#   series matrix  skip to the line after !series_matrix_table_begin first
#   paired columns beta and detection P value side by side; may carry suffixes, alternate, or come in blocks
#   sample match   scan the raw value, the part after the colon and the barcode within it for every !Sample_ field
#
# All functions call masked generics with the base:: prefix (BiocGenerics masks
# grep, grepl, match, order, paste and others; unqualified calls reach the wrong object).
# ==============================================================

MR_BARCODE  <- "[0-9]{8,}_R[0-9]{2}C[0-9]{2}"
MR_BETA_SUF <- "\\.(Noob_)?Beta$|\\.AVG_Beta$|_Beta$"
MR_PVAL_SUF <- "\\.Detection_Pval$|\\.Detection\\.Pval$|\\.pval$|_Pval$"

mr_split <- function(x, sep = "\t")
  base::trimws(base::gsub('^"|"$', "", base::strsplit(x, sep, fixed = TRUE)[[1]]))

mr_header <- function(raw_dir, g, n = 2000L) {
  f <- file.path(raw_dir, g, base::paste0(g, "_series_matrix.txt.gz"))
  if (!file.exists(f)) stop("Header file not found: ", f)
  con <- gzfile(f, "rt"); on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

mr_data_file <- function(raw_dir, g, min_mb = 5) {
  fs <- list.files(file.path(raw_dir, g), pattern = "\\.(txt|csv)\\.gz$", full.names = TRUE)
  fs <- fs[file.info(fs)$size > min_mb * 1024^2]
  if (!length(fs)) stop("No matrix containing data found: ", g)
  fs[base::which.max(file.info(fs)$size)]
}

# Read the matrix and return a data.table (first column = probe ID)
mr_read <- function(f, verbose = TRUE) {
  is_series <- base::grepl("series_matrix", base::basename(f), fixed = TRUE)
  n_skip <- 0L
  if (is_series) {
    con <- gzfile(f, "rt")
    k <- 0L
    repeat {
      ln <- readLines(con, n = 1L, warn = FALSE)
      if (!length(ln)) { close(con); stop("Data start marker not found: ", base::basename(f)) }
      k <- k + 1L
      if (base::grepl("!series_matrix_table_begin", ln, fixed = TRUE)) break
    }
    close(con)
    n_skip <- k
  }
  con <- gzfile(f, "rt")
  head_lines <- readLines(con, n = n_skip + 2L, warn = FALSE)
  close(con)
  two <- head_lines[(n_skip + 1L):(n_skip + 2L)]
  sep <- if (base::grepl("\t", two[1], fixed = TRUE)) "\t" else ","
  h <- mr_split(two[1], sep); d <- mr_split(two[2], sep)
  shifted <- length(h) == length(d) - 1L
  if (verbose) message(sprintf("    delimiter %s; header %d fields, data rows %d fields%s",
                               if (sep == "\t") "tab" else "comma",
                               length(h), length(d),
                               if (shifted) "; header lacks the first field, treated as shifted" else ""))
  dt <- data.table::fread(f, skip = n_skip + 1L, header = FALSE, sep = sep,
                          check.names = FALSE, showProgress = FALSE, fill = TRUE)
  if (is_series) {
    last <- base::as.character(dt[[1]][base::nrow(dt)])
    if (!base::is.na(last) && base::startsWith(last, "!")) dt <- dt[-base::nrow(dt), ]
  }
  nm <- if (shifted) c("probe_id", h) else h
  if (length(nm) != base::ncol(dt))
    stop("Number of column names ", length(nm), " does not match number of columns ", base::ncol(dt), "")
  data.table::setnames(dt, nm)
  dt
}

# Paired-column decision: return the indices of the columns kept and the matching sample names
mr_pick_beta <- function(M, cn, verbose = TRUE) {
  has_suf <- base::mean(base::grepl(MR_BETA_SUF, cn) | base::grepl(MR_PVAL_SUF, cn))
  if (has_suf > 0.8) {
    keep <- base::grepl(MR_BETA_SUF, cn)
    if (verbose) message("  Column names carry beta / detection P value suffixes; beta columns kept: ", base::sum(keep), "")
    return(list(idx = base::which(keep), names = base::sub(MR_BETA_SUF, "", cn[keep])))
  }
  if (base::mean(base::duplicated(cn)) <= 0.4) {
    if (verbose) message("  Single-column matrix")
    return(list(idx = seq_along(cn), names = cn))
  }
  n2 <- length(cn) / 2
  inter <- base::all(cn[base::seq(1, length(cn), 2)] == cn[base::seq(2, length(cn), 2)])
  block <- length(cn) %% 2 == 0 && base::all(cn[1:n2] == cn[(n2 + 1):length(cn)])
  if (verbose) message("  Column names repeat in pairs: alternating ", inter, ", in blocks ", block)
  if (!inter && !block) stop("Layout of paired columns cannot be determined")
  i1 <- if (inter) base::seq(1, length(cn), 2) else 1:n2
  i2 <- if (inter) base::seq(2, length(cn), 2) else (n2 + 1):length(cn)
  k <- base::min(2000L, base::nrow(M))
  m1 <- base::colMeans(M[1:k, i1, drop = FALSE], na.rm = TRUE)
  m2 <- base::colMeans(M[1:k, i2, drop = FALSE], na.rm = TRUE)
  first <- base::mean(m1 > m2, na.rm = TRUE) > 0.5
  if (verbose) message("  The group with the larger mean is taken as beta: ", if (first) "first group" else "second group",
                       sprintf(" (median of means in the two groups %.3f / %.3f)",
                               stats::median(m1, na.rm = TRUE),
                               stats::median(m2, na.rm = TRUE)))
  sel <- if (first) i1 else i2
  list(idx = sel, names = cn[sel])
}

# Sample matching: scan all !Sample_ fields and return the best strategy
mr_map_samples <- function(cols, hd) {
  gsm <- mr_split(base::grep("^!Sample_geo_accession", hd, value = TRUE)[1])[-1]
  ttl <- mr_split(base::grep("^!Sample_title", hd, value = TRUE)[1])[-1]
  cand <- list(list(name = "GSM accession", key = gsm),
               list(name = "sample title", key = ttl))
  bc_ttl <- base::rep(NA_character_, length(ttl))
  hit <- base::regexpr(MR_BARCODE, ttl) > 0
  if (base::any(hit)) {
    bc_ttl[hit] <- base::regmatches(ttl, base::regexpr(MR_BARCODE, ttl))
    cand[[length(cand) + 1]] <- list(name = "barcode in title", key = bc_ttl)
  }
  cand[[length(cand) + 1]] <- list(
    name = "sample_N in title",
    key = base::sub("^\\s*sample[ _]?([0-9]+).*$", "sample_\\1", ttl, ignore.case = TRUE))
  for (t in base::unique(base::sub("\t.*$", "", base::grep("^!Sample_", hd, value = TRUE)))) {
    for (one in base::grep(base::paste0("^", t, "\t"), hd, value = TRUE)) {
      v <- mr_split(one)[-1]
      if (length(v) != length(gsm)) next
      nm <- base::sub("^!Sample_", "", t)
      cand[[length(cand) + 1]] <- list(name = base::paste0(nm, " raw value"), key = v)
      vv <- base::trimws(base::sub("^[^:]*:", "", v))
      if (!identical(vv, v))
        cand[[length(cand) + 1]] <- list(name = base::paste0(nm, " after colon"), key = vv)
      h2 <- base::regexpr(MR_BARCODE, v) > 0
      if (base::any(h2)) {
        k <- base::rep(NA_character_, length(v))
        k[h2] <- base::regmatches(v, base::regexpr(MR_BARCODE, v))
        cand[[length(cand) + 1]] <- list(name = base::paste0(nm, " barcode within"), key = k)
      }
    }
  }
  best <- list(frac = -1); sc <- list()
  for (one in cand) {
    fr <- base::mean(cols %in% one$key)
    sc[[length(sc) + 1]] <- data.frame(name = one$name, frac = fr, stringsAsFactors = FALSE)
    if (fr > best$frac) best <- list(frac = fr, name = one$name, key = one$key)
  }
  best$gsm <- gsm
  best$score <- base::do.call(rbind, sc)
  best
}
