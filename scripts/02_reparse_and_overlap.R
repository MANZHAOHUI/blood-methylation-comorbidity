# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders
# Step 02: re-parse the metadata + check overlap between the PD cohorts
#
# Part A: re-parse the sample characteristics by tag name
#         Step 01 assumed that the !Sample_characteristics_ch1 lines are aligned across samples.
#         That assumption does not hold for GSE72774, GSE87648, GSE152027 and others,
#         so columns were shifted (e.g. 8637 in the age column, pack-years in the smoking column).
#         This part reads "tag: value" sample by sample and assigns values by tag,
#         checks ages for plausibility and counts discarded values without a tag.
#
# Part B: use sample titles to check whether GSE111629 / GSE72774 / GSE145361 share individuals
#         (zero cost, done first)
#
# Part C: the 450K array carries 65 rs probes; use them as a genotype fingerprint to confirm overlap
#         (the remote beta matrix is read as a stream, not saved to disk and not loaded whole into memory)
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/02_reparse_and_overlap.R
# Prerequisite: 01 has been run; headers are cached in output/01_metadata/*_header.txt
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
if (!dir.exists(PROJ_ROOT)) stop("Project root does not exist: ", PROJ_ROOT)

HDR_DIR <- file.path(PROJ_ROOT, "output", "01_metadata")
if (!dir.exists(HDR_DIR)) {
  stop("Output directory of 01 not found; please run 01_extract_geo_metadata.R first: ", HDR_DIR)
}

OUT_DIR <- file.path(PROJ_ROOT, "output", "02_overlap")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

options(timeout = 3600)

GSE_LIST <- c("GSE145361", "GSE111629", "GSE72774",
              "GSE87648",  "GSE112611",
              "GSE152027", "GSE84727",  "GSE80417")

PD_SETS <- c("GSE145361", "GSE111629", "GSE72774")

# Switch for part C: the GSE145361 matrix is very large; running the first two sets first is recommended
RUN_FINGERPRINT   <- TRUE
FINGERPRINT_SETS  <- c("GSE111629", "GSE72774")   # add "GSE145361" when needed

split_fields <- function(line) {
  gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]])
}

read_header <- function(gse) {
  f <- file.path(HDR_DIR, paste0(gse, "_header.txt"))
  if (!file.exists(f)) stop("Cached header missing: ", f)
  readLines(f, warn = FALSE, encoding = "UTF-8")
}

# ==============================================================
# Part A: assign by tag
# ==============================================================
parse_by_tag <- function(lines, gse) {
  acc_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
  if (length(acc_line) == 0L) stop(gse, " header lacks !Sample_geo_accession")
  acc <- split_fields(acc_line[1])[-1]
  n <- length(acc)

  df <- data.frame(geo_accession = acc, stringsAsFactors = FALSE)

  t_line <- grep("^!Sample_title", lines, value = TRUE)
  if (length(t_line)) {
    v <- split_fields(t_line[1])[-1]
    df$sample_title <- if (length(v) == n) v else NA_character_
  }

  ch <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  tag_all <- character(0); val_all <- character(0); idx_all <- integer(0)
  n_untagged <- 0L

  for (l in ch) {
    v <- split_fields(l)[-1]
    if (length(v) != n) v <- c(v, rep(NA_character_, n))[seq_len(n)]
    has <- !is.na(v) & nzchar(v) & grepl(":", v, fixed = TRUE)
    n_untagged <- n_untagged + sum(!is.na(v) & nzchar(v) & !has)
    if (!any(has)) next
    tag_all <- c(tag_all, trimws(sub(":.*$", "", v[has])))
    val_all <- c(val_all, trimws(sub("^[^:]*:", "", v[has])))
    idx_all <- c(idx_all, which(has))
  }

  if (length(tag_all)) {
    tag_all <- make.names(tolower(trimws(tag_all)))
    for (tg in unique(tag_all)) {
      sel <- which(tag_all == tg)
      v <- rep(NA_character_, n)
      dup <- 0L
      for (k in sel) {
        i <- idx_all[k]
        if (is.na(v[i])) v[i] <- val_all[k] else dup <- dup + 1L
      }
      if (dup > 0L) {
        message("  [", gse, "] field ", tg, ": ", dup,
                " repeated values within one sample; the first occurrence was kept")
      }
      df[[tg]] <- v
    }
  }
  attr(df, "n_untagged") <- n_untagged
  df
}

summarise_df <- function(df, gse) {
  out <- c(paste0("### ", gse, "   samples = ", nrow(df)))
  nu <- attr(df, "n_untagged")
  if (!is.null(nu) && nu > 0L) {
    out <- c(out, paste0("  [warn] characteristic values without a tag: ", nu, "; discarded, the header needs manual inspection"))
  }
  for (nm in names(df)) {
    if (nm %in% c("geo_accession", "sample_title")) next
    v <- df[[nm]]
    num <- suppressWarnings(as.numeric(v))
    num_ok <- sum(!is.na(num))
    if (grepl("age", nm, ignore.case = TRUE) && num_ok > 0.5 * length(v)) {
      bad <- sum(!is.na(num) & (num < 0 | num > 110))
      out <- c(out, sprintf("  [age] %s: n=%d  min=%.1f  median=%.1f  max=%.1f%s",
                            nm, num_ok, min(num, na.rm = TRUE),
                            median(num, na.rm = TRUE), max(num, na.rm = TRUE),
                            if (bad > 0) paste0("   <<< ", bad, " implausible values; this column is suspect") else ""))
      next
    }
    u <- unique(v[!is.na(v)])
    if (length(u) <= 15L) {
      tb <- sort(table(v, useNA = "ifany"), decreasing = TRUE)
      out <- c(out, paste0("  [", nm, "] ",
                           paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                                 collapse = "; ")))
    } else if (num_ok > 0.5 * length(v)) {
      out <- c(out, sprintf("  [num] %s: n=%d  min=%.3g  median=%.3g  max=%.3g",
                            nm, num_ok, min(num, na.rm = TRUE),
                            median(num, na.rm = TRUE), max(num, na.rm = TRUE)))
    } else {
      out <- c(out, sprintf("  [text] %s: %d distinct values (first 5: %s)",
                            nm, length(u), paste(head(u, 5), collapse = ", ")))
    }
  }
  c(out, "")
}

meta <- list()
report <- character(0)
for (gse in GSE_LIST) {
  message("\n==== Part A  ", gse, " ====")
  df <- parse_by_tag(read_header(gse), gse)
  meta[[gse]] <- df
  write.csv(df, file.path(OUT_DIR, paste0(gse, "_metadata_bytag.csv")),
            row.names = FALSE)
  report <- c(report, summarise_df(df, gse))
}
writeLines(report, file.path(OUT_DIR, "02A_metadata_report.txt"), useBytes = TRUE)
message("\nPart A finished: ", file.path(OUT_DIR, "02A_metadata_report.txt"))

# ==============================================================
# Part B: overlap of sample titles
# ==============================================================
message("\n==== Part B: comparison of sample titles ====")
titles <- lapply(PD_SETS, function(g) {
  d <- meta[[g]]
  if (is.null(d$sample_title)) character(0) else trimws(d$sample_title)
})
names(titles) <- PD_SETS

b_lines <- "### Part B: overlap of sample titles among the three PD data sets"
for (g in PD_SETS) {
  b_lines <- c(b_lines, sprintf("  %s title examples: %s", g,
                                paste(head(titles[[g]], 5), collapse = " | ")))
}
for (i in 1:(length(PD_SETS) - 1)) {
  for (j in (i + 1):length(PD_SETS)) {
    a <- titles[[PD_SETS[i]]]; b <- titles[[PD_SETS[j]]]
    inter <- intersect(a[nzchar(a)], b[nzchar(b)])
    b_lines <- c(b_lines, sprintf("  %s vs %s : %d identical titles (denominators %d / %d)",
                                  PD_SETS[i], PD_SETS[j], length(inter),
                                  length(a), length(b)))
    if (length(inter)) {
      b_lines <- c(b_lines, paste0("    e.g. ", paste(head(inter, 8), collapse = ", ")))
    }
  }
}
writeLines(b_lines, file.path(OUT_DIR, "02B_title_overlap.txt"), useBytes = TRUE)
cat(paste(b_lines, collapse = "\n"), "\n")

# ==============================================================
# Part C: genotype fingerprint from rs probes
# ==============================================================
SUPP <- list(
  GSE145361 = list(
    u = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145361/suppl/GSE145361_Vallerga2020_NCOMMS_AvgBeta_Matrix-file.txt.gz",
    sep = "\t"),
  GSE111629 = list(
    u = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111629/suppl/GSE111629_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz",
    sep = "\t"),
  GSE72774 = list(
    u = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72774/suppl/GSE72774_datBetaNormalized.csv.gz",
    sep = ",")
)

# Scan the remote matrix as a stream and keep only rows starting with rs
extract_rs <- function(gse) {
  cfg <- SUPP[[gse]]
  if (is.null(cfg)) stop("No matrix URL configured for ", gse, "")
  message("  Streaming ", gse, "... (the whole file is traversed; the time depends on the network speed)")
  con <- gzcon(url(cfg$u, open = "rb"))
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  hdr <- readLines(con, n = 1L, warn = FALSE)
  if (length(hdr) == 0L) stop(gse, ": matrix is empty or the URL is invalid")
  samp <- gsub('^"|"$', "", strsplit(hdr, cfg$sep, fixed = TRUE)[[1]])
  samp <- samp[-1]

  keep <- character(0); nline <- 0L
  repeat {
    chunk <- readLines(con, n = 20000L, warn = FALSE)
    if (length(chunk) == 0L) break
    nline <- nline + length(chunk)
    hit <- chunk[grepl('^"?rs', chunk)]
    if (length(hit)) keep <- c(keep, hit)
  }
  message("    scanned ", nline, " rows, rs probes found: ", length(keep), "")
  if (length(keep) == 0L) {
    warning(gse, ": the matrix contains no rs probes (filtered upstream); the fingerprint method cannot be used for this set")
    return(NULL)
  }

  parts <- strsplit(keep, cfg$sep, fixed = TRUE)
  ids <- gsub('^"|"$', "", vapply(parts, `[`, character(1), 1L))
  mat <- t(vapply(parts, function(p) {
    suppressWarnings(as.numeric(gsub('"', "", p[-1])))
  }, numeric(length(samp))))
  rownames(mat) <- ids
  colnames(mat) <- samp
  mat
}

geno_call <- function(beta) {
  g <- matrix(NA_integer_, nrow(beta), ncol(beta),
              dimnames = dimnames(beta))
  g[beta < 0.25] <- 0L
  g[beta >= 0.25 & beta <= 0.75] <- 1L
  g[beta > 0.75] <- 2L
  g
}

if (RUN_FINGERPRINT && length(FINGERPRINT_SETS) >= 2L) {
  message("\n==== Part C: rs-probe fingerprint comparison ====")
  geno <- list()
  for (g in FINGERPRINT_SETS) {
    b <- extract_rs(g)
    if (!is.null(b)) geno[[g]] <- geno_call(b)
  }
  c_lines <- "### Part C: genotype concordance of rs probes"
  gs <- names(geno)
  if (length(gs) >= 2L) {
    for (i in 1:(length(gs) - 1)) {
      for (j in (i + 1):length(gs)) {
        A <- geno[[gs[i]]]; B <- geno[[gs[j]]]
        common <- intersect(rownames(A), rownames(B))
        c_lines <- c(c_lines, sprintf("  %s vs %s : %d shared rs probes",
                                      gs[i], gs[j], length(common)))
        if (length(common) < 20L) {
          c_lines <- c(c_lines, "    too few shared probes; the result is unreliable")
          next
        }
        A <- A[common, , drop = FALSE]; B <- B[common, , drop = FALSE]
        conc <- matrix(NA_real_, ncol(A), ncol(B),
                       dimnames = list(colnames(A), colnames(B)))
        for (cc in seq_len(ncol(A))) {
          conc[cc, ] <- colMeans(B == A[, cc], na.rm = TRUE)
        }
        hit <- which(conc >= 0.90, arr.ind = TRUE)
        c_lines <- c(c_lines,
                     sprintf("    sample pairs with concordance >= 0.90: %d", nrow(hit)),
                     sprintf("    concordance distribution: median %.3f, 99th percentile %.3f, maximum %.3f",
                             median(conc, na.rm = TRUE),
                             quantile(conc, 0.99, na.rm = TRUE),
                             max(conc, na.rm = TRUE)))
        if (nrow(hit)) {
          pairs_df <- data.frame(
            sample_a = rownames(conc)[hit[, 1]],
            sample_b = colnames(conc)[hit[, 2]],
            concordance = conc[hit],
            stringsAsFactors = FALSE)
          pairs_df <- pairs_df[order(-pairs_df$concordance), ]
          write.csv(pairs_df,
                    file.path(OUT_DIR, sprintf("02C_matched_pairs_%s_%s.csv",
                                               gs[i], gs[j])),
                    row.names = FALSE)
          c_lines <- c(c_lines, "    matched pairs written to CSV")
        }
      }
    }
  }
  writeLines(c_lines, file.path(OUT_DIR, "02C_fingerprint.txt"), useBytes = TRUE)
  cat(paste(c_lines, collapse = "\n"), "\n")
} else {
  message("\nPart C skipped (RUN_FINGERPRINT = FALSE)")
}

message("\nAll output is in: ", OUT_DIR)
message("Next: review 02A_metadata_report.txt, 02B_title_overlap.txt and 02C_fingerprint.txt.")

