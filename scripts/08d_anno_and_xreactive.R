# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 08d: fix the annotation package call + obtain the cross-reactive probe lists
#
# Two things:
# (1) The annotation package test in 08c failed because of a wrong call, not because of the package.
#     Annotation packages of this kind must be attached with library() first; the annotation is then taken from the object of the same name.
# (2) maxprobes still cannot be installed, and remotes hides the real error in the console.
#     Its whole value is two lists of probe IDs, which is not worth a prolonged struggle with the installation:
#     this script tries to obtain the real error and at the same time downloads the list files directly from the repository;
#     either route is sufficient.
#
# Environment: Windows + RStudio  R 4.4.1
# Location: <project root>/scripts/08d_anno_and_xreactive.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
OUT_DIR <- file.path(PROJ_ROOT, "output", "08_setup")
REF_DIR <- file.path(DATA_ROOT, "probe_lists")
for (d in c(OUT_DIR, REF_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create directory: ", d)
}

options(timeout = 1200)

# Console output, warnings, errors and package versions are all recorded in a log file under output/logs/
source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start("08d_anno_and_xreactive", proj_root = PROJ_ROOT)

have <- function(p) requireNamespace(p, quietly = TRUE)
ANNO <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"

# --------------------------------------------------------------
say("=== 1. 450K annotation (using the correct way of attaching) ===")
r <- tryCatch({
  suppressPackageStartupMessages(library(ANNO, character.only = TRUE))
  a <- minfi::getAnnotation(get(ANNO))
  say("  Annotation: ", nrow(a), " rows, ", ncol(a), " columns")
  say("  chrX ", sum(a$chr == "chrX"), ", chrY ", sum(a$chr == "chrY"))
  sc <- intersect(c("Probe_rs", "Probe_maf", "CpG_rs", "CpG_maf",
                    "SBE_rs", "SBE_maf"), colnames(a))
  say("  SNP-related columns: ", paste(sc, collapse = ", "))
  say("  Probe_rs non-empty: ", sum(!is.na(a$Probe_rs)),
      "; CpG_rs non-empty: ", sum(!is.na(a$CpG_rs)),
      "; SBE_rs non-empty: ", sum(!is.na(a$SBE_rs)))
  say("  UCSC_RefGene_Name present: ", "UCSC_RefGene_Name" %in% colnames(a))
  saveRDS(a, file.path(REF_DIR, "anno450k.rds"))
  say("  Saved anno450k.rds")
  "pass"
}, error = function(e) paste0("failed: ", conditionMessage(e)))
say("  ", r)
say("")

# --------------------------------------------------------------
say("=== 2. Reference matrices shipped with EpiDISH ===")
r <- tryCatch({
  d <- utils::data(package = "EpiDISH")$results[, "Item"]
  say("  Data sets: ", paste(d, collapse = ", "))
  for (nm in d) {
    e <- new.env(); utils::data(list = nm, package = "EpiDISH", envir = e)
    o <- get(ls(e)[1], envir = e)
    if (is.matrix(o)) say("    ", nm, ": ", nrow(o), " x ", ncol(o),
                          "  Column names: ", paste(colnames(o), collapse = ", "))
  }
  "pass"
}, error = function(e) paste0("failed: ", conditionMessage(e)))
say("  ", r)
say("")

# --------------------------------------------------------------
say("=== 3. maxprobes: obtain the real installation error ===")
if (!have("maxprobes")) {
  tb <- file.path(tempdir(), "maxprobes.tar.gz")
  r <- tryCatch({
    utils::download.file("https://api.github.com/repos/markgene/maxprobes/tarball/master",
                         tb, mode = "wb", quiet = TRUE)
    say("  tarball downloaded, ", round(file.size(tb) / 1024), " KB")
    out <- system2(file.path(R.home("bin"), "R"),
                   c("CMD", "INSTALL", shQuote(tb)),
                   stdout = TRUE, stderr = TRUE)
    say("  --- raw output of R CMD INSTALL ---")
    for (l in out) say("   | ", l)
    "executed"
  }, error = function(e) paste0("failed: ", conditionMessage(e)))
  say("  ", r)
}
say("  maxprobes now available: ", have("maxprobes"))
say("")

# --------------------------------------------------------------
say("=== 4. Bypass the installation and fetch the cross-reactive probe lists directly ===")
r <- tryCatch({
  api <- "https://api.github.com/repos/markgene/maxprobes/contents/inst/extdata"
  js <- paste(readLines(url(api), warn = FALSE), collapse = "")
  nm <- regmatches(js, gregexpr('"name":"[^"]+"', js))[[1]]
  nm <- sub('^"name":"', "", sub('"$', "", nm))
  dl <- regmatches(js, gregexpr('"download_url":"[^"]+"', js))[[1]]
  dl <- gsub("\\\\/", "/", sub('^"download_url":"', "", sub('"$', "", dl)))
  say("  Files in inst/extdata of the repository: ", length(nm), "")
  for (i in seq_along(nm)) say("    ", nm[i])
  for (i in seq_along(nm)) {
    dest <- file.path(REF_DIR, nm[i])
    if (!file.exists(dest)) utils::download.file(dl[i], dest, mode = "wb", quiet = TRUE)
    say("    ", nm[i], "  ", round(file.size(dest) / 1024, 1), " KB")
    hd <- tryCatch(readLines(dest, n = 3, warn = FALSE), error = function(e) character(0))
    if (length(hd)) say("      first lines: ", paste(substr(hd, 1, 90), collapse = " // "))
  }
  "pass"
}, error = function(e) paste0("failed: ", conditionMessage(e)))
say("  ", r)
say("")

say("=== Summary ===")
say("  Annotation package available: ", have(ANNO))
say("  maxprobes available: ", have("maxprobes"))
say("  Directory of downloaded list files: ", REF_DIR)
say("  Directory contents: ", paste(list.files(REF_DIR), collapse = ", "))

log_stop()
message("The complete console log is in ", file.path(PROJ_ROOT, "output", "logs"))
message("Next: review that log file.")
