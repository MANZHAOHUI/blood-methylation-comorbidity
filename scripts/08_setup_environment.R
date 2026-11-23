# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 08: set up and verify the analysis environment
#
# Only packages that are actually needed are installed. Each is tested minimally right after installation,
# and a list at the end states which succeeded, which failed and which step a failure affects.
#
# Already available: BiocManager, data.table, limma
# Needed:
#   EpiDISH                              cell fraction estimation + CellDMC (no substitute)
#   IlluminaHumanMethylation450kanno...  probe annotation: chromosome, gene, SNP flags (D2 and functional annotation)
#   maxprobes                            list of cross-reactive probes (Chen 2013 / Pidsley 2016)
#   EpiSmokEr                            DNAm smoking score (SSc of Elliott 2014, continuous)
#   methylclock or dnaMethyAge           DNAm age (60% of ages are missing in GSE152027, a proxy is needed)
#
# Not installed for now: wateRmelon (no normalisation from IDAT), sva (batch strategy undecided)
#
# Environment: Windows + RStudio  R 4.4.1
# Location: <project root>/scripts/08_setup_environment.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
OUT_DIR <- file.path(PROJ_ROOT, "output", "08_setup")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

# Mirror in China, clearly faster when international links are slow. Comment out these two lines if the mirror is unavailable.
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options(timeout = 1200)

log <- character(0)
say <- function(...) { m <- paste0(...); message(m); log <<- c(log, m) }

have <- function(p) requireNamespace(p, quietly = TRUE)

need_cran <- function(p) {
  if (have(p)) return(invisible(TRUE))
  say("  Installing CRAN package ", p)
  try(install.packages(p), silent = TRUE)
  have(p)
}

need_bioc <- function(p) {
  if (have(p)) return(invisible(TRUE))
  say("  Installing Bioconductor package ", p, " (many dependencies, may take a while)")
  try(BiocManager::install(p, ask = FALSE, update = FALSE), silent = TRUE)
  have(p)
}

need_github <- function(repo, p) {
  if (have(p)) return(invisible(TRUE))
  if (!have("remotes")) need_cran("remotes")
  if (!have("remotes")) { say("  remotes could not be installed; cannot install ", p); return(FALSE) }
  say("  Installing from GitHub: ", repo)
  try(remotes::install_github(repo, upgrade = "never"), silent = TRUE)
  have(p)
}

say("R version: ", R.version.string)
say("")

status <- list()

# --------------------------------------------------------------
say("[1] EpiDISH - cell fractions and CellDMC")
status$EpiDISH <- isTRUE(need_bioc("EpiDISH"))
if (status$EpiDISH) {
  ok <- try({
    data("centDHSbloodDMC.m", package = "EpiDISH", envir = environment())
    ref <- get("centDHSbloodDMC.m", envir = environment())
    say("  Reference matrix centDHSbloodDMC.m: ", nrow(ref), " probes x ", ncol(ref), " cell types")
    say("  Cell types: ", paste(colnames(ref), collapse = ", "))
    TRUE
  }, silent = TRUE)
  status$EpiDISH <- !inherits(ok, "try-error")
}
say("  Result: ", if (isTRUE(status$EpiDISH)) "available" else "not available")
say("")

# --------------------------------------------------------------
say("[2] 450K probe annotation")
ANNO <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
status$anno <- isTRUE(need_bioc(ANNO))
if (status$anno) {
  ok <- try({
    a <- minfi::getAnnotation(get(ANNO, envir = asNamespace(ANNO)))
    say("  Annotation rows: ", nrow(a), "; columns: ", paste(head(colnames(a), 12), collapse = ", "))
    say("  Sex-chromosome probes: chrX ", sum(a$chr == "chrX"), ", chrY ", sum(a$chr == "chrY"))
    snp_cols <- intersect(c("Probe_rs", "CpG_rs", "SBE_rs"), colnames(a))
    say("  Available SNP flag columns: ", paste(snp_cols, collapse = ", "))
    TRUE
  }, silent = TRUE)
  status$anno <- !inherits(ok, "try-error")
}
say("  Result: ", if (isTRUE(status$anno)) "available" else "not available")
say("")

# --------------------------------------------------------------
say("[3] maxprobes - cross-reactive probes")
status$maxprobes <- isTRUE(need_github("markgene/maxprobes", "maxprobes"))
if (status$maxprobes) {
  ok <- try({
    x450 <- maxprobes::xreactive_probes(array_type = "450K")
    xepic <- maxprobes::xreactive_probes(array_type = "EPIC")
    say("  450K cross-reactive probes: ", length(x450), "; EPIC ", length(xepic))
    TRUE
  }, silent = TRUE)
  status$maxprobes <- !inherits(ok, "try-error")
}
say("  Result: ", if (isTRUE(status$maxprobes)) "available" else "not available")
say("")

# --------------------------------------------------------------
say("[4] EpiSmokEr - DNAm smoking score")
status$EpiSmokEr <- isTRUE(need_github("sailalithabollepalli/EpiSmokEr", "EpiSmokEr"))
if (status$EpiSmokEr) {
  ok <- try({
    e <- new.env()
    data("Illig_data", package = "EpiSmokEr", envir = e)
    obj <- get(ls(e)[1], envir = e)
    say("  SSc weight table: ", nrow(obj), " rows, columns: ", paste(colnames(obj), collapse = ", "))
    TRUE
  }, silent = TRUE)
  if (inherits(ok, "try-error")) {
    say("  The weight table has an unexpected name; listing the data sets of the package for checking: ")
    say("  ", paste(utils::data(package = "EpiSmokEr")$results[, "Item"], collapse = ", "))
  }
}
say("  Result: ", if (isTRUE(status$EpiSmokEr)) "installed" else "not available")
say("")

# --------------------------------------------------------------
say("[5] DNAm age - try methylclock from Bioconductor first, then dnaMethyAge")
status$clock <- FALSE
clock_pkg <- NA_character_
if (need_bioc("methylclock")) { status$clock <- TRUE; clock_pkg <- "methylclock" }
if (!status$clock) {
  if (need_github("yiluyucheng/dnaMethyAge", "dnaMethyAge")) {
    status$clock <- TRUE; clock_pkg <- "dnaMethyAge"
  }
}
say("  Result: ", if (status$clock) paste0("available (", clock_pkg, ")") else "neither could be installed")
say("")

# --------------------------------------------------------------
say("========== Summary ==========")
impact <- c(
  EpiDISH    = "cell fractions and CellDMC; without it no cell-type attribution is possible",
  anno       = "probe filtering (sex chromosomes, SNP probes) and gene annotation",
  maxprobes  = "filtering of cross-reactive probes; without it the step can only be declared as not done in the limitations",
  EpiSmokEr  = "DNAm smoking score; without it smoking cannot be adjusted for in cohorts lacking a smoking field",
  clock      = "DNAm age; without it the missing ages of GSE152027 (60%) cannot be filled"
)
for (k in names(impact)) {
  say(sprintf("  %-10s %-6s  %s", k,
              if (isTRUE(status[[k]])) "available" else "not available", impact[[k]]))
}
say("")
say("sessionInfo summary: ")
si <- utils::sessionInfo()
say("  ", si$R.version$version.string, "  Platform ", si$platform)

writeLines(log, file.path(OUT_DIR, "08_setup_report.txt"), useBytes = TRUE)
message("\nReport: ", file.path(OUT_DIR, "08_setup_report.txt"))
message("Next: review 08_setup_report.txt.")
