# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 08c: install maxprobes + test each package in actual use
#
# The result of 08b had an ordering problem: maxprobes was attempted before minfi and the 450K annotation
# package were installed, so it failed; those two packages were then pulled in as dependencies of
# EpiSmokEr. The environment is now different and a simple retry is enough.
#
# In addition, "available" in the 08b summary only means that requireNamespace finds the package,
# not that it works. This script makes one real call to each package.
#
# Environment: Windows + RStudio  R 4.4.1
# Location: <project root>/scripts/08c_verify_and_retry.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
OUT_DIR <- file.path(PROJ_ROOT, "output", "08_setup")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

options(timeout = 1200)
log <- character(0)
say <- function(...) { m <- paste0(...); message(m); log <<- c(log, m) }
have <- function(p) requireNamespace(p, quietly = TRUE)

check <- function(label, expr) {
  say("---- ", label, " ----")
  r <- tryCatch({ force(expr); "pass" },
                error = function(e) paste0("failed: ", conditionMessage(e)))
  say("  ", r); invisible(r)
}

ANNO <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"

# --------------------------------------------------------------
say("=== Installing maxprobes (minfi and the annotation package are now in place) ===")
if (!have("maxprobes")) {
  say("  minfi available: ", have("minfi"), "; annotation package available: ", have(ANNO))
  r <- tryCatch({
    remotes::install_github("markgene/maxprobes", upgrade = "never", quiet = FALSE)
    "Installation command executed"
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  say("  ", r)
}
say("  maxprobes now available: ", have("maxprobes"))
say("")

# --------------------------------------------------------------
say("=== Testing each package ===")

check("EpiDISH reference matrix and RPC call", {
  e <- new.env()
  data("centDHSbloodDMC.m", package = "EpiDISH", envir = e)
  ref <- get("centDHSbloodDMC.m", envir = e)
  say("  Reference matrix ", nrow(ref), " probes x ", ncol(ref), " types")
  say("  Cell types: ", paste(colnames(ref), collapse = ", "))
  set.seed(1)
  toy <- matrix(runif(nrow(ref) * 4), nrow = nrow(ref),
                dimnames = list(rownames(ref), paste0("s", 1:4)))
  est <- EpiDISH::epidish(beta.m = toy, ref.m = ref, method = "RPC")$estF
  say("  Dimensions of test output ", nrow(est), " x ", ncol(est),
      ", row sums about ", sprintf("%.3f", mean(rowSums(est))))
  say("  CellDMC exists: ", "CellDMC" %in% getNamespaceExports("EpiDISH"))
})

check("450K annotation", {
  a <- minfi::getAnnotation(getFromNamespace(ANNO, ANNO))
  say("  Annotation: ", nrow(a), " rows")
  say("  chrX ", sum(a$chr == "chrX"), ", chrY ", sum(a$chr == "chrY"))
  sc <- intersect(c("Probe_rs", "CpG_rs", "SBE_rs", "Probe_maf"), colnames(a))
  say("  SNP-related columns: ", paste(sc, collapse = ", "))
  say("  Probes with Probe_rs: ", sum(!is.na(a$Probe_rs)))
  say("  Gene column UCSC_RefGene_Name present: ",
      "UCSC_RefGene_Name" %in% colnames(a))
})

check("maxprobes cross-reactive probe lists", {
  if (!have("maxprobes")) stop("maxprobes is not installed")
  x450 <- maxprobes::xreactive_probes(array_type = "450K")
  say("  450K cross-reactive probes: ", length(x450))
  say("  EPIC ", length(maxprobes::xreactive_probes(array_type = "EPIC")))
})

check("EpiSmokEr weight table", {
  d <- utils::data(package = "EpiSmokEr")$results[, "Item"]
  say("  Data sets in the package: ", paste(d, collapse = ", "))
  say("  Exported functions: ", paste(getNamespaceExports("EpiSmokEr"), collapse = ", "))
})

check("dnaMethyAge clock list", {
  ac <- dnaMethyAge::availableClock()
  say("  Available clocks: ", length(ac), "", paste(ac, collapse = ", "))
})

# --------------------------------------------------------------
say("")
say("=== Summary ===")
for (p in c("EpiDISH", ANNO, "minfi", "maxprobes", "EpiSmokEr",
            "dnaMethyAge", "limma", "data.table")) {
  say(sprintf("  %-46s %s", p, if (have(p)) "installed" else "not installed"))
}
say("")
say("Library paths: ", .libPaths()[1])
if (have("renv")) say("Note: this is an renv project; remember renv::snapshot() to record versions after installing")

writeLines(log, file.path(OUT_DIR, "08c_verify_report.txt"), useBytes = TRUE)
message("\nReport: ", file.path(OUT_DIR, "08c_verify_report.txt"))
message("Next: review 08c_verify_report.txt.")
