# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Step 08b: diagnosis of failed installations
#
# Step 08 wrapped the installations in try(silent = TRUE); failures were only reported as "not available"
# and the original error text was lost - exactly the kind of silent failure that should be avoided.
# This script swallows no errors: every step prints the full error text and writes it to the report.
#
# Clue: dnaMethyAge was installed successfully from GitHub, so remotes works,
# GitHub is reachable and the build toolchain is fine. The failure of the other two GitHub packages
# is therefore most likely caused by the Bioconductor packages they depend on, and all Bioconductor
# installations failed as well. Bioconductor is therefore diagnosed first.
#
# Environment: Windows + RStudio  R 4.4.1
# Location: <project root>/scripts/08b_diagnose_install.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
OUT_DIR <- file.path(PROJ_ROOT, "output", "08_setup")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create directory: ", OUT_DIR)

options(timeout = 1200)

log <- character(0)
say <- function(...) { m <- paste0(...); message(m); log <<- c(log, m) }

# Evaluate the expression and keep the complete information whether it succeeds or not
run <- function(label, expr) {
  say("---- ", label, " ----")
  res <- tryCatch(withCallingHandlers(
      { force(expr); "OK" },
      warning = function(w) {
        say("  [WARNING] ", conditionMessage(w)); invokeRestart("muffleWarning")
      }),
    error = function(e) paste0("ERROR: ", conditionMessage(e)))
  say("  Result: ", res)
  invisible(res)
}

have <- function(p) requireNamespace(p, quietly = TRUE)

# --------------------------------------------------------------
say("=== Basic information ===")
say("R: ", R.version.string)
say("Platform: ", R.version$platform)
say("Library paths: ", paste(.libPaths(), collapse = " | "))
say("Library path writable: ",
    paste(vapply(.libPaths(), function(p) as.character(file.access(p, 2) == 0),
                 character(1)), collapse = " | "))
say("Current repos: ", paste(getOption("repos"), collapse = " | "))
say("BiocManager installed: ", have("BiocManager"))
if (have("BiocManager")) {
  say("BiocManager version: ", as.character(utils::packageVersion("BiocManager")))
  run("BiocManager::version()", say("  Bioc version: ", as.character(BiocManager::version())))
  run("BiocManager::repositories()",
      say("  ", paste(BiocManager::repositories(), collapse = " | ")))
}
say("")

# --------------------------------------------------------------
say("=== Test 1: default Bioconductor repository (no mirror) ===")
options(BioC_mirror = NULL)
options(repos = c(CRAN = "https://cloud.r-project.org"))
run("Reading the Bioc software package index", {
  db <- utils::available.packages(repos = BiocManager::repositories()["BioCsoft"])
  say("  Index entries: ", nrow(db))
  say("  Contains EpiDISH: ", "EpiDISH" %in% rownames(db))
})
run("Reading the Bioc annotation package index", {
  db2 <- utils::available.packages(repos = BiocManager::repositories()["BioCann"])
  say("  Index entries: ", nrow(db2))
  say("  Contains the 450K annotation package: ",
      "IlluminaHumanMethylation450kanno.ilmn12.hg19" %in% rownames(db2))
})
say("")

say("=== Test 2: install EpiDISH (default repository) ===")
if (!have("EpiDISH")) run("BiocManager::install(EpiDISH)",
                          BiocManager::install("EpiDISH", ask = FALSE, update = FALSE))
say("  EpiDISH now available: ", have("EpiDISH"))
say("")

say("=== Test 3: install the 450K annotation package (default repository) ===")
ANNO <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
if (!have(ANNO)) run(paste0("BiocManager::install(", ANNO, ")"),
                     BiocManager::install(ANNO, ask = FALSE, update = FALSE))
say("  Annotation package now available: ", have(ANNO))
say("")

# --------------------------------------------------------------
say("=== Test 4: GitHub channel (dnaMethyAge succeeded, so the channel itself works) ===")
say("  remotes installed: ", have("remotes"))
if (!have("EpiDISH")) {
  run("remotes::install_github(sjczheng/EpiDISH)",
      remotes::install_github("sjczheng/EpiDISH", upgrade = "never"))
  say("  EpiDISH now available: ", have("EpiDISH"))
}
if (!have("maxprobes")) {
  run("remotes::install_github(markgene/maxprobes)",
      remotes::install_github("markgene/maxprobes", upgrade = "never"))
  say("  maxprobes now available: ", have("maxprobes"))
}
if (!have("EpiSmokEr")) {
  run("remotes::install_github(sailalithabollepalli/EpiSmokEr)",
      remotes::install_github("sailalithabollepalli/EpiSmokEr", upgrade = "never"))
  say("  EpiSmokEr now available: ", have("EpiSmokEr"))
}
say("")

# --------------------------------------------------------------
say("=== Summary of the current state ===")
for (p in c("EpiDISH", ANNO, "maxprobes", "EpiSmokEr", "dnaMethyAge",
            "minfi", "limma", "data.table")) {
  say(sprintf("  %-46s %s", p, if (have(p)) "available" else "not available"))
}

# List of clocks available in dnaMethyAge; decides how DNAm age is computed
if (have("dnaMethyAge")) {
  say("")
  say("=== Clocks provided by dnaMethyAge ===")
  run("Listing available clocks", {
    nm <- ls(asNamespace("dnaMethyAge"))
    say("  Exported functions: ", paste(getNamespaceExports("dnaMethyAge"), collapse = ", "))
  })
}

writeLines(log, file.path(OUT_DIR, "08b_diagnose_report.txt"), useBytes = TRUE)
message("\nReport: ", file.path(OUT_DIR, "08b_diagnose_report.txt"))
message("Next: review the complete 08b_diagnose_report.txt; the original error text matters.")
