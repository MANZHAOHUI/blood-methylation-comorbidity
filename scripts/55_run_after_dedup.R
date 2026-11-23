# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 55: after de-duplication of GSE87648, re-run all downstream steps in order
#
# Run 54 first to replace the model of GSE87648 with the de-duplicated version, then re-run the downstream steps. Every step that reads *_design.rds or
# *_ewas_final.rds must be re-run so that the result tables agree with the new model. This step is only
# a scheduler that sources the step scripts one by one in a fixed order; each step runs in its own new environment
# and memory is released afterwards; if any step fails the scheduler stops.
#
# Order and dependencies:
#   54 de-duplicate GSE87648 and rebuild the model (safe to repeat)   about 2 min
#   39 concordance matrix (45 pairs, 1,000 permutations)              about 160 min
#   40 leukocyte composition axis                                     about 11 min
#   41 upstream audit (sex inference, Fig. S4a)                       about 5 min
#   42 sex mismatch in GSE87648                                       about 1 min
#   43 bimodality quality index (now reads the exclusion list of 54)  about 2 min
#   44 quality covariate sensitivity                                  about 60 min
#   45 within-disease positive controls (10,000 permutations)         about 190 min
#   46 detection limits (depends on 39 and 45)                        about 1 min
#   47 region level (2,000 permutations)                              about 90 min
#   50 overlap approach                                               about 2 min
#   51 sex scale (45 pairs, 500 permutations)                         about 90 min
#   52 diagnosis of GSE198904 (reference cohorts include GSE87648)    about 13 min
#   53 manuscript figures (produced only if the self-checks pass)     about 5 min
#   About 10.5 hours in total; running overnight is recommended.
#
# Environment: Windows + RStudio. Can be interrupted at any time (Esc in the RStudio console or Stop);
# sourcing the script again continues from the unfinished step; finished steps are skipped on the basis of their logs.
# Location: <project root>/scripts/55_run_after_dedup.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
SCRIPT_NAME <- "55_run_after_dedup"
SCR_DIR <- file.path(PROJ_ROOT, "scripts")
OUT_DIR <- file.path(PROJ_ROOT, "output", "55_run")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(OUT_DIR)) stop("Failed to create the output directory: ", OUT_DIR)

# The scheduler itself opens no log file (each step script writes its own formal log under output/logs);
# it only prints progress to the console and writes the run time of each step to output/55_run/run_timing.csv.
say <- function(...) cat(..., "\n", sep = "")

STEPS <- c("54_dedup_gse87648", "39_concordance_matrix", "40_composition_matrix", "41_upstream_audit",
           "42_sex_mismatch_gse87648", "43_bimodality_qc", "44_quality_covariate_sensitivity",
           "45_positive_controls_stable", "46_detection_limits_all_pairs",
           "47_region_level_robustness", "50_overlap_screen", "51_sex_positive_control",
           "52_gse198904_diagnosis", "53_manuscript_figures")


# Each step is sourced within this session in its own new environment; gc() releases large matrices between steps.
# (Rscript child processes are not used: the scripts run in an renv project and a child process would not find the renv library.)
# Interruptible and resumable: a step with a log that was written after the exclusion list of step 54 and ends with "Status: completed normally"
# is considered finished and skipped; an interrupted step has no such log and is re-run.
LOG_DIR <- file.path(PROJ_ROOT, "output", "logs")
EX_FILE <- file.path(PROJ_ROOT, "output", "54_dedup", "excluded_samples.csv")
step_done <- function(s) {
  if (!file.exists(EX_FILE)) return(FALSE)
  t_ex <- file.mtime(EX_FILE)
  logs <- list.files(LOG_DIR, pattern = base::paste0("^", s, "_[0-9]{8}_[0-9]{6}\\.log$"), full.names = TRUE)
  if (!length(logs)) return(FALSE)
  logs <- logs[file.mtime(logs) >= t_ex - 60]
  if (!length(logs)) return(FALSE)
  for (lg in logs) {
    tx <- readLines(lg, warn = FALSE, encoding = "UTF-8")
    if (base::any(base::grepl("Status: completed normally", tx, fixed = TRUE))) return(TRUE)
  }
  FALSE
}

timing <- list()
for (s in STEPS) {
  f <- file.path(SCR_DIR, base::paste0(s, ".R"))
  if (!file.exists(f)) stop("Script not found: ", f)
  if (s != "54_dedup_gse87648" && step_done(s)) {
    say("[SKIPPED] ", s, ": a complete log from after de-duplication exists")
    timing[[s]] <- data.frame(step = s, ok = TRUE, minutes = 0, stringsAsFactors = FALSE)
    next
  }
  say("================ Running ", s, "  ", format(Sys.time(), "%H:%M:%S"), " ================")
  t0 <- Sys.time()
  env <- new.env(parent = globalenv())
  ok <- tryCatch({ source(f, local = env, encoding = "UTF-8", echo = FALSE); TRUE },
                 error = function(e) { message("Step ", s, " failed with an error: ", conditionMessage(e)); FALSE })
  rm(env); invisible(gc())
  # When a sub-script fails, its log and sink may still be open; they are closed here so that the console resumes output
  if (!ok && exists(".LOGENV") && !is.null(.LOGENV$con)) try(log_stop(), silent = TRUE)
  mins <- base::as.numeric(difftime(Sys.time(), t0, units = "mins"))
  timing[[s]] <- data.frame(step = s, ok = ok, minutes = round(mins, 1), stringsAsFactors = FALSE)
  if (!ok) {
    say("[FAILED] ", s, ", elapsed ", sprintf("%.1f", mins), " min; stopping. See the log of this step under output/logs")
    write.csv(base::do.call(rbind, timing), file.path(OUT_DIR, "run_timing.csv"), row.names = FALSE)
    stop("Step ", s, " failed")
  }
  say("[DONE] ", s, "  elapsed ", sprintf("%.1f", mins), " min")
  say("")
  if (s == "54_dedup_gse87648") {
    ex <- file.path(PROJ_ROOT, "output", "54_dedup", "excluded_samples.csv")
    if (!file.exists(ex)) stop("Step 54 did not write an exclusion list")
    say("Exclusion list of step 54: ")
    print(read.csv(ex, stringsAsFactors = FALSE)[, c("gse", "gsm", "kept_partner", "group_dropped")])
    print(read.csv(file.path(PROJ_ROOT, "output", "54_dedup", "before_after_summary.csv"), stringsAsFactors = FALSE))
    say("")
  }
}
TM <- base::do.call(rbind, timing)
write.csv(TM, file.path(OUT_DIR, "run_timing.csv"), row.names = FALSE)
print(TM)
say("")
say("All ", length(STEPS), " steps finished, ", sprintf("%.1f", base::sum(TM$minutes) / 60), " hours in total.")
