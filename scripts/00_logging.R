# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Shared logging component
#
# Usage: three lines at the top of every script
#   source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
#   log_start("09_probe_filtering")
#   ...script body; all output goes through say(), not message()...
#   log_stop()
#
# Recorded: the complete console output, warnings, errors with the call stack, run time,
# versions of R and of key packages, and the path and run arguments of the script itself.
# Logs are written to <project root>/output/logs/<name>_<timestamp>.log.
#
# Two design points:
#   (1) Standard output uses sink(split = TRUE), so it goes to both console and log,
#       and the progress of long jobs remains visible in real time.
#   (2) When a script stops on an error, the options(error=) hook writes the call stack to the log
#       and closes all sinks. Otherwise the sink is left dangling and the console stays silent,
#       which looks as if R were broken.
# ==============================================================

.LOGENV <- new.env(parent = emptyenv())

log_start <- function(name,
                      proj_root = get0("PROJ_ROOT", ifnotfound =
                        "C:/blood_methylation_project"),
                      capture_messages = TRUE) {
  if (!is.null(.LOGENV$con)) log_stop(quiet = TRUE)

  dir <- file.path(proj_root, "output", "logs")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  if (!dir.exists(dir)) stop("Failed to create log directory: ", dir)

  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  path <- file.path(dir, paste0(name, "_", stamp, ".log"))
  con <- file(path, open = "wt", encoding = "UTF-8")

  .LOGENV$con <- con
  .LOGENV$path <- path
  .LOGENV$name <- name
  .LOGENV$t0 <- Sys.time()
  .LOGENV$msg_sunk <- isTRUE(capture_messages)
  .LOGENV$old_error <- getOption("error")

  message("Log: ", path)

  # Clear warnings left over from the previous script. warnings() reads the session-level last.warning,
  # so when several scripts run in one R session, old warnings would be attributed to the new script.
  try(suppressWarnings(rm("last.warning", envir = baseenv())), silent = TRUE)
  assign("last.warning", NULL, envir = globalenv())
  try(rm("last.warning", envir = globalenv()), silent = TRUE)

  sink(con, split = TRUE)
  if (isTRUE(capture_messages)) sink(con, type = "message")

  cat("==============================================\n")
  cat("Script: ", name, "\n", sep = "")
  cat("Start: ", format(.LOGENV$t0, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("R    : ", R.version.string, "  ", R.version$platform, "\n", sep = "")
  cat("Working directory: ", getwd(), "\n", sep = "")
  cat("Library paths: ", .libPaths()[1], "\n", sep = "")
  base_fns <- ls(baseenv())
  shadow <- character(0)
  for (env in setdiff(search(), c(".GlobalEnv", "package:base"))) {
    if (!grepl("^package:", env)) next
    hit <- intersect(ls(env), base_fns)
    if (length(hit)) shadow <- c(shadow, paste0(env, " -> ",
                                                paste(hit, collapse = ", ")))
  }
  if (length(shadow)) {
    cat("Loaded packages that mask base functions (calling these function names unqualified reaches the wrong object):\n")
    for (x in shadow) cat("  ", x, "\n", sep = "")
  } else {
    cat("No package masks base functions\n")
  }
  cat("==============================================\n\n")

  options(error = function() {
    cat("\n---------- Script stopped on an error ----------\n")
    cat("Error message: ", geterrmessage(), sep = "")
    # traceback() cannot be used here: it reads the call stack of the previous error; the current one is only written after the hook returns.
    # sys.calls() gives the call stack of the current error; the last frame is this hook itself and is dropped.
    tb <- tryCatch({
      sc <- sys.calls()
      sc <- sc[-length(sc)]
      if (!length(sc)) "(top-level call, no call stack)" else {
        ln <- vapply(sc, function(cl) {
          tx <- paste(deparse(cl, width.cutoff = 200L, nlines = 2L), collapse = " ")
          sr <- attr(cl, "srcref")
          if (!is.null(sr)) tx <- paste0(tx, " at ", basename(attr(sr, "srcfile")$filename), "#", sr[1L])
          if (nchar(tx) > 300L) tx <- paste0(substr(tx, 1L, 300L), " ...")
          tx
        }, character(1))
        paste(sprintf("%d: %s", rev(seq_along(ln)), rev(ln)), collapse = "\n")
      }
    }, error = function(e) "(call stack not available)")
    cat("Call stack:\n", tb, "\n", sep = "")
    log_stop(status = "stopped on error")
  })

  invisible(path)
}

# Output goes through say(): with the sink open it reaches both console and log
say <- function(...) {
  cat(paste0(...), "\n", sep = "")
  invisible(NULL)
}

log_stop <- function(status = "completed normally", quiet = FALSE) {
  if (is.null(.LOGENV$con)) return(invisible(NULL))

  el <- difftime(Sys.time(), .LOGENV$t0, units = "mins")
  cat("\n==============================================\n")
  cat("Status: ", status, "\n", sep = "")
  cat("End: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "   elapsed ", sprintf("%.1f", as.numeric(el)), " min\n", sep = "")

  w <- tryCatch(warnings(), error = function(e) NULL)
  if (length(w)) {
    cat("Warnings: ", length(w), "\n", sep = "")
    nm <- names(w)
    for (i in seq_len(min(30L, length(nm)))) cat("  - ", nm[i], "\n", sep = "")
    if (length(nm) > 30L) cat("  (remaining ", length(nm) - 30L, " omitted)\n", sep = "")
  } else {
    cat("Warnings: none\n")
  }

  pk <- c("EpiDISH", "minfi", "limma", "data.table", "EpiSmokEr",
          "dnaMethyAge", "maxprobes",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  cat("Versions of key packages:\n")
  for (p in pk) {
    v <- tryCatch(as.character(utils::packageVersion(p)),
                  error = function(e) NA_character_)
    cat("  ", p, "  ", if (is.na(v)) "not installed" else v, "\n", sep = "")
  }
  cat("==============================================\n")

  # Close the sinks level by level; the order must not be reversed
  if (isTRUE(.LOGENV$msg_sunk)) try(sink(type = "message"), silent = TRUE)
  while (sink.number() > 0L) try(sink(), silent = TRUE)
  try(close(.LOGENV$con), silent = TRUE)

  options(error = .LOGENV$old_error)
  path <- .LOGENV$path
  .LOGENV$con <- NULL

  if (!quiet) message("Log written: ", path)
  invisible(path)
}
