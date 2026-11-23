# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 34: covariates of the six additional cohorts
#
# Identical to what step 11 did for the first six sets; three items:
#   (1) twelve leukocyte fractions from EpiDISH (reference matrix cent12CT450k.m)
#   (2) DNAm age (dnaMethyAge, HorvathS2013)
#   (3) DNAm smoking score (SSc of EpiSmokEr)
#
# The key convention of 11 is kept: always use the raw matrix without probe selection.
# Some of the sites needed by the clock and the reference matrix are among the probes removed by D2
# (sex chromosomes, SNPs, cross-reactive), so a filtered matrix would lack sites.
# The raw matrices of the first round are stored under rds/; those of this round are re-read from the archives in raw/,
# using the shared component 00_matrix_reader.R rather than another copy of the reading logic.
#
# EpiSmokEr is not attached: it exports a function named all, which would mask base::all,
# unqualified all() in later scripts of the same session would reach the wrong object. All calls use EpiSmokEr::.
#
# Sets are processed and saved one by one; a failure midway does not affect cohorts already finished.
#
# Environment: Windows + RStudio. Peak memory occurs with GSE338167 (850 thousand probes x 621 samples)
# Location: <project root>/scripts/34_covariates_panel2.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "34_covariates_panel2"

RAW_DIR <- file.path(DATA_ROOT, "raw")
SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
COV_DIR <- file.path(DATA_ROOT, "covariates")
OUT_DIR <- file.path(PROJ_ROOT, "output", "34_covariates")
for (d in c(COV_DIR, OUT_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(data.table))

FORCE <- TRUE      # Age and smoking score were not obtained in the previous run and must be recomputed, hence TRUE
GSES <- c("GSE84727", "GSE80417", "GSE198904", "GSE273970", "GSE338167", "GSE144858")
EPICV2 <- "GSE338167"

# --------------------------------------------------------------
# The internal data of dnaMethyAge and EpiSmokEr must be loaded first
#
# The functions of both packages depend on their lazily loaded data objects (probeAnnotation21kdatMethUsed
# of Horvath, Illig_data of EpiSmokEr, and others). When functions are called
# with the package name and double colon, these objects are not loaded automatically, which gives
# "object ... not found".
#
# library() is not a solution: EpiSmokEr exports a function named all,
# which masks base::all once attached, so unqualified all() in later scripts of the same session fails.
#
# Instead, all data sets of both packages are loaded into the global environment. The chain of environments in which package functions look up variables
# ends at the global environment, so the objects are found without attaching any package.
# --------------------------------------------------------------
ensure_pkg_data <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) { say("  [NOTE] package not found: ", pkg); return(character(0)) }
  ds <- tryCatch(utils::data(package = pkg)$results[, "Item"], error = function(e) character(0))
  ds <- base::unique(base::sub("\\s.*$", "", ds))
  ok <- character(0)
  for (d in ds) {
    r <- tryCatch({ utils::data(list = d, package = pkg, envir = globalenv()); TRUE },
                  error = function(e) FALSE, warning = function(w) FALSE)
    if (isTRUE(r) && exists(d, envir = globalenv())) ok <- c(ok, d)
  }
  say("  ", pkg, ": data sets ", length(ds), " (",
      base::paste(ds, collapse = ", "), ")")
  say("    loaded successfully: ", length(ok),
      if (length(ok)) base::paste0(" (", base::paste(utils::head(ok, 8), collapse = ", "),
                                   if (length(ok) > 8) " ..." else "", ")") else "")
  ok
}
say("=== Loading the internal data of the two packages ===")
# For EpiSmokEr the global-environment route is verified (Illig_data in place, smoking scores obtained for all six sets)
loaded_smk <- ensure_pkg_data("EpiSmokEr")

# It does not work for dnaMethyAge: the function calls data("HorvathS2013") internally, and data() without the package
# argument only searches attached packages, so the package must be attached.
# EpiSmokEr must not be attached (its exported all masks base::all); dnaMethyAge has no such problem.
say("  dnaMethyAge needs to be attached (it fetches clock data internally with data() without a package name)")
before <- search()
att <- tryCatch({ library(dnaMethyAge); TRUE },
                error = function(e) { say("  [FAIL] attach failed: ", conditionMessage(e)); FALSE })
if (att) {
  say("  Attached; new search path entries: ",
      base::paste(base::setdiff(search(), before), collapse = ", "))
  cf <- base::intersect(base::ls("package:dnaMethyAge"), base::ls("package:base"))
  say("  Exported objects with the same name as base objects: ", if (length(cf)) base::paste(cf, collapse = ", ") else "none")
}
for (need in c("probeAnnotation21kdatMethUsed", "Illig_data")) {
  say("  Key object ", need, ": ",
      if (exists(need)) "in place" else "[not in place, but can be fetched by the package itself if attach succeeded]")
}
say("")

e <- new.env()
utils::data("cent12CT450k.m", package = "EpiDISH", envir = e)
REF <- base::get("cent12CT450k.m", envir = e)
say("EpiDISH reference matrix: ", base::nrow(REF), " probes x ", base::ncol(REF), " types: ",
    base::paste(base::colnames(REF), collapse = ", "))
say("")

summ <- list()
for (g in GSES) {
  out_csv <- file.path(COV_DIR, base::paste0(g, "_covariates.csv"))
  say("==== ", g, " ====")
  if (file.exists(out_csv) && !FORCE) { say("  Exists, skipped (set FORCE to rebuild)"); next }

  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  want <- base::colnames(obj$beta)
  ph <- obj$pheno
  rm(obj); invisible(gc())
  say("  Target samples: ", length(want), " samples")

  hd <- mr_header(RAW_DIR, g)
  f <- mr_data_file(RAW_DIR, g)
  say("  Reading raw matrix ", base::basename(f),
      " (", sprintf("%.2f GB", file.info(f)$size / 1024^3), ")")
  t0 <- Sys.time()
  dt <- mr_read(f)
  probes <- base::as.character(dt[[1]])
  cn <- base::names(dt)[-1]
  M <- base::as.matrix(dt[, -1, with = FALSE])
  rm(dt); invisible(gc())
  pick <- mr_pick_beta(M, cn)
  M <- M[, pick$idx, drop = FALSE]
  cn <- pick$names
  mp <- mr_map_samples(cn, hd)
  say("  Mapping strategy: ", mp$name, ", success rate ", sprintf("%.3f", mp$frac))
  if (mp$frac < 0.95) { say("  [FAIL] mapping failed"); next }
  base::colnames(M) <- mp$gsm[base::match(cn, mp$key)]

  if (g %in% EPICV2) {
    bs <- base::sub("_.*$", "", probes)
    nd <- base::sum(base::duplicated(bs))
    say("  EPIC v2: suffixes stripped and duplicated probes merged: ", nd, " records")
    if (nd > 0) {
      s <- base::rowsum(M, group = bs, na.rm = TRUE)
      n <- base::rowsum((!base::is.na(M)) * 1, group = bs, na.rm = TRUE)
      M <- s / n; M[!base::is.finite(M)] <- NA
      rm(s, n); invisible(gc())
    } else base::rownames(M) <- bs
  } else base::rownames(M) <- probes

  keep <- base::which(base::colnames(M) %in% want)
  M <- M[, keep, drop = FALSE]
  say("  Obtained ", base::ncol(M), " samples, probes ", base::nrow(M),
      "; time ", sprintf("%.1f", base::as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  if (base::ncol(M) < length(want))
    say("  [NOTE] fewer than the target number of samples, difference ", length(want) - base::ncol(M), " samples")

  res <- data.frame(gsm = base::colnames(M), stringsAsFactors = FALSE)

  # ---- cell fractions ----
  ov <- base::intersect(base::rownames(REF), base::rownames(M))
  say("  Reference matrix probes found: ", length(ov), " / ", base::nrow(REF))
  if (length(ov) < 0.8 * base::nrow(REF)) {
    say("  [NOTE] low hit rate; fraction estimates may be unstable")
  }
  est <- EpiDISH::epidish(beta.m = M, ref.m = REF, method = "RPC")$estF
  est <- est[base::match(res$gsm, base::rownames(est)), , drop = FALSE]
  base::colnames(est) <- base::paste0("ct_", base::colnames(est))
  res <- cbind(res, base::as.data.frame(est))
  say("  Median row sum ", sprintf("%.3f", stats::median(base::rowSums(est), na.rm = TRUE)),
      "; median neutrophil fraction ", sprintf("%.3f", stats::median(est[, "ct_Neu"], na.rm = TRUE)))

  # ---- DNAm age ----
  ma <- tryCatch(dnaMethyAge::methyAge(M, clock = "HorvathS2013", do_plot = FALSE),
                 error = function(err) { say("  [NOTE] age computation failed: ", conditionMessage(err)); NULL })
  if (!is.null(ma)) {
    col <- base::intersect(c("mAge", "Age", "HorvathS2013"), base::names(ma))[1]
    idc <- base::intersect(c("Sample", "sample", "ID"), base::names(ma))[1]
    res$dnam_age_horvath <- base::as.numeric(ma[[col]])[base::match(res$gsm, ma[[idc]])]
    say("  DNAm age: median ", sprintf("%.1f", stats::median(res$dnam_age_horvath, na.rm = TRUE)),
        ", range ", sprintf("%.1f to %.1f", base::min(res$dnam_age_horvath, na.rm = TRUE),
                           base::max(res$dnam_age_horvath, na.rm = TRUE)))
    if ("age_clean" %in% base::names(ph)) {
      a <- base::suppressWarnings(base::as.numeric(ph$age_clean))[base::match(res$gsm, ph$gsm)]
      ok <- base::is.finite(a) & base::is.finite(res$dnam_age_horvath)
      if (base::sum(ok) > 30)
        say("  Correlation with chronological age: r = ", sprintf("%.3f", stats::cor(a[ok], res$dnam_age_horvath[ok])),
            ", mean difference ", sprintf("%.1f", base::mean(res$dnam_age_horvath[ok] - a[ok])), " years")
    }
  } else res$dnam_age_horvath <- NA_real_

  # ---- smoking score ----
  ss <- tryCatch(EpiSmokEr::SSc(dataset = M),
                 error = function(err) { say("  [NOTE] smoking score failed: ", conditionMessage(err)); NULL })
  if (!is.null(ss)) {
    sv <- if (is.data.frame(ss)) ss[[base::ncol(ss)]] else base::as.numeric(ss)
    nmv <- if (is.data.frame(ss)) base::rownames(ss) else base::names(ss)
    res$dnam_smoking <- base::as.numeric(sv)[base::match(res$gsm, nmv)]
    say("  Smoking score: median ", sprintf("%.3f", stats::median(res$dnam_smoking, na.rm = TRUE)),
        ", missing ", base::sum(base::is.na(res$dnam_smoking)), " samples")
  } else res$dnam_smoking <- NA_real_

  write.csv(res, out_csv, row.names = FALSE)
  say("  Saved ", base::basename(out_csv))
  summ[[length(summ) + 1]] <- data.frame(
    gse = g, n = base::nrow(res),
    ref_hit = length(ov),
    age_ok = base::sum(!base::is.na(res$dnam_age_horvath)),
    smoke_ok = base::sum(!base::is.na(res$dnam_smoking)),
    stringsAsFactors = FALSE)
  rm(M); invisible(gc())
  say("")
}

if (length(summ)) {
  tab <- base::do.call(rbind, summ)
  write.csv(tab, file.path(OUT_DIR, "covariate_summary.csv"), row.names = FALSE)
  say("=== Summary ==="); print(tab)
}
log_stop()
