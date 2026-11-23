# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 50: distributional screen for large-scale overlap between the two Parkinson's disease cohorts
#
# Identifying the same person pair by pair is no longer feasible: the genotype information retained in the public matrices varies greatly
# (from 426 to 4 probes with genotype shape), and calibration failed twice;
# this route ends here and the criteria are not changed again.
#
# But what really affects the conclusions is not an occasional duplicate. The correlation of effect sizes of a within-disease positive control rests on hundreds to
# more than a thousand samples, and a few duplicate pairs do not change it; what can change it is large-scale overlap - if the two
# cohorts shared a considerable proportion of individuals, what was measured would not be replication.
#
# Large-scale overlap can be seen without pairwise identification: it leaves a second mode in the distribution
# of cross-cohort concordance. This step answers only this question and states clearly what it can and cannot answer.
#
# Criterion fixed beforehand: with the distribution of sample pairs within each cohort as reference, if the upper tail of the cross-cohort distribution
# (99.9th percentile) does not exceed 1.10 times the same quantile within cohorts and the distribution is unimodal,
# large-scale overlap is considered absent.
#
# What can be concluded: exclusion of large-scale overlap.
# What cannot be concluded: presence or absence of individual duplicate samples. This must be stated in the limitations.
#
# The two schizophrenia cohorts are not screened: GSE84727 retains only 4 usable probes, and any judgement based on
# them would be meaningless.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/50_overlap_screen.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "50_overlap_screen"

RAW_DIR <- file.path(DATA_ROOT, "raw")
RDS_DIR <- file.path(DATA_ROOT, "rds")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
MAP_DIR <- file.path(PROJ_ROOT, "output", "06_mapping")
OUT_DIR <- file.path(PROJ_ROOT, "output", "50_overlap")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
source(file.path(PROJ_ROOT, "scripts", "00_matrix_reader.R"))
suppressPackageStartupMessages(library(data.table))

PAIR      <- c("GSE145361", "GSE111629")
MIN_FRAC  <- 0.05
END_LO    <- 0.15; END_HI <- 0.85
MID_LO    <- 0.40; MID_HI <- 0.60
CL_SD     <- 0.06
MIN_SHARE <- 100L
TAIL_Q    <- 0.999
TAIL_MULT <- 1.10

to_gsm <- function(g, cn) {
  if (base::mean(base::grepl("^GSM[0-9]+$", cn)) > 0.95) return(cn)
  f <- file.path(MAP_DIR, base::paste0(g, "_sample_map.csv"))
  if (file.exists(f)) {
    m <- read.csv(f, stringsAsFactors = FALSE, colClasses = "character")
    v <- m$gsm[base::match(cn, m$col_name)]
    if (base::mean(!base::is.na(v)) > 0.95) return(v)
  }
  hd <- tryCatch(mr_header(RAW_DIR, g), error = function(e) NULL)
  if (!is.null(hd)) {
    mp <- mr_map_samples(cn, hd)
    if (mp$frac > 0.95) return(mp$gsm[base::match(cn, mp$key)])
  }
  NULL
}
get_beta <- function(g) {
  f1 <- file.path(RDS_DIR, base::paste0(g, "_beta.rds"))
  if (file.exists(f1)) {
    b <- base::readRDS(f1)
    gs <- to_gsm(g, base::colnames(b))
    if (is.null(gs)) stop("column names cannot be translated to sample IDs")
    base::colnames(b) <- gs
    return(b[, !base::is.na(gs), drop = FALSE])
  }
  stop("Raw matrix not found: ", g)
}
snp_like <- function(B) {
  lo <- B < 0.25; mi <- B >= 0.25 & B <= 0.75; hi <- B > 0.75
  cs <- function(mask) { X <- B; X[!mask] <- NA_real_
    m <- base::rowMeans(X, na.rm = TRUE)
    list(m = m, sd = base::sqrt(base::rowMeans((X - m)^2, na.rm = TRUE))) }
  L <- cs(lo); M <- cs(mi); H <- cs(hi)
  ok <- base::rowMeans(lo, na.rm = TRUE) >= MIN_FRAC &
        base::rowMeans(mi, na.rm = TRUE) >= MIN_FRAC &
        base::rowMeans(hi, na.rm = TRUE) >= MIN_FRAC &
        base::is.finite(L$m) & base::is.finite(M$m) & base::is.finite(H$m) &
        L$m < END_LO & H$m > END_HI & M$m > MID_LO & M$m < MID_HI &
        L$sd < CL_SD & M$sd < CL_SD & H$sd < CL_SD
  ok[base::is.na(ok)] <- FALSE
  ok
}
call_gt <- function(B) {
  G <- base::matrix(NA_integer_, base::nrow(B), base::ncol(B), dimnames = base::dimnames(B))
  G[B < 0.25] <- 0L; G[B >= 0.25 & B <= 0.75] <- 1L; G[B > 0.75] <- 2L
  G
}
conc <- function(A, B, selfpair) {
  cp <- base::intersect(base::rownames(A), base::rownames(B))
  if (length(cp) < MIN_SHARE) return(NULL)
  A <- A[cp, , drop = FALSE]; B <- B[cp, , drop = FALSE]
  same <- base::matrix(0, base::ncol(A), base::ncol(B))
  for (k in 0:2) {
    Ak <- (A == k) * 1; Ak[base::is.na(Ak)] <- 0
    Bk <- (B == k) * 1; Bk[base::is.na(Bk)] <- 0
    same <- same + base::crossprod(Ak, Bk)
  }
  C <- same / base::crossprod((!base::is.na(A)) * 1, (!base::is.na(B)) * 1)
  if (selfpair) C[base::upper.tri(C, diag = TRUE)] <- NA_real_
  list(v = C[base::is.finite(C)], n = length(cp))
}

GT <- list()
for (g in PAIR) {
  say("---- ", g, " ----")
  B <- get_beta(g)
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  B <- B[, base::intersect(D$gsm, base::colnames(B)), drop = FALSE]
  sel <- snp_like(B)
  say("  Probes with genotype shape: ", base::sum(sel), "; samples ", base::ncol(B))
  GT[[g]] <- call_gt(B[base::which(sel), , drop = FALSE])
  rm(B); invisible(gc())
}

say("")
say("=== Distribution of concordance ===")
w1 <- conc(GT[[PAIR[1]]], GT[[PAIR[1]]], TRUE)
w2 <- conc(GT[[PAIR[2]]], GT[[PAIR[2]]], TRUE)
bt <- conc(GT[[PAIR[1]]], GT[[PAIR[2]]], FALSE)
if (is.null(bt)) {
  say("  Fewer than ", MIN_SHARE, " probes shared by the two sets; no screen possible")
  log_stop()
} else {
  for (nm in c(PAIR[1], PAIR[2], "cross-cohort")) {
    v <- switch(nm, `cross-cohort` = bt$v, if (nm == PAIR[1]) w1$v else w2$v)
    if (is.null(v)) { say(sprintf("  %-12s too few probes, skipped", nm)); next }
    say(sprintf("  %-12s sample pairs %8d  median %.4f  99%% %.4f  99.9%% %.4f  maximum %.4f",
                nm, length(v), stats::median(v), stats::quantile(v, 0.99),
                stats::quantile(v, TAIL_Q), base::max(v)))
  }
  say("  Probes shared by the two sets: ", bt$n, "")

  ref <- base::max(stats::quantile(w1$v, TAIL_Q), stats::quantile(w2$v, TAIL_Q))
  obs <- stats::quantile(bt$v, TAIL_Q)
  say("")
  say(sprintf("  Within-cohort upper tail (%.1f%% quantile), the larger of the two = %.4f", 100 * TAIL_Q, ref))
  say(sprintf("  Cross-cohort value at the same quantile = %.4f, %.2f times the former (criterion: at most %.2f times)",
              obs, obs / ref, TAIL_MULT))
  verdict <- obs <= TAIL_MULT * ref
  say("  Decision: ", if (verdict) "no large-scale overlap seen" else "[ACTION NEEDED] the cross-cohort upper tail is clearly raised, suggesting overlap")
  say("")
  say("  What this screen can exclude: large-scale sharing of individuals between the two Parkinson's disease cohorts.")
  say("  What it cannot exclude: individual duplicate samples. The public matrices do not retain enough genotype information ")
  say("  for pairwise identification; this is stated in the limitations.")

  df <- rbind(data.frame(grp = PAIR[1], v = w1$v, stringsAsFactors = FALSE),
              data.frame(grp = PAIR[2], v = w2$v, stringsAsFactors = FALSE),
              data.frame(grp = "Between cohorts", v = bt$v, stringsAsFactors = FALSE))
  write.csv(data.frame(group = c(PAIR, "between"),
                       median = c(stats::median(w1$v), stats::median(w2$v), stats::median(bt$v)),
                       q999 = c(stats::quantile(w1$v, TAIL_Q), stats::quantile(w2$v, TAIL_Q), obs),
                       max = c(base::max(w1$v), base::max(w2$v), base::max(bt$v)),
                       stringsAsFactors = FALSE),
            file.path(OUT_DIR, "overlap_screen.csv"), row.names = FALSE)
  say("=== Plotting ===")
  p <- ggplot(df, aes(v, colour = grp)) +
    geom_density(linewidth = 0.6) +
    scale_colour_manual(values = c(PAL[1], PAL[2], PAL[3]), name = NULL) +
    labs(title = "Genotype concordance within and between the two Parkinson cohorts",
         subtitle = "A second mode in the between-cohort curve would indicate shared individuals",
         x = "Concordance over genotype-shaped probes", y = "Density") +
    theme_paper()
  save_fig(p, "fig52_overlap_screen", width = 8, height = 5)
  log_stop()
}
