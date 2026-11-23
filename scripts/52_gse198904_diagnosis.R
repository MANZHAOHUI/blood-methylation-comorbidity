# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 52: why GSE198904 does not pass the sex positive control
#
# In the results of 51 the median sex concordance among nine cohorts is 0.12 to 0.20, whereas the median between GSE198904
# and the other nine is only +0.0007; all negative values among the forty-five pairs involve it. Without it,
# the median of the remaining thirty-six pairs is +0.1845 and the minimum +0.0493, all positive without exception.
#
# A real signal reproduced by all cohorts disappears or even reverses in this one set - this is what a positive control
# is for: it identified a problematic cohort. This cohort is one of the two for depression,
# the within-disease positive control of depression fails (-0.0089), and the two depression cohorts do not even agree on sex
# (+0.0418). The negative conclusion for depression may therefore be caused by a defect of this cohort rather than be
# a property of depression. This must be clarified.
#
# Two explanations to be distinguished:
#   (A) model: the cohort uses 50 methylation principal components. Only components associated with case-control status were excluded;
#       those associated with sex were kept, so the sex effect is absorbed by the components.
#   (B) data: sex is wrongly annotated, or samples are misaligned.
#   Misaligned samples under (B) are unlikely, however - DNAm age and recorded age correlate at 0.926 in this cohort,
#   which would be impossible if columns and samples were shifted.
#
# Test: estimate the sex effect of this cohort at K = 0, 5, 10, 20, 50 and compare it with the other
# nine sets. If concordance clearly recovers as K decreases, it is (A); if it stays near zero, it is the sex annotation
# problem of (B). The correlation of each of the fifty components with sex is reported as well.
#
# The disease effect is checked in the same way: does the within-disease concordance of depression change with K? Concordance is a rank
# correlation and unaffected by the inflation factor, so inflation at small K does not interfere with this comparison.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/52_gse198904_diagnosis.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
DATA_ROOT <- "D:/blood_methylation_data"
SCRIPT_NAME <- "52_gse198904_diagnosis"

SET_DIR <- file.path(DATA_ROOT, "analysis_sets")
RES_DIR <- file.path(DATA_ROOT, "ewas_final")
OUT_DIR <- file.path(PROJ_ROOT, "output", "52_gse198904")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

TARGET <- "GSE198904"
REFS <- c("GSE145361", "GSE111629", "GSE84727", "GSE80417", "GSE338167",
          "GSE144858", "GSE273970", "GSE87648", "GSE112611")
K_GRID <- c(0L, 5L, 10L, 20L, 50L)
BETA_CLIP <- 0.001

load_M <- function(g) {
  D <- base::readRDS(file.path(RES_DIR, base::paste0(g, "_design.rds")))
  obj <- base::readRDS(file.path(SET_DIR, base::paste0(g, "_analysis_filtered.rds")))
  b <- obj$beta[D$probes, D$gsm, drop = FALSE]
  rm(obj); invisible(gc())
  b[!base::is.finite(b)] <- NA
  b[b < BETA_CLIP] <- BETA_CLIP; b[b > 1 - BETA_CLIP] <- 1 - BETA_CLIP
  list(M = base::log2(b / (1 - b)), D = D)
}
fwl <- function(Mt, X0, ex) {
  q0 <- qr(X0)
  gr <- qr.resid(q0, ex)
  base::as.numeric(base::crossprod(gr, qr.resid(q0, Mt))) / base::sum(gr^2)
}

say("=== 1. Sex effect of the other nine sets (reference) ===")
REF <- list()
for (g in REFS) {
  L <- load_M(g); Dd <- L$D$design
  sc <- base::grep("^sex", base::colnames(Dd), value = TRUE)
  X0 <- Dd[, base::setdiff(base::colnames(Dd), sc), drop = FALSE]
  REF[[g]] <- list(eff = fwl(base::t(L$M), X0, base::as.numeric(Dd[, sc])),
                   probes = L$D$probes)
  say("  ", g, " finished")
  rm(L); invisible(gc())
}

say("")
say("=== 2. Design and principal components of ", TARGET, " ===")
L <- load_M(TARGET); Dd <- L$D$design; Mt <- base::t(L$M)
rm(L); invisible(gc())
sc <- base::grep("^sex", base::colnames(Dd), value = TRUE)
mp <- base::grep("^mpc", base::colnames(Dd), value = TRUE)
sexv <- base::as.numeric(Dd[, sc])
say("  Sex term '", sc, "': number coded 1: ", base::sum(sexv == 1), " / ", length(sexv),
    " (", sprintf("%.1f%%", 100 * base::mean(sexv == 1)), ")")
say("  Methylation principal components in the design: ", length(mp), "")
rs <- base::sapply(mp, function(k) stats::cor(Dd[, k], sexv))
o <- base::order(-base::abs(rs))
say("  Five principal components most strongly associated with sex: ")
for (k in utils::head(o, 5))
  say(sprintf("    %-6s r = %+.3f", mp[k], rs[k]))
say(sprintf("  principal components with |r| > 0.3: %d; with |r| > 0.5: %d",
            base::sum(base::abs(rs) > 0.3), base::sum(base::abs(rs) > 0.5)))
base_cov <- base::setdiff(base::colnames(Dd), c(sc, mp))

say("")
say("=== 3. Sex effect as a function of K ===")
res <- list()
for (K in K_GRID) {
  keep <- c(base_cov, utils::head(mp, K))
  X0 <- Dd[, keep, drop = FALSE]
  q <- qr(X0)
  if (q$rank < base::ncol(X0)) X0 <- X0[, base::sort(q$pivot[base::seq_len(q$rank)]), drop = FALSE]
  e <- fwl(Mt, X0, sexv)
  rr <- base::vapply(REFS, function(g) {
    cm <- base::intersect(REF[[g]]$probes, base::colnames(Mt))
    stats::cor(REF[[g]]$eff[base::match(cm, REF[[g]]$probes)],
               e[base::match(cm, base::colnames(Mt))], method = "spearman")
  }, numeric(1))
  say(sprintf("  K=%2d  concordance with the sex effect of the nine sets: median %+.4f  minimum %+.4f  maximum %+.4f",
              K, stats::median(rr), base::min(rr), base::max(rr)))
  res[[length(res) + 1]] <- data.frame(K = K, median = stats::median(rr),
    min = base::min(rr), max = base::max(rr), stringsAsFactors = FALSE)
}
rt <- base::do.call(rbind, res)
write.csv(rt, file.path(OUT_DIR, "sex_concordance_by_K.csv"), row.names = FALSE)

say("")
say("=== 4. Disease effect as a function of K (within-disease concordance of depression) ===")
EA <- base::readRDS(file.path(RES_DIR, "GSE273970_ewas_final.rds"))
gv <- base::as.numeric(Dd[, "groupcase"])
dres <- list()
for (K in K_GRID) {
  keep <- base::setdiff(c(base_cov, sc, utils::head(mp, K)), "groupcase")
  X0 <- Dd[, keep, drop = FALSE]
  q <- qr(X0)
  if (q$rank < base::ncol(X0)) X0 <- X0[, base::sort(q$pivot[base::seq_len(q$rank)]), drop = FALSE]
  e <- fwl(Mt, X0, gv)
  cm <- base::intersect(EA$probe, base::colnames(Mt))
  r <- stats::cor(EA$effect[base::match(cm, EA$probe)], e[base::match(cm, base::colnames(Mt))],
                  method = "spearman")
  say(sprintf("  K=%2d  concordance of the disease effect with GSE273970 %+.4f", K, r))
  dres[[length(dres) + 1]] <- data.frame(K = K, rho = r, stringsAsFactors = FALSE)
}
write.csv(base::do.call(rbind, dres), file.path(OUT_DIR, "disease_concordance_by_K.csv"),
          row.names = FALSE)

say("")
say("=== Interpretation ===")
rise <- rt$median[rt$K == 0] - rt$median[rt$K == 50]
say(sprintf("  From K = 50 to K = 0 the median sex concordance changes by %+.4f", rise))
if (rt$median[rt$K == 0] >= 0.10) {
  say("  At K=0 the sex concordance returns to the level of the other cohorts: the principal components absorbed the sex effect (explanation A).")
  say("  The data are sound; the problem lies in the model.")
} else {
  say("  Even without principal components the sex concordance stays far below the other cohorts: the components are not the cause.")
  say("  The sex annotation of this cohort itself is suspect (explanation B) and needs to be checked against the original publication.")
}
log_stop()
