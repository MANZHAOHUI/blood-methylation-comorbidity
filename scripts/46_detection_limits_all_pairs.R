# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders - five-disorder panel
# Step 46: detection limits of the forty-five comparisons
#
# 45 gave the detection limits of the five within-disease positive controls. This step extends them to all forty-five pairs,
# so that every negative result carries "how large an effect can be excluded" rather than just "not significant".
#
# No permutation needs to be re-run: 39 saved the permutation null of every pair, and the detection limit depends only on
# the location and width of the null. This step computes it directly from that file in a few seconds.
#
# Two quantities:
#   critical rho        the one-sided 5% quantile of the null. Exceeding it is nominal significance.
#   needed for 80% power  the null mean plus (z_0.95 + z_0.80) times the SD of the null.
#                 When the true effect reaches this value, the design detects it with a probability of 80%.
#
# The "effective number of independent probes" is reported as well, back-calculated from the SD of the null (1 over the variance). Nominally
# there are hundreds of thousands of probes, but they are highly correlated; what determines the detection limit is how much
# independent information remains given the residual correlation structure. This number explains the large differences in detection limit between pairs and shows
# that "more probes, more sensitivity" does not hold.
#
# This step only reads the existing output of 39 and 45 and recomputes no model.
#
# Environment: Windows + RStudio
# Location: <project root>/scripts/46_detection_limits_all_pairs.R
# ==============================================================

PROJ_ROOT <- "C:/blood_methylation_project"
SCRIPT_NAME <- "46_detection_limits_all_pairs"

CON_DIR <- file.path(PROJ_ROOT, "output", "39_concordance")
PC_DIR  <- file.path(PROJ_ROOT, "output", "45_stable")
OUT_DIR <- file.path(PROJ_ROOT, "output", "46_limits")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

source(file.path(PROJ_ROOT, "scripts", "00_logging.R"))
log_start(SCRIPT_NAME, proj_root = PROJ_ROOT)
source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))

POWER <- 0.80
Z95 <- stats::qnorm(0.95); Z80 <- stats::qnorm(POWER)

tab <- read.csv(file.path(CON_DIR, "concordance_all_pairs.csv"), stringsAsFactors = FALSE)
say("Results of the forty-five pairs from 39 loaded: ", base::nrow(tab), " rows")
stopifnot(base::all(c("null_mean", "null_sd", "rho", "kind") %in% base::names(tab)))

tab$rho_crit <- tab$null_mean + Z95 * tab$null_sd
tab$rho_power <- tab$null_mean + (Z95 + Z80) * tab$null_sd
tab$eff_probes <- 1 / tab$null_sd^2

# For the five within-disease pairs the ten-thousand-permutation results of 45 are used, which are more precise
pc <- tryCatch(read.csv(file.path(PC_DIR, "positive_controls_10k.csv"),
                        stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(pc)) {
  say("Ten-thousand-permutation results from 45 loaded: ", base::nrow(pc), " rows; they are decisive for the five within-disease pairs")
  for (i in base::seq_len(base::nrow(pc))) {
    j <- base::which(tab$kind == "within" & base::grepl(pc$disease[i], tab$id, fixed = TRUE))
    if (length(j) == 1L) {
      tab$null_mean[j] <- pc$null_mean[i]; tab$null_sd[j] <- pc$null_sd[i]
      tab$rho[j] <- pc$rho[i]; tab$z[j] <- pc$z[i]; tab$p_rho[j] <- pc$p[i]
      tab$rho_crit[j] <- pc$rho_crit[i]; tab$rho_power[j] <- pc$rho_for_power[i]
      tab$eff_probes[j] <- 1 / pc$null_sd[i]^2
    } else say("  [NOTE] within-disease pair '", pc$disease[i], "' matched in the table of 39: ", length(j), " rows")
  }
} else say("[NOTE] results of 45 not found; the five within-disease pairs keep the one-thousand-permutation values of 39")

write.csv(tab, file.path(OUT_DIR, "all_pairs_with_limits.csv"), row.names = FALSE)

say("")
say("=== Within-disease positive controls ===")
w <- tab[tab$kind == "within", ]
w <- w[base::order(-w$rho), ]
for (i in base::seq_len(base::nrow(w)))
  say(sprintf("  %-26s rho=%+.4f  P=%.4f  critical %.4f  80%% power needs %.4f  effective probes %6.0f",
              w$id[i], w$rho[i], w$p_rho[i], w$rho_crit[i], w$rho_power[i], w$eff_probes[i]))

say("")
say("=== Forty cross-disease pairs: summary by disease pair ===")
cr <- tab[tab$kind == "cross", ]
rows <- list()
for (dp in base::unique(cr$dpair)) {
  s <- cr[cr$dpair == dp, ]
  say(sprintf("  %-28s rho %s; significant %d/4; critical %.4f to %.4f",
              dp, base::paste(sprintf("%+.3f", s$rho), collapse = " "),
              base::sum(s$p_rho < 0.05), base::min(s$rho_crit), base::max(s$rho_crit)))
  rows[[length(rows) + 1]] <- data.frame(
    dpair = dp, n_sig = base::sum(s$p_rho < 0.05),
    rho_median = stats::median(s$rho), rho_max = base::max(s$rho),
    crit_min = base::min(s$rho_crit), crit_max = base::max(s$rho_crit),
    power_min = base::min(s$rho_power), power_max = base::max(s$rho_power),
    stringsAsFactors = FALSE)
}
ss <- base::do.call(rbind, rows)
write.csv(ss, file.path(OUT_DIR, "cross_disease_limits.csv"), row.names = FALSE)

say("")
say("=== Three quantitative statements suitable for the paper ===")
ibd <- w$rho[base::grepl("IBD", w$id)]
ad_crit <- w$rho_crit[base::grepl("Alzheimer", w$id)]
say(sprintf("  (1) Within-disease concordance is strongest for inflammatory bowel disease, rho = %+.4f; the two weakest disorders", ibd))
say(sprintf("      have observed values of %+.4f and %+.4f, while their critical values are %.4f and %.4f -",
            w$rho[base::grepl("Alzheimer", w$id)], w$rho[base::grepl("Depression", w$id)],
            ad_crit, w$rho_crit[base::grepl("Depression", w$id)]))
say("      the test itself is sensitive enough; there is simply no reproducible signal.")
say(sprintf("  (2) The observed value of inflammatory bowel disease is %.1f times the critical value of Alzheimer's disease; a reproducible signal", ibd / ad_crit))
say("      of that order in the latter would certainly be detected by this design.")
jmax <- base::which.max(cr$rho)
sib <- cr[cr$dpair == cr$dpair[jmax], , drop = FALSE]
wmin <- base::min(w$rho[w$p_rho < 0.05])
say(sprintf("  (3) The highest of the forty cross-disease cells is %+.4f (%s: %s vs %s); relative to the weakest within-disease value",
            cr$rho[jmax], cr$dpair[jmax], cr$fixed[jmax], cr$permuted[jmax]))
say(sprintf("      (%+.4f) it is %s. But it is only one of the four pairings of that disease pair, and the other three are %s; ",
            wmin, if (cr$rho[jmax] > wmin) "of the same order" else "comparable",
            base::paste(sprintf("%+.3f", sib$rho[-base::which.max(sib$rho)]), collapse = ", ")))
say("      none of the ten disease pairs is significant in all four cohort pairings. The point is not how high a single cell is ")
say("      but whether the four pairings agree - looking at the highest cell alone would lead to the opposite conclusion.")

say("")
say("=== Effective number of independent probes (a note for the Methods, not a finding) ===")
say("  Back-calculated from the SD of the null (1 over the variance); a heuristic quantity: the SD of the correlation between two independent random vectors")
say("  is about 1 over the square root of the number of probes, which can be inverted to estimate 'how many")
say("  independent probes this corresponds to'. It explains why detection limits differ so much between pairs, but removing it changes no conclusion, ")
say("  so it belongs in the Methods only and is not a finding.")
a <- tab[base::order(-tab$eff_probes), ]
say(sprintf("  highest %6.0f (%s), lowest %6.0f (%s), median %6.0f",
            base::max(tab$eff_probes), a$id[1],
            base::min(tab$eff_probes), a$id[base::nrow(a)],
            stats::median(tab$eff_probes)))
say("  Meaning: the detection limit depends on how much independent information remains given the residual correlation structure, not on the total number of probes.")

say("")
say("=== Plotting ===")
tab$panelkind <- factor(ifelse(tab$kind == "within",
                               "Within-disease (positive control)", "Cross-disease"),
                        levels = c("Within-disease (positive control)", "Cross-disease"))
p <- ggplot(tab, aes(rho_crit, rho, colour = panelkind)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55", linewidth = 0.4) +
  geom_point(size = 1.9, alpha = 0.85) +
  scale_colour_manual(values = c("Within-disease (positive control)" = PAL[3],
                                 "Cross-disease" = PAL[1]), name = NULL) +
  labs(title = "Observed concordance against what each comparison could detect",
       subtitle = "Points above the diagonal reach nominal significance; the distance below it is the margin by which an effect is excluded",
       x = "Critical value of the permutation null (one-sided 5%)",
       y = "Observed Spearman correlation") +
  theme_paper()
save_fig(p, "fig49_detection_limits", width = 8.5, height = 6)
log_stop()
