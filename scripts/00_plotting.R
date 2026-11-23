# ==============================================================
# Shared blood DNA methylation signatures of comorbid disorders (PD x IBD)
# Shared plotting component
#
# Usage:
#   source(file.path(PROJ_ROOT, "scripts", "00_plotting.R"))
#   p <- ggplot(...) + theme_paper()
#   save_fig(p, "fig01_cohort_overview", width = 7, height = 4.5)
#
# Every call produces two files:
#   <name>.png  - 600 dpi, rendered with ragg (text and line quality clearly better than the default device)
#   <name>.pdf  - vector, via cairo_pdf; no loss on submission or enlargement
# and appends the file name, size, generating script and time to figures/figure_manifest.csv,
# so that it can later be checked that every figure cited in the text exists and is the latest version.
#
# All text in figures is in English. The manuscript goes to an English-language journal and the figures go straight into it;
# mixed-language figures cannot be used, and embedding CJK fonts in PDF often causes problems.
# ==============================================================

.PLOTENV <- new.env(parent = emptyenv())

# Colour-blind-friendly palette (Okabe-Ito)
PAL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
         "#E69F00", "#56B4E9", "#F0E442", "#999999")

.need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE),
                       logical(1))]
  if (length(miss)) {
    message("Installing plotting dependencies: ", paste(miss, collapse = ", "))
    install.packages(miss)
  }
  still <- pkgs[!vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE),
                        logical(1))]
  if (length(still)) stop("The following plotting dependencies could not be installed; figures cannot be produced: ",
                          paste(still, collapse = ", "))
}
.need(c("ggplot2", "ragg", "scales", "patchwork"))

suppressPackageStartupMessages({
  library(ggplot2)
})

theme_paper <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      panel.border = element_rect(linewidth = 0.4, colour = "grey40"),
      axis.text = element_text(colour = "grey20"),
      axis.title = element_text(face = "plain"),
      strip.background = element_rect(fill = "grey94", colour = "grey40",
                                      linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      legend.key = element_blank(),
      legend.background = element_blank(),
      plot.title = element_text(face = "bold", size = base_size + 1.5),
      plot.subtitle = element_text(colour = "grey30", size = base_size - 0.5),
      plot.caption = element_text(colour = "grey45", size = base_size - 2,
                                  hjust = 0)
    )
}

fig_dir <- function(proj_root = get0("PROJ_ROOT", ifnotfound = getwd())) {
  d <- file.path(proj_root, "output", "figures")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  if (!dir.exists(d)) stop("Failed to create figure directory: ", d)
  d
}

# p may be a ggplot object or a function without arguments (for base R graphics)
save_fig <- function(p, name, width = 7, height = 5, dpi = 600,
                     proj_root = get0("PROJ_ROOT", ifnotfound = getwd()),
                     script = get0("SCRIPT_NAME", ifnotfound = NA_character_)) {
  d <- fig_dir(proj_root)
  png_path <- file.path(d, paste0(name, ".png"))
  pdf_path <- file.path(d, paste0(name, ".pdf"))

  draw <- function() if (is.function(p)) p() else print(p)

  ragg::agg_png(png_path, width = width, height = height, units = "in",
                res = dpi, background = "white")
  draw(); grDevices::dev.off()

  grDevices::cairo_pdf(pdf_path, width = width, height = height,
                       onefile = FALSE)
  draw(); grDevices::dev.off()

  row <- data.frame(
    figure = name,
    png_kb = round(file.size(png_path) / 1024, 1),
    pdf_kb = round(file.size(pdf_path) / 1024, 1),
    width_in = width, height_in = height, dpi = dpi,
    script = script,
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE)
  man <- file.path(d, "figure_manifest.csv")
  if (file.exists(man)) {
    old <- read.csv(man, stringsAsFactors = FALSE)
    old <- old[old$figure != name, , drop = FALSE]   # For figures with the same name the latest one wins
    row <- rbind(old, row)
  }
  write.csv(row, man, row.names = FALSE)

  cat("  [figure] ", name, "  ", width, "x", height, " in @", dpi, " dpi",
      "  png ", row$png_kb[nrow(row)], " KB\n", sep = "")
  invisible(c(png = png_path, pdf = pdf_path))
}

# Large scatter plots are rasterised first and overlaid with vector axes, to avoid hundreds of thousands of points in the PDF
# (falls back to an ordinary point layer without error when ggrastr is not installed)
geom_point_big <- function(...) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterise(ggplot2::geom_point(...), dpi = 600)
  } else {
    ggplot2::geom_point(...)
  }
}
