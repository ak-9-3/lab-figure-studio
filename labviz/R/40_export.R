# ── Downloads ─────────────────────────────────────────────────────────────
# Two kinds of export:
#   1. the picture (PNG / PDF / SVG — vector formats stay editable)
#   2. a ZIP holding the user's data + a standalone R script that rebuilds
#      the exact figure they are looking at, so they can keep tweaking it
#      long after they close the browser tab.

save_plot <- function(file, plot, width, height, dpi = 450) {
  dev <- tools::file_ext(file)
  if (dev == "svg" && requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(file, plot, width = width, height = height, device = svglite::svglite)
  } else if (dev == "pdf") {
    ggplot2::ggsave(file, plot, width = width, height = height,
                    device = pdf_device())
  } else {
    ggplot2::ggsave(file, plot, width = width, height = height, dpi = dpi,
                    bg = "white", device = png_device())
  }
}

# Turn an R value into source code that recreates it.
deparse1v <- function(x) paste(deparse(x, width.cutoff = 500), collapse = "")

# Build a directory somewhere temporary, then zip whatever ended up in it.
zip_dir <- function(zipfile, build) {
  tmp <- file.path(tempdir(), paste0("export_", as.integer(runif(1, 1e6, 1e7))))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  build(tmp)
  old <- setwd(tmp); on.exit(setwd(old), add = TRUE, after = FALSE)
  files <- list.files(".", recursive = TRUE)
  # utils::zip() shells out to a zip binary, which does not exist in webR.
  # The zip package is pure C (miniz) and works everywhere.
  if (requireNamespace("zip", quietly = TRUE)) zip::zip(zipfile, files)
  else utils::zip(zipfile, files = files, flags = "-r9Xq")
}

copy_helpers <- function(tmp, helpers) {
  if (!length(helpers)) return(invisible())
  dir.create(file.path(tmp, "R"), showWarnings = FALSE)
  file.copy(helpers, file.path(tmp, "R", basename(helpers)), overwrite = TRUE)
}

# Bundle: data.csv + analysis.R + the shared helper files the script sources.
export_bundle <- function(zipfile, data, script, helpers = NULL, extra = list()) {
  zip_dir(zipfile, function(tmp) {
    utils::write.csv(data, file.path(tmp, "data.csv"), row.names = FALSE)
    writeLines(script, file.path(tmp, "analysis.R"))
    copy_helpers(tmp, helpers)
    for (nm in names(extra)) writeLines(extra[[nm]], file.path(tmp, nm))
  })
}

safe_name <- function(x) gsub("_+", "_", gsub("[^A-Za-z0-9._-]", "_", x))

# One figure per level of `batch_col`, plus per-level and pooled tables.
# This is the button that replaces "can you make me a figure" entirely: hand
# somebody a sheet with twenty experiments in it and they get twenty figures.
lda_batch_bundle <- function(zipfile, data, batch_col, mode, groups, cols, opts,
                             helpers = NULL, script = NULL, progress = NULL) {
  levels_ <- unique(as.character(data[[batch_col]]))
  zip_dir(zipfile, function(tmp) {
    for (d in c("Plots", "Tables", "Inputs")) dir.create(file.path(tmp, d))
    utils::write.csv(data, file.path(tmp, "Inputs", "data.csv"), row.names = FALSE)
    copy_helpers(tmp, helpers)
    if (!is.null(script)) writeLines(script, file.path(tmp, "analysis.R"))

    all_est <- list(); all_cmp <- list(); skipped <- character()

    for (i in seq_along(levels_)) {
      lv <- levels_[i]
      if (is.function(progress)) progress(i / length(levels_), lv)
      d <- data[as.character(data[[batch_col]]) == lv, , drop = FALSE]
      g <- intersect(groups, unique(as.character(d$Group)))
      fit <- tryCatch(lda_fit(d, g, mode = mode, conf = opts$conf),
                      error = function(e) NULL)
      if (is.null(fit) || length(g) < 1) { skipped <- c(skipped, lv); next }

      cmp <- tryCatch(
        lda_pairwise(fit, lda_pairs(g, if (opts$ref %in% g) opts$ref else g[1],
                                    opts$comp), adjust = opts$adjust),
        error = function(e) NULL)

      sub <- if (nzchar(opts$bar_sub)) paste0(opts$bar_sub, "  -  ", lv) else lv
      p_bar <- lda_bar_plot(fit, d, cols, cmp,
                            title = opts$bar_title, subtitle = sub,
                            show_dots = opts$dots, show_labels = opts$labels,
                            ylab = opts$ylab)
      p_crv <- if (identical(mode, "wells"))
        lda_elda_plot(fit, d, cols, opts$drc_title, lv)
      else
        lda_drc_plot(fit, d, cols, opts$drc_title, lv, x_step = opts$x_step)

      stem <- safe_name(lv)
      for (ext in c("pdf", "png"))
        save_plot(file.path(tmp, "Plots", sprintf("%s_frequency.%s", stem, ext)),
                  p_bar, opts$w, opts$h, 300)
      for (ext in c("pdf", "png"))
        save_plot(file.path(tmp, "Plots", sprintf("%s_curve.%s", stem, ext)),
                  p_crv, opts$w, opts$h, 300)

      est <- fit$est; est[[batch_col]] <- lv
      est$gof_p <- fit$gof_p
      utils::write.csv(est, file.path(tmp, "Tables",
                                      sprintf("%s_frequencies.csv", stem)),
                       row.names = FALSE)
      all_est[[lv]] <- est
      if (!is.null(cmp) && nrow(cmp)) {
        cmp[[batch_col]] <- lv
        utils::write.csv(cmp, file.path(tmp, "Tables",
                                        sprintf("%s_comparisons.csv", stem)),
                         row.names = FALSE)
        all_cmp[[lv]] <- cmp
      }
    }

    if (length(all_est))
      utils::write.csv(do.call(rbind, all_est),
                       file.path(tmp, "Tables", "ALL_frequencies.csv"),
                       row.names = FALSE)
    if (length(all_cmp))
      utils::write.csv(do.call(rbind, all_cmp),
                       file.path(tmp, "Tables", "ALL_comparisons.csv"),
                       row.names = FALSE)

    writeLines(c(
      "Limiting-dilution batch export",
      paste("Generated:", format(Sys.time())),
      paste("Split by:", batch_col),
      paste("Levels:", length(levels_)),
      if (length(skipped))
        paste("Skipped (model would not fit):", paste(skipped, collapse = ", ")),
      "",
      "Plots/   one frequency bar chart and one curve per level (PDF + PNG)",
      "Tables/  per-level frequencies and comparisons, plus pooled ALL_* files",
      "Inputs/  the data exactly as it was analysed",
      "R/       the fitting and plotting functions",
      "analysis.R  rerun or edit everything in RStudio",
      "",
      "PDF is vector: every element stays editable in Illustrator or Inkscape.",
      "gof_p in the frequency tables is the goodness-of-fit p-value; small",
      "values mean one frequency per group does not describe that level well."
    ), file.path(tmp, "README.txt"))
  })
}

lda_script <- function(mode, groups, cols, opts) {
  helper_files <- c("00_theme.R", "02_lda_stats.R", "03_lda_plots.R")
  paste0(
'# ═══════════════════════════════════════════════════════════════════════
#  Limiting-dilution analysis — reproducible script
#  Generated by the lab figure app on ', format(Sys.Date()), '
#
#  HOW TO USE
#    1. Keep this file next to data.csv and the R/ folder.
#    2. Open in RStudio and hit Source. Figures land in figures/.
#    3. Everything below is plain ggplot2 — change anything you like.
# ═══════════════════════════════════════════════════════════════════════

library(ggplot2); library(scales)
for (f in list.files("R", full.names = TRUE)) source(f)

df   <- read.csv("data.csv", stringsAsFactors = FALSE)
mode <- "', mode, '"

groups <- ', deparse1v(groups), '
cols   <- ', deparse1v(cols), '

# ── Fit ────────────────────────────────────────────────────────────────
lda <- lda_fit(df, groups, mode = mode, conf = ', opts$conf, ')
print(lda$est)

# ── Group comparisons (Wald on the log-frequency scale) ────────────────
pairs <- ', deparse1v(opts$pairs), '
comps <- lda_pairwise(lda, pairs, adjust = "', opts$adjust, '")
print(comps)

# ── Figures ────────────────────────────────────────────────────────────
dir.create("figures", showWarnings = FALSE)

p_bar <- lda_bar_plot(lda, df, cols, comps,
                      title    = ', deparse1v(opts$bar_title), ',
                      subtitle = ', deparse1v(opts$bar_sub), ',
                      show_dots = ', opts$dots, ', show_labels = ', opts$labels, ',
                      ylab = ', deparse1v(opts$ylab), ')
ggsave("figures/LDA_frequency.pdf", p_bar,
       width = ', opts$w, ', height = ', opts$h, ', device = pdf_device())

p_curve <- ', if (mode == "cells")
  paste0('lda_drc_plot(lda, df, cols,
                     title    = ', deparse1v(opts$drc_title), ',
                     subtitle = ', deparse1v(opts$drc_sub), ',
                     x_step   = ', opts$x_step, ')')
else
  paste0('lda_elda_plot(lda, df, cols,
                      title    = ', deparse1v(opts$drc_title), ',
                      subtitle = ', deparse1v(opts$drc_sub), ')'), '
ggsave("figures/LDA_curve.pdf", p_curve,
       width = ', opts$w, ', height = ', opts$h, ', device = pdf_device())

write.csv(lda$est, "figures/frequencies.csv", row.names = FALSE)
if (!is.null(comps)) write.csv(comps, "figures/comparisons.csv", row.names = FALSE)
')
}
