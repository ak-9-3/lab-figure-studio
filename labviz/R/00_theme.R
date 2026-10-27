# ── Shared look & feel ────────────────────────────────────────────────────
# One theme, one palette, used by every module so all figures from the app
# belong to the same visual family.

theme_lab <- function(base_size = 16) {
  ggplot2::theme_classic(base_size = base_size) %+replace%
    ggplot2::theme(
      text         = ggplot2::element_text(family = "sans", colour = "black"),
      axis.line    = ggplot2::element_line(colour = "black", linewidth = 0.8),
      axis.ticks   = ggplot2::element_line(colour = "black", linewidth = 0.55),
      axis.title.x = ggplot2::element_text(size = base_size + 1, face = "bold",
                                           margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(size = base_size + 1, face = "bold",
                                           angle = 90, vjust = 0.5,
                                           margin = ggplot2::margin(r = 10)),
      axis.text    = ggplot2::element_text(size = base_size - 2, colour = "black"),
      plot.title   = ggplot2::element_text(size = base_size + 3, face = "bold",
                                           hjust = 0.5,
                                           margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = base_size - 4, colour = "grey35",
                                            hjust = 0.5,
                                            margin = ggplot2::margin(b = 12)),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = base_size - 3, face = "bold"),
      legend.position = "top",
      plot.margin  = ggplot2::margin(10, 16, 10, 12)
    )
}

# Colour-blind-safe defaults; grey first so controls/NTC land on grey.
PALETTES <- list(
  `Lab default` = c("#9B9B9B", "#D55E00", "#0072B2", "#009E73",
                    "#CC79A7", "#E69F00", "#56B4E9", "#3D405B"),
  Okabe_Ito     = c("#999999", "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
  Muted         = c("#AAAAAA", "#88CCEE", "#CC6677", "#DDCC77",
                    "#117733", "#332288", "#AA4499", "#44AA99"),
  Viridis       = c("#440154", "#46327E", "#365C8D", "#277F8E",
                    "#1FA187", "#4AC16D", "#A0DA39", "#FDE725")
)

pal_for <- function(groups, name = "Lab default") {
  p <- PALETTES[[name]] %||% PALETTES[[1]]
  stats::setNames(p[((seq_along(groups) - 1) %% length(p)) + 1], groups)
}

# p-value -> stars, the convention used in the lab's figures.
p_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "ns"))))
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 4)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# TRUE when running inside webR/shinylive, where there is no shell, no system
# zip binary, and no cairo. Several helpers branch on this.
is_wasm <- function() {
  identical(R.version$os, "emscripten") ||
    grepl("emscripten", R.version$platform, fixed = TRUE)
}

# cairo_pdf gives nicer font handling, but it is absent in webR and on some
# macOS R builds it reports itself as available and then fails at draw time.
# Probe it once for real, cache the answer, fall back to the base PDF device.
pdf_device <- local({
  cached <- NULL
  function() {
    if (is.null(cached)) {
      ok <- FALSE
      if (!is_wasm()) {
        f  <- tempfile(fileext = ".pdf")
        ok <- tryCatch({
          grDevices::cairo_pdf(f); graphics::plot.new(); grDevices::dev.off()
          file.exists(f) && file.size(f) > 500
        }, warning = function(w) FALSE, error = function(e) FALSE)
        while (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE)
        unlink(f)
      }
      cached <<- if (isTRUE(ok)) grDevices::cairo_pdf else grDevices::pdf
    }
    cached
  }
})

# PNG downloads use R's own png(). ragg would need textshaping, and the web
# build downloads every referenced package up front — a 13 MB Unicode library
# on every first visit for no visible difference.
png_device <- function() grDevices::png

# A colour picker using the browser's own <input type="color">. The
# colourpicker package did the same job but pulled in R Markdown, knitr and
# shinyjs — several megabytes of download for one widget.
colorInput <- function(inputId, label, value = "#000000") {
  div(class = "shiny-input-container mb-2",
      tags$label(class = "control-label small", `for` = inputId, label),
      tags$input(id = inputId, type = "color", value = value,
                 class = "lab-color form-control form-control-color w-100"))
}

COLOR_INPUT_JS <- "
$(function () {
  var b = new Shiny.InputBinding();
  $.extend(b, {
    find: function (scope) { return $(scope).find('input[type=color].lab-color'); },
    getValue: function (el) { return el.value; },
    setValue: function (el, v) { el.value = v; },
    subscribe: function (el, cb) { $(el).on('input.lab change.lab', function () { cb(true); }); },
    unsubscribe: function (el) { $(el).off('.lab'); },
    receiveMessage: function (el, d) { if (d.value) { el.value = d.value; $(el).trigger('change'); } }
  });
  Shiny.inputBindings.register(b, 'lab.color');
});"

# Plain HTML tables: numbers are already rounded, so print them as written.
as_display <- function(d) {
  d[] <- lapply(d, function(x) if (is.numeric(x)) as.character(x) else x)
  d
}
