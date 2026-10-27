# ── LDA figures ───────────────────────────────────────────────────────────

lighten <- function(col, amount = 0.85) {
  rgbv <- grDevices::col2rgb(col) / 255
  mix  <- rgbv + (1 - rgbv) * amount
  grDevices::rgb(mix[1, ], mix[2, ], mix[3, ])
}

# Significance brackets, stacked shortest-span first so they never collide.
bracket_layers <- function(comps, groups, y0, step, tip, size = 6) {
  if (is.null(comps) || !nrow(comps)) return(list())
  comps$xi <- match(comps$A, groups)
  comps$xj <- match(comps$B, groups)
  comps$span <- abs(comps$xj - comps$xi)
  comps <- comps[order(comps$span, pmin(comps$xi, comps$xj)), ]
  comps$y <- y0 + step * (seq_len(nrow(comps)) - 1)

  seg <- do.call(rbind, lapply(seq_len(nrow(comps)), function(k) {
    r <- comps[k, ]; x1 <- min(r$xi, r$xj); x2 <- max(r$xi, r$xj)
    data.frame(x = c(x1, x1, x2), xend = c(x2, x1, x2),
               y = c(r$y, r$y, r$y), yend = c(r$y, r$y - tip, r$y - tip))
  }))
  list(
    ggplot2::geom_segment(data = seg,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.6, colour = "black", inherit.aes = FALSE),
    ggplot2::geom_text(data = comps,
      ggplot2::aes(x = (xi + xj) / 2, y = y + tip * 0.35, label = stars),
      size = size, fontface = "bold", vjust = 0, inherit.aes = FALSE)
  )
}

lda_bar_plot <- function(lda, df, cols, comps = NULL,
                         title = NULL, subtitle = NULL,
                         show_dots = TRUE, show_labels = TRUE,
                         ylab = "Organoid-forming frequency") {
  est <- lda$est
  g   <- lda$groups
  obs <- lda_observed(df[df$Group %in% g, , drop = FALSE], lda$mode)
  obs$Group <- factor(obs$Group, levels = g)

  top  <- max(c(est$CI_high, if (show_dots) obs$rate), na.rm = TRUE)
  n_br <- if (is.null(comps)) 0 else nrow(comps)
  step <- top * 0.13
  y0   <- top * 1.20
  ymax <- if (n_br) y0 + step * n_br else top * 1.15

  p <- ggplot2::ggplot(est, ggplot2::aes(Group, Frequency)) +
    ggplot2::geom_col(ggplot2::aes(fill = Group, colour = Group),
                      width = 0.6, linewidth = 1.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_low, ymax = CI_high),
                           width = 0.13, linewidth = 0.9, colour = "black")

  if (show_dots)
    p <- p + ggplot2::geom_point(
      data = obs, ggplot2::aes(Group, rate), inherit.aes = FALSE,
      shape = 21, size = 3.2, stroke = 1.1, colour = "grey25", fill = NA,
      position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1))

  if (show_labels)
    p <- p + ggplot2::geom_text(
      ggplot2::aes(y = CI_high + top * 0.045,
                   label = paste0("1/", round(1 / Frequency))),
      size = 5.4, fontface = "bold", vjust = 0)

  p <- p + bracket_layers(comps, g, y0, step, step * 0.28)

  p +
    ggplot2::scale_fill_manual(values = lighten(cols[g], 0.87), guide = "none") +
    ggplot2::scale_colour_manual(values = cols[g], guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, ymax),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = ylab, title = title, subtitle = subtitle) +
    theme_lab() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 15, face = "bold"))
}

# "cells" layout: organoids counted vs cells seeded, line through the origin.
lda_drc_plot <- function(lda, df, cols, title = NULL, subtitle = NULL,
                         x_step = NULL, ribbon = TRUE) {
  g   <- lda$groups
  est <- lda$est
  raw <- df[df$Group %in% g, , drop = FALSE]
  raw$Group <- factor(raw$Group, levels = g)

  x_max <- max(raw$CellCount, na.rm = TRUE) * 1.06
  xs    <- seq(0, x_max, length.out = 400)
  pred  <- do.call(rbind, lapply(seq_along(g), function(i) {
    data.frame(CellCount = xs, Group = g[i],
               fit   = est$Frequency[i] * xs,
               ci_lo = est$CI_low[i]    * xs,
               ci_hi = est$CI_high[i]   * xs)
  }))
  pred$Group <- factor(pred$Group, levels = g)
  y_max <- max(c(raw$Organoids, pred$ci_hi), na.rm = TRUE) * 1.05
  if (is.null(x_step) || !is.finite(x_step) || x_step <= 0)
    x_step <- diff(range(pretty(c(0, x_max), n = 5)))[1] / 5

  p <- ggplot2::ggplot()
  if (ribbon)
    p <- p + ggplot2::geom_ribbon(
      data = pred, ggplot2::aes(CellCount, ymin = ci_lo, ymax = ci_hi, fill = Group),
      alpha = 0.18, colour = NA) +
      ggplot2::geom_line(data = pred, ggplot2::aes(CellCount, ci_lo, colour = Group),
                         linewidth = 0.5, linetype = "dashed", alpha = 0.85) +
      ggplot2::geom_line(data = pred, ggplot2::aes(CellCount, ci_hi, colour = Group),
                         linewidth = 0.5, linetype = "dashed", alpha = 0.85)

  p +
    ggplot2::geom_line(data = pred, ggplot2::aes(CellCount, fit, colour = Group),
                       linewidth = 1.5) +
    ggplot2::geom_point(data = raw, ggplot2::aes(CellCount, Organoids, colour = Group),
                        shape = 16, size = 3.4) +
    ggplot2::scale_colour_manual(values = cols[g]) +
    ggplot2::scale_fill_manual(values = cols[g], guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, x_max),
                                breaks = seq(0, x_max, by = x_step),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0, y_max),
                                breaks = scales::pretty_breaks(6),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(x = "Cells seeded", y = "Organoids formed",
                  title = title, subtitle = subtitle) +
    theme_lab()
}

# "wells" layout: the classic ELDA plot — log fraction of negative wells.
lda_elda_plot <- function(lda, df, cols, title = NULL, subtitle = NULL) {
  g   <- lda$groups
  est <- lda$est
  raw <- df[df$Group %in% g, , drop = FALSE]
  raw$Group <- factor(raw$Group, levels = g)
  raw$neg   <- pmax((raw$Tested - raw$Positive) / raw$Tested, 1e-3)
  raw$logneg <- log(raw$neg)

  x_max <- max(raw$Dose, na.rm = TRUE) * 1.06
  xs    <- seq(0, x_max, length.out = 200)
  pred  <- do.call(rbind, lapply(seq_along(g), function(i)
    data.frame(Dose = xs, Group = g[i],
               fit = -est$Frequency[i] * xs,
               lo  = -est$CI_high[i]   * xs,
               hi  = -est$CI_low[i]    * xs)))
  pred$Group <- factor(pred$Group, levels = g)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = pred,
      ggplot2::aes(Dose, ymin = lo, ymax = hi, fill = Group), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(data = pred, ggplot2::aes(Dose, fit, colour = Group), linewidth = 1.4) +
    ggplot2::geom_point(data = raw, ggplot2::aes(Dose, logneg, colour = Group),
                        shape = 16, size = 3.2) +
    ggplot2::scale_colour_manual(values = cols[g]) +
    ggplot2::scale_fill_manual(values = cols[g], guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, x_max),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(x = "Cells per well", y = "log fraction of negative wells",
                  title = title, subtitle = subtitle) +
    theme_lab()
}
