# ── LDA figures ───────────────────────────────────────────────────────────

lighten <- function(col, amount = 0.85) {
  rgbv <- grDevices::col2rgb(col) / 255
  mix  <- rgbv + (1 - rgbv) * amount
  stats::setNames(grDevices::rgb(mix[1, ], mix[2, ], mix[3, ]), names(col))
}

# ── Panels ────────────────────────────────────────────────────────────────
# Lab group names usually carry two things: "Ad16-4oxoRA" is cell line Ad16
# treated with 4oxoRA. Printing twelve of those side by side under one axis is
# unreadable, so when every name splits the same way — and the second half
# repeats across the first — draw one panel per line with the treatment under
# each bar. Returns NULL when the names do not have that shape.
group_panels <- function(groups) {
  g <- as.character(groups)
  if (length(g) < 4) return(NULL)
  for (sep in c("-", "_", " | ", ": ", " ")) {
    pos <- regexpr(sep, g, fixed = TRUE)
    if (any(pos < 1)) next
    pre <- trimws(substr(g, 1, pos - 1))
    lab <- trimws(substr(g, pos + nchar(sep), nchar(g)))
    if (!all(nzchar(pre) & nzchar(lab))) next
    n_pan <- length(unique(pre))
    if (n_pan < 2 || n_pan == length(g) || !any(duplicated(lab))) next
    return(data.frame(Group = g,
                      Panel = factor(pre, levels = unique(pre)),
                      Label = factor(lab, levels = unique(lab)),
                      stringsAsFactors = FALSE))
  }
  NULL
}

use_panels <- function(groups, panels = "auto") {
  if (isFALSE(panels)) return(NULL)
  group_panels(groups)
}

# Default colours: when groups form panels, the same treatment gets the same
# colour in every panel (DMSO grey everywhere, ATRA orange everywhere).
lda_default_cols <- function(groups, palette = "Lab default") {
  p <- group_panels(groups)
  if (is.null(p)) return(pal_for(groups, palette))
  lab <- levels(p$Label)
  pl  <- pal_for(lab, palette)
  stats::setNames(unname(pl[as.character(p$Label)]), groups)
}

# Text shrinks as bars get narrower, so labels never run into each other.
label_size <- function(k) if (k <= 6) 5.4 else if (k <= 10) 4.4 else if (k <= 16) 3.6 else 3

# ── Significance marks ────────────────────────────────────────────────────
# A handful of comparisons → brackets, stacked shortest-span first.
# More than that → stars printed above each bar (vs. the reference), because
# fifteen stacked brackets is the jumble this is trying to avoid.
MAX_BRACKETS <- 4

plan_marks <- function(comps, layout) {
  if (is.null(comps) || !nrow(comps)) return(NULL)
  comps$pi <- layout$Panel[match(comps$A, layout$Group)]
  comps$pj <- layout$Panel[match(comps$B, layout$Group)]
  comps <- comps[!is.na(comps$pi) & !is.na(comps$pj) &
                 as.character(comps$pi) == as.character(comps$pj), , drop = FALSE]
  if (!nrow(comps)) return(NULL)
  comps$Panel <- comps$pi
  comps$xi <- layout$x[match(comps$A, layout$Group)]
  comps$xj <- layout$x[match(comps$B, layout$Group)]

  per_panel <- max(table(as.character(comps$Panel)))
  if (per_panel <= MAX_BRACKETS) {
    comps$kind <- "bracket"
  } else {
    # Keep only comparisons against each panel's most-compared group.
    keep <- unlist(lapply(split(seq_len(nrow(comps)), as.character(comps$Panel)), function(i) {
      hub <- names(which.max(table(c(comps$A[i], comps$B[i]))))
      i[comps$A[i] == hub | comps$B[i] == hub]
    }))
    comps <- comps[sort(keep), , drop = FALSE]
    hub <- vapply(split(seq_len(nrow(comps)), as.character(comps$Panel)), function(i)
      names(which.max(table(c(comps$A[i], comps$B[i])))), character(1))
    comps$hub   <- hub[as.character(comps$Panel)]
    comps$other <- ifelse(comps$A == comps$hub, comps$B, comps$A)
    comps$x     <- layout$x[match(comps$other, layout$Group)]
    comps$kind  <- "star"
  }
  comps
}

bracket_layers <- function(marks, y0, step, tip, size = 6) {
  if (is.null(marks) || !nrow(marks) || marks$kind[1] != "bracket") return(list())
  marks$span <- abs(marks$xj - marks$xi)
  marks <- do.call(rbind, lapply(split(marks, as.character(marks$Panel)), function(m) {
    m <- m[order(m$span, pmin(m$xi, m$xj)), , drop = FALSE]
    m$y <- y0 + step * (seq_len(nrow(m)) - 1)
    m
  }))
  seg <- do.call(rbind, lapply(seq_len(nrow(marks)), function(k) {
    r <- marks[k, ]; x1 <- min(r$xi, r$xj); x2 <- max(r$xi, r$xj)
    data.frame(Panel = r$Panel, x = c(x1, x1, x2), xend = c(x2, x1, x2),
               y = r$y, yend = c(r$y, r$y - tip, r$y - tip))
  }))
  list(
    ggplot2::geom_segment(data = seg,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.6, colour = "black", inherit.aes = FALSE),
    ggplot2::geom_text(data = marks,
      ggplot2::aes(x = (xi + xj) / 2, y = y + tip * 0.35, label = stars),
      size = size, fontface = "bold", vjust = 0, inherit.aes = FALSE)
  )
}

star_layers <- function(marks, y, size = 5) {
  if (is.null(marks) || !nrow(marks) || marks$kind[1] != "star") return(list())
  marks$y <- y
  list(ggplot2::geom_text(data = marks,
    ggplot2::aes(x = x, y = y, label = stars),
    size = size, fontface = "bold", vjust = 0, inherit.aes = FALSE))
}

# ── Bar chart ─────────────────────────────────────────────────────────────
lda_bar_plot <- function(lda, df, cols, comps = NULL,
                         title = NULL, subtitle = NULL,
                         show_dots = TRUE, show_labels = TRUE,
                         ylab = "Organoid-forming frequency", panels = "auto") {
  est <- lda$est
  g   <- lda$groups
  obs <- lda_observed(df[df$Group %in% g, , drop = FALSE], lda$mode)

  pan <- use_panels(g, panels)
  layout <- if (is.null(pan))
    data.frame(Group = g, Panel = factor("all"), Label = factor(g, levels = g),
               stringsAsFactors = FALSE)
  else pan
  # x position inside each panel, in the order the labels first appear
  layout$x <- ave(as.integer(layout$Label), layout$Panel,
                  FUN = function(v) match(v, sort(unique(v))))
  k <- max(table(layout$Panel))                # bars in the widest panel

  est <- merge(est, layout, by = "Group", sort = FALSE)
  obs <- merge(obs, layout, by = "Group", sort = FALSE)
  est$Group <- factor(est$Group, levels = g)

  dot_top <- if (show_dots) tapply(obs$rate, obs$Group, max) else NULL
  est$label_y <- pmax(est$CI_high,
                      if (show_dots) dot_top[as.character(est$Group)] else -Inf,
                      na.rm = TRUE)

  top   <- max(c(est$CI_high, if (show_dots) obs$rate), na.rm = TRUE)
  marks <- plan_marks(comps, layout)
  n_br  <- if (is.null(marks) || marks$kind[1] != "bracket") 0
           else max(table(as.character(marks$Panel)))
  step  <- top * 0.13
  y0    <- top * 1.20
  ymax  <- if (n_br) y0 + step * n_br
           else if (!is.null(marks)) top * 1.30 else top * 1.15

  p <- ggplot2::ggplot(est, ggplot2::aes(Label, Frequency)) +
    ggplot2::geom_col(ggplot2::aes(fill = Group, colour = Group),
                      width = 0.6, linewidth = if (k > 10) 0.9 else 1.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_low, ymax = CI_high),
                           width = 0.13, linewidth = 0.9, colour = "black")

  if (show_dots)
    p <- p + ggplot2::geom_point(
      data = obs, ggplot2::aes(Label, rate), inherit.aes = FALSE,
      shape = 21, size = if (k > 10) 2.2 else 3.2, stroke = 1, colour = "grey25", fill = NA,
      position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1))

  if (show_labels)
    p <- p + ggplot2::geom_text(
      ggplot2::aes(y = label_y + top * 0.04,
                   label = paste0("1/", round(1 / Frequency))),
      size = label_size(k), fontface = "bold", vjust = 0)

  p <- p + bracket_layers(marks, y0, step, step * 0.28)
  p <- p + star_layers(marks, top * 1.16, size = label_size(k))

  star_note <- if (!is.null(marks) && marks$kind[1] == "star")
    sprintf("Stars: vs %s. Every comparison is in the table.",
            paste(unique(if (is.null(pan)) marks$hub
                         else as.character(pan$Label[match(marks$hub, pan$Group)])),
                  collapse = " / "))

  p <- p +
    ggplot2::scale_fill_manual(values = lighten(cols[g], 0.87), guide = "none") +
    ggplot2::scale_colour_manual(values = cols[g], guide = "none") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(limits = c(0, ymax),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = ylab, title = title, subtitle = subtitle,
                  caption = star_note) +
    theme_lab() +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 10, colour = "grey35"))

  # One panel per line, each showing only its own treatments. Brackets and
  # stars use the bar's position inside its panel (layout$x), which is how a
  # free discrete scale numbers them.
  if (!is.null(pan))
    p <- p + ggplot2::facet_grid(~ Panel, scales = "free_x", space = "free_x")

  # Tilt the labels once, all together, they would run into each other.
  chars  <- sum(nchar(as.character(layout$Label))) + 2 * nrow(layout)
  angled <- k > 6 || chars > 32
  p + ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = if (k > 12) 11 else if (k > 6) 13 else 15, face = "bold",
      angle = if (angled) 45 else 0, hjust = if (angled) 1 else 0.5,
      vjust = if (angled) 1 else 0.5),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 16, face = "bold"),
    panel.spacing.x = ggplot2::unit(1.2, "lines"))
}

# "cells" layout: organoids counted vs cells seeded, line through the origin.
# Groups that form panels get one panel each; more than six unpanelled groups
# get one small panel per group, because nine overlapping ribbons say nothing.
lda_drc_plot <- function(lda, df, cols, title = NULL, subtitle = NULL,
                         x_step = NULL, ribbon = TRUE, panels = "auto") {
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

  pan <- use_panels(g, panels)
  facet <- NULL
  key <- "Group"; key_cols <- cols[g]
  if (!is.null(pan)) {
    raw$Panel  <- pan$Panel[match(raw$Group, pan$Group)]
    pred$Panel <- pan$Panel[match(pred$Group, pan$Group)]
    # Same treatment, same colour in every panel? Then the legend needs only
    # the treatment names — the panel title already says which line it is.
    by_lab <- tapply(unname(cols[pan$Group]), pan$Label, function(v) length(unique(v)) == 1)
    if (all(by_lab)) {
      key <- "Label"
      raw$Label  <- pan$Label[match(raw$Group, pan$Group)]
      pred$Label <- pan$Label[match(pred$Group, pan$Group)]
      key_cols <- stats::setNames(unname(cols[pan$Group[match(levels(pan$Label), pan$Label)]]),
                                  levels(pan$Label))
    }
    facet <- ggplot2::facet_wrap(~ Panel, nrow = if (nlevels(pan$Panel) > 3) 2 else 1)
  } else if (length(g) > 6) {
    facet <- ggplot2::facet_wrap(~ Group)
  }

  p <- ggplot2::ggplot()
  if (ribbon)
    p <- p + ggplot2::geom_ribbon(
      data = pred, ggplot2::aes(CellCount, ymin = ci_lo, ymax = ci_hi, fill = .data[[key]], group = Group),
      alpha = 0.18, colour = NA) +
      ggplot2::geom_line(data = pred, ggplot2::aes(CellCount, ci_lo, colour = .data[[key]], group = Group),
                         linewidth = 0.5, linetype = "dashed", alpha = 0.85) +
      ggplot2::geom_line(data = pred, ggplot2::aes(CellCount, ci_hi, colour = .data[[key]], group = Group),
                         linewidth = 0.5, linetype = "dashed", alpha = 0.85)

  n_leg <- length(key_cols)
  p <- p +
    ggplot2::geom_line(data = pred, ggplot2::aes(CellCount, fit, colour = .data[[key]], group = Group),
                       linewidth = 1.5) +
    ggplot2::geom_point(data = raw, ggplot2::aes(CellCount, Organoids, colour = .data[[key]], group = Group),
                        shape = 16, size = 3.4) +
    ggplot2::scale_colour_manual(values = key_cols,
      guide = if (!is.null(facet) && is.null(pan)) "none"
              else ggplot2::guide_legend(nrow = if (n_leg > 3) ceiling(n_leg / 3) else 1)) +
    ggplot2::scale_fill_manual(values = key_cols, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, x_max),
                                breaks = seq(0, x_max, by = x_step),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0, y_max),
                                breaks = scales::pretty_breaks(6),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(x = "Cells seeded", y = "Organoids formed",
                  title = title, subtitle = subtitle) +
    theme_lab()
  if (!is.null(facet))
    p <- p + facet + ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = if (length(g) > 9 && is.null(pan)) 10 else 14,
                                         face = "bold"),
      legend.text = ggplot2::element_text(size = if (n_leg > 6) 10 else 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 11),
      panel.spacing = ggplot2::unit(1, "lines"))
  p
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
    ggplot2::scale_colour_manual(values = cols[g],
      guide = ggplot2::guide_legend(nrow = if (length(g) > 3) ceiling(length(g) / 3) else 1)) +
    ggplot2::scale_fill_manual(values = cols[g], guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, x_max),
                                expand = ggplot2::expansion(mult = c(0, 0.01))) +
    ggplot2::labs(x = "Cells per well", y = "log fraction of negative wells",
                  title = title, subtitle = subtitle) +
    theme_lab()
}
