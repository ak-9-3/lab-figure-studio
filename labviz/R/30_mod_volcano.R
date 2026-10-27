# ── Module: volcano plot ──────────────────────────────────────────────────
# Takes a differential-expression results table (DESeq2 / edgeR / limma /
# anything with a fold change and a p-value) and draws a labelled volcano.

VOLC_SPEC <- list(
  Gene   = list(label = "Gene / feature", type = "character",
                aliases = c("gene", "genesymbol", "symbol", "id", "feature", "name", "x")),
  log2FC = list(label = "log2 fold change", type = "numeric",
                aliases = c("log2foldchange", "log2fc", "logfc", "lfc", "fc", "coef")),
  PValue = list(label = "p-value", type = "numeric",
                aliases = c("pvalue", "pval", "p", "praw", "pvalues")),
  padj   = list(label = "adjusted p-value", type = "numeric", optional = TRUE,
                aliases = c("padj", "fdr", "qvalue", "adjpval", "adjpvalue", "bh"))
)

VOLC_EXAMPLE <- local({
  set.seed(11)
  n <- 4000
  g <- sprintf("Gene%04d", seq_len(n))
  lfc <- stats::rnorm(n, 0, 0.6)
  hit <- sample(n, 160)
  lfc[hit] <- lfc[hit] + sample(c(-1, 1), 160, TRUE) * runif(160, 1.2, 4)
  p <- stats::runif(n)^3
  p[hit] <- p[hit] / 4000
  data.frame(Gene = g, log2FC = round(lfc, 3), PValue = p,
             padj = stats::p.adjust(p, "BH"), stringsAsFactors = FALSE)
})

volcano_prep <- function(df, use_adj, fc_cut, p_cut) {
  df$pval_used <- if (use_adj && "padj" %in% names(df)) df$padj else df$PValue
  df <- df[is.finite(df$log2FC) & is.finite(df$pval_used), , drop = FALSE]
  df$pval_used[df$pval_used <= 0] <- .Machine$double.xmin
  df$y <- -log10(df$pval_used)
  df$Status <- ifelse(df$pval_used < p_cut & df$log2FC >=  fc_cut, "Up",
               ifelse(df$pval_used < p_cut & df$log2FC <= -fc_cut, "Down", "NS"))
  df$Status <- factor(df$Status, levels = c("Down", "NS", "Up"))
  df
}

volcano_plot <- function(df, cols, fc_cut, p_cut, label_genes = character(),
                         title = NULL, subtitle = NULL, use_adj = TRUE,
                         pt_size = 1.6, alpha = 0.75, y_cap = NA,
                         show_counts = TRUE) {
  if (is.finite(y_cap) && y_cap > 0) {
    df$capped <- df$y > y_cap
    df$y <- pmin(df$y, y_cap)
  } else df$capped <- FALSE

  ylab <- if (use_adj) expression(-log[10]~italic(p)[adj]) else expression(-log[10]~italic(p))
  lim  <- max(abs(df$log2FC), na.rm = TRUE) * 1.05

  p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, y)) +
    ggplot2::geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed",
                        colour = "grey55", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed",
                        colour = "grey55", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(colour = Status, shape = capped),
                        size = pt_size, alpha = alpha) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    ggplot2::scale_colour_manual(values = cols, breaks = c("Down", "Up")) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim)) +
    ggplot2::labs(x = expression(log[2]~fold~change), y = ylab,
                  title = title, subtitle = subtitle) +
    theme_lab() +
    ggplot2::theme(legend.position = "top")

  if (show_counts) {
    n_up <- sum(df$Status == "Up"); n_dn <- sum(df$Status == "Down")
    ytop <- max(df$y, na.rm = TRUE)
    p <- p + ggplot2::annotate("text", x = c(-lim * 0.92, lim * 0.92), y = ytop,
                               label = c(paste0("▼ ", n_dn), paste0("▲ ", n_up)),
                               hjust = c(0, 1), vjust = 1, fontface = "bold",
                               size = 4.6, colour = c(cols[["Down"]], cols[["Up"]]))
  }

  lab <- df[df$Gene %in% label_genes, , drop = FALSE]
  if (nrow(lab)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, ggplot2::aes(label = Gene), size = 3.8, fontface = "bold",
        min.segment.length = 0, box.padding = 0.4, max.overlaps = Inf,
        segment.colour = "grey45", colour = "black", seed = 1)
    } else {
      p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(label = Gene),
                                  size = 3.6, vjust = -0.8, fontface = "bold")
    }
  }
  p
}

volcanoSidebarUI <- function(id) {
  ns <- NS(id)
  accordion(
    open = c("1 · Data", "2 · Thresholds"),
    accordion_panel("1 · Data", icon = bsicons::bs_icon("file-earmark-spreadsheet"),
      dataInputUI(ns("d"), VOLC_SPEC, "Example DE results"),
      div(class = "form-text", "Any DESeq2 / edgeR / limma results table works.")),
    accordion_panel("2 · Thresholds", icon = bsicons::bs_icon("rulers"),
      checkboxInput(ns("use_adj"), "Use adjusted p-value", TRUE),
      numericInput(ns("p_cut"), "p cutoff", 0.05, 0, 1, 0.01),
      numericInput(ns("fc_cut"), "|log2 FC| cutoff", 1, 0, 20, 0.25)),
    accordion_panel("3 · Labels", icon = bsicons::bs_icon("tag"),
      radioButtons(ns("label_mode"), NULL,
        choiceNames  = c("None", "Top by p-value", "Top by fold change", "My gene list"),
        choiceValues = c("none", "p", "fc", "custom"), selected = "p"),
      conditionalPanel(sprintf("input['%s'] != 'none' && input['%s'] != 'custom'",
                               ns("label_mode"), ns("label_mode")),
        numericInput(ns("top_n"), "How many", 15, 1, 100, 1)),
      conditionalPanel(sprintf("input['%s'] == 'custom'", ns("label_mode")),
        textAreaInput(ns("gene_list"), NULL, rows = 5,
                      placeholder = "One gene per line (or comma separated)"))),
    accordion_panel("4 · Figure", icon = bsicons::bs_icon("sliders"),
      textInput(ns("title"), "Title", "Differential expression"),
      textInput(ns("sub"), "Subtitle", ""),
      colorInput(ns("col_up"),   "Up",   "#D55E00"),
      colorInput(ns("col_down"), "Down", "#0072B2"),
      colorInput(ns("col_ns"),   "Not significant", "#C8C8C8"),
      sliderInput(ns("pt"), "Point size", 0.4, 4, 1.6, 0.1),
      sliderInput(ns("alpha"), "Opacity", 0.1, 1, 0.75, 0.05),
      numericInput(ns("ycap"), "Cap y axis at (0 = off)", 0, 0, 400, 5),
      checkboxInput(ns("counts"), "Show up/down counts", TRUE),
      sliderInput(ns("w"), "Width (in)", 3, 16, 7, 0.5),
      sliderInput(ns("h"), "Height (in)", 3, 16, 6.5, 0.5),
      sliderInput(ns("dpi"), "PNG resolution", 100, 600, 450, 50))
  )
}

volcanoUI <- function(id) {
  ns <- NS(id)
  tagList(
    card(full_screen = TRUE, card_header("Volcano"),
      plotOutput(ns("plot"), height = "620px"),
      card_footer(class = "d-flex gap-2 flex-wrap",
        downloadButton(ns("pdf"), "PDF", class = "btn-sm btn-outline-primary"),
        downloadButton(ns("svg"), "SVG", class = "btn-sm btn-outline-primary"),
        downloadButton(ns("png"), "PNG", class = "btn-sm btn-outline-primary"))),
    card(card_header("Significant hits"), uiOutput(ns("tbl_note")), div(class = "table-responsive", style = "max-height:480px; overflow-y:auto;", tableOutput(ns("tbl"))),
      card_footer(downloadButton(ns("csv"), "Download hits (CSV)",
                                 class = "btn-sm btn-outline-primary"))))
}

volcanoServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    dat <- dataInputServer("d", VOLC_SPEC, VOLC_EXAMPLE)

    prepped <- reactive({
      req(input$p_cut, input$fc_cut)
      volcano_prep(dat(), isTRUE(input$use_adj), input$fc_cut, input$p_cut)
    })

    labels <- reactive({
      d <- prepped(); sig <- d[d$Status != "NS", , drop = FALSE]
      switch(input$label_mode,
        none = character(),
        p    = utils::head(sig$Gene[order(sig$pval_used)], input$top_n),
        fc   = utils::head(sig$Gene[order(-abs(sig$log2FC))], input$top_n),
        custom = {
          g <- unlist(strsplit(input$gene_list %||% "", "[,\n\r\t ]+"))
          g <- trimws(g); g[nzchar(g)]
        })
    })

    cols <- reactive(c(Down = input$col_down, NS = input$col_ns, Up = input$col_up))

    plot_obj <- reactive({
      cap <- if (is.null(input$ycap) || input$ycap <= 0) NA else input$ycap
      volcano_plot(prepped(), cols(), input$fc_cut, input$p_cut, labels(),
                   input$title, input$sub, isTRUE(input$use_adj),
                   input$pt, input$alpha, cap, isTRUE(input$counts))
    })

    output$plot <- renderPlot(plot_obj(), res = 96)

    hits <- reactive({
      d <- prepped(); d <- d[d$Status != "NS", , drop = FALSE]
      d <- d[order(d$pval_used), ]
      keep <- intersect(c("Gene", "log2FC", "PValue", "padj", "Status"), names(d))
      d[, keep, drop = FALSE]
    })
    output$tbl <- renderTable({
      d <- hits()
      validate(need(nrow(d) > 0, "Nothing passes the current thresholds."))
      for (k in intersect(c("log2FC", "PValue", "padj"), names(d)))
        d[[k]] <- signif(d[[k]], 3)
      as_display(utils::head(d, 100))
    }, striped = TRUE, spacing = "s", width = "100%", align = "l")
    output$tbl_note <- renderUI({
      n <- nrow(hits()); if (n <= 100) return(NULL)
      div(class = "form-text mb-2",
          sprintf("Showing the top 100 of %d hits by p-value. The CSV has all of them.", n))
    })
    output$csv <- downloadHandler("volcano_hits.csv",
      function(f) utils::write.csv(hits(), f, row.names = FALSE))
    dl <- function(ext) downloadHandler(paste0("volcano.", ext),
      function(f) save_plot(f, plot_obj(), input$w, input$h, input$dpi))
    output$pdf <- dl("pdf"); output$svg <- dl("svg"); output$png <- dl("png")
  })
}
