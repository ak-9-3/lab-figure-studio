# ── Module: RT-qPCR (ΔΔCt) ────────────────────────────────────────────────
# Technical replicates are averaged first; statistics are run on ΔCt values
# (never on fold changes — those are log-normal and t-tests on them are wrong).

QPCR_SPEC <- list(
  Gene      = list(label = "Gene / target", type = "character",
                   aliases = c("gene", "target", "assay", "primer", "detector")),
  Group     = list(label = "Group / condition", type = "character",
                   aliases = c("group", "condition", "treatment", "sample", "genotype")),
  Replicate = list(label = "Biological replicate", type = "character",
                   aliases = c("replicate", "biorep", "rep", "experiment", "batch")),
  Ct        = list(label = "Ct / Cq value", type = "numeric",
                   aliases = c("ct", "cq", "cycle", "ctmean", "cqmean"))
)

QPCR_EXAMPLE <- local({
  set.seed(7)
  grid <- expand.grid(Gene = c("Actb", "Lgr5", "Olfm4", "Axin2"),
                      Group = c("Control", "Treated"),
                      Replicate = c("R1", "R2", "R3"),
                      Tech = 1:2, stringsAsFactors = FALSE)
  base <- c(Actb = 18, Lgr5 = 26, Olfm4 = 25, Axin2 = 27)
  shift <- c(Actb = 0, Lgr5 = -1.6, Olfm4 = -2.1, Axin2 = 0.2)
  grid$Ct <- round(base[grid$Gene] +
                   ifelse(grid$Group == "Treated", shift[grid$Gene], 0) +
                   stats::rnorm(nrow(grid), 0, 0.25), 2)
  grid[, c("Gene", "Group", "Replicate", "Ct")]
})

qpcr_analyse <- function(df, hk, ref, targets) {
  # 1. collapse technical replicates
  agg <- stats::aggregate(Ct ~ Gene + Group + Replicate, data = df, FUN = mean)
  # 2. ΔCt against the housekeeping gene, within each Group × Replicate
  hkv <- agg[agg$Gene == hk, c("Group", "Replicate", "Ct")]
  names(hkv)[3] <- "Ct_hk"
  m <- merge(agg[agg$Gene %in% targets, ], hkv, by = c("Group", "Replicate"))
  if (!nrow(m))
    stop("No rows left — is the housekeeping gene measured in every replicate?")
  m$dCt <- m$Ct - m$Ct_hk
  # 3. ΔΔCt against the reference group, per gene
  refm <- stats::aggregate(dCt ~ Gene, data = m[m$Group == ref, ], FUN = mean)
  names(refm)[2] <- "dCt_ref"
  m <- merge(m, refm, by = "Gene")
  m$ddCt <- m$dCt - m$dCt_ref
  m$Fold <- 2^(-m$ddCt)
  m
}

qpcr_summary <- function(m) {
  s <- do.call(rbind, lapply(split(m, list(m$Gene, m$Group), drop = TRUE), function(d) {
    n <- nrow(d)
    data.frame(Gene = d$Gene[1], Group = d$Group[1], n = n,
               Fold = mean(d$Fold),
               SEM  = if (n > 1) stats::sd(d$Fold) / sqrt(n) else 0,
               stringsAsFactors = FALSE)
  }))
  rownames(s) <- NULL
  s
}

qpcr_tests <- function(m, ref) {
  groups <- setdiff(unique(m$Group), ref)
  if (!length(groups)) return(NULL)
  out <- do.call(rbind, lapply(unique(m$Gene), function(g) {
    do.call(rbind, lapply(groups, function(gr) {
      a <- m$dCt[m$Gene == g & m$Group == ref]
      b <- m$dCt[m$Gene == g & m$Group == gr]
      p <- if (length(a) > 1 && length(b) > 1)
             tryCatch(stats::t.test(b, a)$p.value, error = function(e) NA_real_) else NA_real_
      data.frame(Gene = g, A = ref, B = gr, Fold = 2^(-(mean(b) - mean(a))),
                 p = p, stringsAsFactors = FALSE)
    }))
  }))
  out$p_adj <- stats::p.adjust(out$p, "BH")
  out$stars <- p_stars(out$p_adj)
  out
}

qpcr_plot <- function(s, m, tests, cols, ref, title, subtitle, ylab,
                      log_y = FALSE, show_dots = TRUE, ncol = NULL) {
  s$Group <- factor(s$Group, levels = names(cols))
  m$Group <- factor(m$Group, levels = names(cols))
  p <- ggplot2::ggplot(s, ggplot2::aes(Group, Fold)) +
    ggplot2::geom_col(ggplot2::aes(fill = Group, colour = Group),
                      width = 0.6, linewidth = 1.1) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(Fold - SEM, 0), ymax = Fold + SEM),
                           width = 0.14, linewidth = 0.8, colour = "black")
  if (show_dots)
    p <- p + ggplot2::geom_point(data = m, ggplot2::aes(Group, Fold),
                                 inherit.aes = FALSE, shape = 21, size = 2.6,
                                 stroke = 1, colour = "grey25", fill = NA,
                                 position = ggplot2::position_jitter(0.07, 0, seed = 3))
  if (!is.null(tests) && nrow(tests)) {
    # stars sit above whichever is taller: the error bar or the highest dot
    a <- stats::aggregate(cbind(hi = Fold + SEM) ~ Gene, data = s, FUN = max)
    b <- stats::aggregate(cbind(hi = Fold)       ~ Gene, data = m, FUN = max)
    tops <- merge(a, b, by = "Gene")
    tops$top <- pmax(tops$hi.x, tops$hi.y)
    lab <- merge(tests, tops[, c("Gene", "top")], by = "Gene")
    p <- p + ggplot2::geom_text(
      data = lab, ggplot2::aes(x = B, y = top * 1.06, label = stars),
      inherit.aes = FALSE, size = 5, fontface = "bold", vjust = 0)
  }
  p <- p +
    ggplot2::facet_wrap(~ Gene, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = lighten(cols, 0.87), guide = "none") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::labs(x = NULL, y = ylab, title = title, subtitle = subtitle) +
    theme_lab(14) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold.italic", size = 15),
                   panel.spacing = ggplot2::unit(1.1, "lines"))
  if (log_y)
    p + ggplot2::scale_y_continuous(trans = "log2", expand = ggplot2::expansion(mult = c(0.05, 0.15)))
  else
    p + ggplot2::scale_y_continuous(limits = c(0, NA),
                                    expand = ggplot2::expansion(mult = c(0, 0.20)))
}

qpcrSidebarUI <- function(id) {
  ns <- NS(id)
  accordion(
    open = c("1 · Data", "2 · Design"),
    accordion_panel("1 · Data", icon = bsicons::bs_icon("file-earmark-spreadsheet"),
      dataInputUI(ns("d"), QPCR_SPEC),
      div(class = "form-text",
          "Long format: one row per well. Technical replicates are averaged automatically.")),
    accordion_panel("2 · Design", icon = bsicons::bs_icon("diagram-3"),
      selectInput(ns("hk"), "Housekeeping gene", choices = NULL),
      selectInput(ns("ref"), "Reference group", choices = NULL),
      selectizeInput(ns("targets"), "Target genes", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button", "drag_drop"))),
      selectizeInput(ns("groups"), "Groups (drag to order)", choices = NULL,
                     multiple = TRUE,
                     options = list(plugins = list("remove_button", "drag_drop")))),
    accordion_panel("3 · Figure", icon = bsicons::bs_icon("sliders"),
      textInput(ns("title"), "Title", "RT-qPCR"),
      textInput(ns("sub"), "Subtitle", "mean ± SEM  -  dots = biological replicates"),
      textInput(ns("ylab"), "Y axis label", "Relative expression (2^-ΔΔCt)"),
      selectInput(ns("palette"), "Palette", choices = names(PALETTES)),
      uiOutput(ns("colour_ui")),
      checkboxInput(ns("dots"), "Show replicate points", TRUE),
      checkboxInput(ns("logy"), "Log2 y-axis", FALSE),
      numericInput(ns("ncol"), "Facet columns", 2, 1, 8, 1),
      sliderInput(ns("w"), "Width (in)", 3, 16, 8, 0.5),
      sliderInput(ns("h"), "Height (in)", 3, 16, 6, 0.5),
      sliderInput(ns("dpi"), "PNG resolution", 100, 600, 450, 50))
  )
}

qpcrUI <- function(id) {
  ns <- NS(id)
  tagList(
    card(full_screen = TRUE, card_header("Relative expression"),
      plotOutput(ns("plot"), height = "600px"),
      card_footer(class = "d-flex gap-2 flex-wrap",
        downloadButton(ns("pdf"), "PDF", class = "btn-sm btn-outline-primary"),
        downloadButton(ns("svg"), "SVG", class = "btn-sm btn-outline-primary"),
        downloadButton(ns("png"), "PNG", class = "btn-sm btn-outline-primary"))),
    layout_columns(col_widths = c(7, 5),
      card(card_header("Per-replicate values"), DTOutput(ns("tbl")),
           card_footer(downloadButton(ns("csv"), "CSV", class = "btn-sm btn-outline-primary"))),
      card(card_header("Welch t-test on ΔCt (BH adjusted)"), DTOutput(ns("stats")))))
}

qpcrServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    dat <- dataInputServer("d", QPCR_SPEC, QPCR_EXAMPLE)

    genes  <- reactive(sort(unique(dat()$Gene)))
    groups <- reactive(unique(as.character(dat()$Group)))

    observeEvent(genes(), {
      g  <- genes()
      hk <- g[grepl("actb|gapdh|b2m|hprt|tbp|rplp0|18s", tolower(g))]
      updateSelectInput(session, "hk", choices = g,
                        selected = if (length(hk)) hk[1] else g[1])
    }, ignoreNULL = TRUE)

    observeEvent(list(genes(), input$hk), {
      tg <- setdiff(genes(), input$hk %||% "")
      updateSelectizeInput(session, "targets", choices = tg, selected = tg,
                           options = list(plugins = list("remove_button", "drag_drop")))
    }, ignoreNULL = TRUE)

    observeEvent(groups(), {
      g <- groups()
      updateSelectInput(session, "ref", choices = g, selected = g[1])
      updateSelectizeInput(session, "groups", choices = g, selected = g,
                           options = list(plugins = list("remove_button", "drag_drop")))
    }, ignoreNULL = TRUE)

    output$colour_ui <- renderUI({
      g <- req(input$groups); d <- pal_for(g, input$palette)
      lapply(g, function(k)
        colourInput(ns(paste0("col_", make.names(k))), k, value = unname(d[k])))
    })
    cols <- reactive({
      g <- req(input$groups); d <- pal_for(g, input$palette)
      stats::setNames(vapply(g, function(k) {
        x <- input[[paste0("col_", make.names(k))]]
        if (is.null(x) || !nzchar(x)) unname(d[k]) else x }, character(1)), g)
    })

    outputOptions(output, "colour_ui", suspendWhenHidden = FALSE)

    m <- reactive({
      req(input$hk, input$ref, input$targets, input$groups)
      d <- dat(); d <- d[d$Group %in% input$groups | d$Gene == input$hk, ]
      out <- qpcr_analyse(d, input$hk, input$ref, input$targets)
      out <- out[out$Group %in% input$groups, ]
      out$Gene <- factor(out$Gene, levels = input$targets)
      out
    })
    summ  <- reactive(qpcr_summary(m()))
    tests <- reactive(qpcr_tests(m(), input$ref))

    plot_obj <- reactive(
      qpcr_plot(summ(), m(), tests(), cols(), input$ref, input$title, input$sub,
                input$ylab, input$logy, input$dots, input$ncol))

    output$plot <- renderPlot(plot_obj(), res = 96)
    output$tbl <- renderDT({
      d <- m()[, c("Gene", "Group", "Replicate", "Ct", "Ct_hk", "dCt", "ddCt", "Fold")]
      d[, 4:8] <- lapply(d[, 4:8], function(x) round(x, 3))
      datatable(d, rownames = FALSE, class = "compact stripe",
                options = list(pageLength = 10, scrollX = TRUE))
    })
    output$stats <- renderDT({
      t <- tests(); validate(need(!is.null(t), "Only one group — nothing to test."))
      d <- data.frame(Gene = t$Gene, Comparison = paste(t$B, "vs", t$A),
                      Fold = signif(t$Fold, 3), p = fmt_p(t$p),
                      `p adj` = fmt_p(t$p_adj), sig = t$stars, check.names = FALSE)
      datatable(d, rownames = FALSE, class = "compact stripe",
                options = list(dom = "t", paging = FALSE, scrollX = TRUE))
    })
    output$csv <- downloadHandler("qPCR_ddCt.csv",
      function(f) utils::write.csv(m(), f, row.names = FALSE))
    dl <- function(ext) downloadHandler(paste0("qPCR.", ext),
      function(f) save_plot(f, plot_obj(), input$w, input$h, input$dpi))
    output$pdf <- dl("pdf"); output$svg <- dl("svg"); output$png <- dl("png")
  })
}
