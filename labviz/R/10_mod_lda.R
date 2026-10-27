# ── Module: limiting-dilution assay ───────────────────────────────────────

LDA_SPEC_CELLS <- list(
  # Many sheets have two seeding columns: what the protocol said (100, 200…)
  # and what was actually counted into the well (87, 179…). The counted one
  # goes in the model, so its names come first.
  CellCount = list(label = "Cells seeded (the counted number, if you have both)", type = "numeric",
                   aliases = c("cellscounted", "counted", "actualcells", "actual", "cellcount",
                               "cells", "cellsseeded", "seeded", "input", "n")),
  Organoids = list(label = "Organoids formed", type = "numeric",
                   aliases = c("organoid", "colony", "colonies", "sphere", "spheroid",
                               "formed", "positive", "count")),
  Group     = list(label = "Group / condition", type = "character",
                   aliases = c("group", "condition", "treatment", "sample", "sgrna", "guide",
                               "genotype", "name")),
  NominalCells = list(label = "Protocol seeding number (only used to check the counted one)",
                      type = "numeric", optional = TRUE,
                      aliases = c("cellsseeded", "seeded", "seeding", "protocol", "nominal",
                                  "target", "planned", "intended")),
  Replicate = list(label = "Replicate", type = "character", optional = TRUE,
                   aliases = c("replicate", "rep", "biorep")),
  Batch     = list(label = "Extra column (cell line, screen, timepoint)", type = "character", optional = TRUE,
                   aliases = c("experiment", "screen", "batch", "plate", "panel",
                               "cellline", "celline", "line", "date", "target"))
)

LDA_SPEC_WELLS <- list(
  Dose     = list(label = "Cells per well", type = "numeric",
                  aliases = c("dose", "cellsperwell", "cells", "seeded")),
  Tested   = list(label = "Wells tested", type = "numeric",
                  aliases = c("tested", "wells", "total", "n")),
  Positive = list(label = "Wells positive", type = "numeric",
                  aliases = c("positive", "pos", "responding", "growth")),
  Group    = list(label = "Group / condition", type = "character",
                  aliases = c("group", "condition", "treatment", "sgrna", "genotype")),
  Batch    = list(label = "Extra column (cell line, screen, timepoint)", type = "character", optional = TRUE,
                  aliases = c("experiment", "screen", "batch", "plate", "panel",
                              "cellline", "celline", "line", "date", "target"))
)

LDA_EXAMPLE_CELLS <- data.frame(
  CellCount = c(137, 289, 512, 668,  71, 152, 316, 534,
                 68, 143, 251, 321, 145, 287, 520, 690),
  Organoids = c( 15,  33,  57,  74,   3,   5,  11,  19,
                  3,   5,  10,  13,  17,  31,  59,  77),
  Group     = rep(c("Cpox", "Alas1", "NTC", "Cpox"), each = 4),
  Replicate = rep(c("R1", "R1", "R1", "R2"), each = 4),
  Batch     = rep("Screen 1", 16),
  stringsAsFactors = FALSE
)

LDA_EXAMPLE_WELLS <- data.frame(
  Dose     = rep(c(10, 30, 100, 300), 3),
  Tested   = rep(24, 12),
  Positive = c(3, 8, 17, 23,  1, 3, 9, 18,  6, 14, 22, 24),
  Group    = rep(c("Control", "KO", "Rescue"), each = 4),
  stringsAsFactors = FALSE
)

# ---- sidebar --------------------------------------------------------------
ldaSidebarUI <- function(id) {
  ns <- NS(id)
  accordion(
    open = c("1 · Data", "2 · Groups & colours"),
    accordion_panel(
      "1 · Data", icon = bsicons::bs_icon("file-earmark-spreadsheet"),
      radioButtons(ns("mode"), "Assay layout", inline = FALSE,
        choiceNames = list(
          HTML("<b>Counted organoids</b><br><span class='text-muted small'>cells seeded &rarr; organoids counted</span>"),
          HTML("<b>Scored wells</b><br><span class='text-muted small'>wells +/- at each dose (classic ELDA)</span>")),
        choiceValues = c("cells", "wells")),
      conditionalPanel(sprintf("input['%s'] == 'cells'", ns("mode")),
                       dataInputUI(ns("d_cells"), LDA_SPEC_CELLS, protocol = TRUE)),
      conditionalPanel(sprintf("input['%s'] == 'wells'", ns("mode")),
                       dataInputUI(ns("d_wells"), LDA_SPEC_WELLS)),
      hr(),
      selectInput(ns("batch_by"), "Split figures by", choices = NULL),
      div(class = "form-text",
          "One figure per level of this column, all in a single ZIP.")
    ),
    accordion_panel(
      "2 · Groups & colours", icon = bsicons::bs_icon("palette"),
      selectizeInput(ns("groups"), "Groups to plot", choices = NULL, multiple = TRUE,
                     options = list(plugins = list("remove_button", "drag_drop"))),
      div(class = "form-text mb-2", "Drag to reorder — order sets the x-axis."),
      selectInput(ns("palette"), "Palette", choices = names(PALETTES)),
      uiOutput(ns("colour_ui"))
    ),
    accordion_panel(
      "3 · Statistics", icon = bsicons::bs_icon("calculator"),
      selectInput(ns("ref"), "Reference group", choices = NULL),
      radioButtons(ns("comp"), "Comparisons", inline = TRUE,
                   choiceNames  = c("none", "vs ref", "all pairs"),
                   choiceValues = c("none", "ref", "all"), selected = "ref"),
      selectInput(ns("adjust"), "Multiple-testing correction",
                  choices = c("none", "BH", "holm", "bonferroni"), selected = "none"),
      sliderInput(ns("conf"), "Confidence level", 0.80, 0.99, 0.95, 0.01)
    ),
    accordion_panel(
      "4 · Figure", icon = bsicons::bs_icon("sliders"),
      textInput(ns("bar_title"), "Bar plot title", "Stem cell frequency"),
      textInput(ns("bar_sub"),   "Bar plot subtitle",
                "95% Wald CI  -  dots = observed rate per dose"),
      textInput(ns("drc_title"), "Curve title", "Dose–response"),
      textInput(ns("drc_sub"),   "Curve subtitle", "shading = 95% CI"),
      textInput(ns("ylab"), "Frequency axis label", "Organoid-forming frequency"),
      checkboxInput(ns("dots"),   "Show observed points", TRUE),
      checkboxInput(ns("labels"), "Show 1/N labels", TRUE),
      checkboxInput(ns("ribbon"), "Show CI ribbon on curve", TRUE),
      numericInput(ns("x_step"), "Curve x-axis tick step", 150, 1, 10000, 10),
      sliderInput(ns("w"), "Width (in)",  3, 14, 5.5, 0.5),
      sliderInput(ns("h"), "Height (in)", 3, 14, 6.0, 0.5),
      sliderInput(ns("dpi"), "PNG resolution", 100, 600, 450, 50)
    )
  )
}

# ---- main panel -----------------------------------------------------------
ldaUI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("nominal_warning")),
    uiOutput(ns("showing_bar")),
    layout_columns(
      id = ns("figs"), col_widths = c(6, 6),
      card(full_screen = TRUE,
        card_header("Stem cell frequency"),
        plotOutput(ns("bar"), height = "520px"),
        card_footer(class = "d-flex gap-2 flex-wrap",
          downloadButton(ns("bar_pdf"), "PDF", class = "btn-sm btn-outline-primary"),
          downloadButton(ns("bar_svg"), "SVG", class = "btn-sm btn-outline-primary"),
          downloadButton(ns("bar_png"), "PNG", class = "btn-sm btn-outline-primary"))),
      card(full_screen = TRUE,
        card_header(textOutput(ns("curve_header"), inline = TRUE)),
        plotOutput(ns("curve"), height = "520px"),
        card_footer(class = "d-flex gap-2 flex-wrap",
          downloadButton(ns("crv_pdf"), "PDF", class = "btn-sm btn-outline-primary"),
          downloadButton(ns("crv_svg"), "SVG", class = "btn-sm btn-outline-primary"),
          downloadButton(ns("crv_png"), "PNG", class = "btn-sm btn-outline-primary")))
    ),
    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Frequencies"), div(class = "table-responsive", tableOutput(ns("freq"))),
           card_footer(downloadButton(ns("freq_csv"), "CSV",
                                      class = "btn-sm btn-outline-primary"))),
      card(card_header("Group comparisons"), div(class = "table-responsive", tableOutput(ns("comp_tbl"))),
           uiOutput(ns("gof")),
           card_footer(downloadButton(ns("comp_csv"), "CSV",
                                      class = "btn-sm btn-outline-primary")))
    ),
    card(
      card_header("Take it with you"),
      card_body(
        p(class = "mb-2",
          "Download a ZIP with your data, the analysis code, and a script that ",
          "rebuilds these exact figures. Open it in RStudio and edit anything."),
        div(class = "d-flex gap-2 flex-wrap align-items-center",
          downloadButton(ns("bundle"), "Download data + R code (.zip)",
                         class = "btn-primary"),
          uiOutput(ns("batch_button"), inline = TRUE)),
        uiOutput(ns("batch_note")))
    )
  )
}

# ---- server ---------------------------------------------------------------
ldaServer <- function(id, helper_dir) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    d_cells <- dataInputServer("d_cells", LDA_SPEC_CELLS, LDA_EXAMPLE_CELLS,
                               protocol = TRUE)
    d_wells <- dataInputServer("d_wells", LDA_SPEC_WELLS, LDA_EXAMPLE_WELLS)
    dat <- reactive(if (identical(input$mode, "wells")) d_wells() else d_cells())

    all_groups <- reactive(unique(as.character(dat()$Group)))

    # ── Batching ──────────────────────────────────────────────────────────
    # Any character column other than Group can split the data into separate
    # figures — one per cell line, screen, timepoint, whatever they called it.
    batch_choices <- reactive({
      d <- dat()
      chr <- names(d)[vapply(d, is.character, logical(1))]
      setdiff(chr, "Group")
    })
    # A sheet holding several experiments must not be pooled: comparing a
    # knockdown against a drug control from a different experiment is
    # meaningless. Split by default whenever the experiment column varies.
    observeEvent(batch_choices(), {
      ch <- batch_choices(); d <- dat()
      auto <- if ("Batch" %in% ch && length(unique(d$Batch)) >= 2) "Batch" else ""
      updateSelectInput(session, "batch_by",
                        choices = c("\u2014 don't split \u2014" = "", ch),
                        selected = auto)
    }, ignoreNULL = FALSE)

    batch_col <- reactive({
      b <- input$batch_by
      if (is.null(b) || !nzchar(b) || !b %in% names(dat())) NULL else b
    })
    batch_levels <- reactive({
      b <- batch_col(); if (is.null(b)) return(NULL)
      unique(as.character(dat()[[b]]))
    })
    observeEvent(batch_levels(), {
      lv <- batch_levels()
      updateSelectInput(session, "batch_show", choices = lv, selected = lv[1])
    }, ignoreNULL = TRUE)

    # If the sheet carried both a protocol dose and a real count, show what
    # choosing the wrong one would do. It is never a rounding error.
    output$nominal_warning <- renderUI({
      d <- tryCatch(cur_data(), error = function(e) NULL)
      if (is.null(d) || !"NominalCells" %in% names(d)) return(NULL)
      cmp <- tryCatch(nominal_cost(d, sel_groups(), input$mode),
                      error = function(e) NULL)
      if (is.null(cmp)) return(NULL)
      card(class = "mb-3 border-warning",
        card_body(class = "py-2",
          div(class = "fw-bold mb-1",
              bsicons::bs_icon("exclamation-triangle-fill"),
              " This sheet had two seeding columns"),
          div(class = "small text-muted mb-2",
              "Fitting the protocol numbers instead of the cells you counted ",
              "would move every frequency:"),
          renderTable(cmp, width = "100%", striped = TRUE, spacing = "xs", digits = 4)))
    })
    outputOptions(output, "nominal_warning", suspendWhenHidden = FALSE)

    output$showing_bar <- renderUI({
      lv <- batch_levels()
      if (is.null(lv) || length(lv) < 2) return(NULL)
      div(class = "card mb-3", div(class = "card-body py-2 d-flex align-items-center gap-3",
        span(class = "fw-bold", bsicons::bs_icon("collection"), " Showing"),
        div(style = "min-width:220px",
            selectInput(ns("batch_show"), NULL, choices = lv, selected = lv[1],
                        width = "100%")),
        span(class = "text-muted small",
             sprintf("Your data holds %d separate experiments, so each gets its own figure. All %d go into the batch ZIP.",
                     length(lv), length(lv)))))
    })
    outputOptions(output, "showing_bar", suspendWhenHidden = FALSE)

    # Everything downstream works on the slice currently being shown.
    cur_data <- reactive({
      b <- batch_col(); d <- dat()
      if (is.null(b)) return(d)
      lv <- input$batch_show %||% batch_levels()[1]
      d[as.character(d[[b]]) == lv, , drop = FALSE]
    })

    # Static inputs updated from the data, rather than renderUI — dynamic UI
    # inside a collapsed accordion has a start-up race that leaves the app
    # blank until you touch a control.
    observeEvent(all_groups(), {
      g <- all_groups()
      updateSelectizeInput(session, "groups", choices = g, selected = g,
                           options = list(plugins = list("remove_button", "drag_drop")))
      updateSelectInput(session, "ref", choices = g, selected = g[1])
    }, ignoreNULL = TRUE)

    sel_groups <- reactive({
      g <- input$groups
      validate(need(length(g) >= 1, "Pick at least one group."))
      g
    })

    # A batch level may not contain every group; fit only what is there.
    fit_groups <- reactive({
      g <- intersect(sel_groups(), unique(as.character(cur_data()$Group)))
      validate(need(length(g) >= 1, "No selected group appears in this slice."))
      g
    })

    output$colour_ui <- renderUI({
      g <- sel_groups(); d <- lda_default_cols(g, input$palette)
      lapply(g, function(k)
        colorInput(ns(paste0("col_", make.names(k))), k, value = unname(d[k])))
    })

    cols <- reactive({
      g <- sel_groups(); d <- lda_default_cols(g, input$palette)
      v <- vapply(g, function(k) {
        x <- input[[paste0("col_", make.names(k))]]
        if (is.null(x) || !nzchar(x)) unname(d[k]) else x
      }, character(1))
      stats::setNames(v, g)
    })

    # More bars need a wider figure; widen the download to fit (the slider
    # still lets them pick anything).
    observeEvent(sel_groups(), {
      g <- sel_groups(); pan <- group_panels(g)
      slots <- length(g) + if (is.null(pan)) 0 else nlevels(pan$Panel) - 1
      w <- if (slots <= 6) 5.5 else min(14, round((1.8 + 0.55 * slots) * 2) / 2)
      updateSliderInput(session, "w", value = w)
      # Many bars: give each figure the full page width instead of half.
      session$sendCustomMessage("lab-wide", list(id = ns("figs"), wide = slots > 8))
    })

    # accordion panels start collapsed; without this their inputs never exist
    outputOptions(output, "colour_ui", suspendWhenHidden = FALSE)

    lda <- reactive({
      validate(need(nrow(cur_data()) >= 2, "Need at least two rows of data."))
      # Turn a fitting failure into something a person can act on, rather than
      # letting a raw R error land in all four output panels.
      out <- tryCatch(
        lda_fit(cur_data(), fit_groups(), mode = input$mode, conf = input$conf),
        error = function(e) conditionMessage(e))
      validate(need(!is.character(out), out))
      out
    })

    pair_list <- reactive({
      g <- fit_groups()
      lda_pairs(g, if (isTRUE(input$ref %in% g)) input$ref else g[1],
                input$comp %||% "ref")
    })

    comps <- reactive({
      p <- pair_list(); if (!length(p)) return(NULL)
      lda_pairwise(lda(), p, adjust = input$adjust)
    })

    bar_plot <- reactive({
      lda_bar_plot(lda(), cur_data(), cols(), comps(),
                   title = input$bar_title, subtitle = input$bar_sub,
                   show_dots = input$dots, show_labels = input$labels,
                   ylab = input$ylab)
    })

    curve_plot <- reactive({
      if (identical(input$mode, "wells"))
        lda_elda_plot(lda(), cur_data(), cols(), input$drc_title, input$drc_sub)
      else
        lda_drc_plot(lda(), cur_data(), cols(), input$drc_title, input$drc_sub,
                     x_step = input$x_step, ribbon = input$ribbon)
    })

    output$curve_header <- renderText(
      if (identical(input$mode, "wells")) "ELDA plot" else "Dose–response curve")
    output$bar   <- renderPlot(bar_plot(),   res = 96)
    output$curve <- renderPlot(curve_plot(), res = 96)

    output$freq <- renderTable({
      d <- lda()$est[, c("Group", "n_rows", "Frequency", "CI_low", "CI_high", "1 in")]
      d$Frequency <- signif(d$Frequency, 4)
      d$CI_low    <- signif(d$CI_low, 4)
      d$CI_high   <- signif(d$CI_high, 4)
      d$`1 in`    <- round(d$`1 in`, 1)
      as_display(d)
    }, striped = TRUE, spacing = "s", width = "100%", align = "l")

    output$comp_tbl <- renderTable({
      d <- comps()
      validate(need(!is.null(d), "No comparisons selected."))
      out <- data.frame(Comparison = paste(d$B, "vs", d$A),
                        `Fold change` = signif(d$Ratio, 3),
                        `95% CI` = sprintf("%.2f–%.2f", d$CI_low, d$CI_high),
                        p = fmt_p(d$p), `p adj` = fmt_p(d$p_adj), sig = d$stars,
                        check.names = FALSE)
      as_display(out)
    }, striped = TRUE, spacing = "s", width = "100%", align = "l")

    output$gof <- renderUI({
      l <- lda()
      if (is.na(l$gof_p) || l$df_resid < 1) return(NULL)
      ok <- l$gof_p > 0.05
      div(class = paste("alert py-1 px-2 m-2 small",
                        if (ok) "alert-light" else "alert-warning"),
          sprintf("Goodness of fit: deviance %.1f on %d df, p = %s. %s",
                  l$deviance, l$df_resid, fmt_p(l$gof_p),
                  if (ok) "Single-hit model fits."
                  else "Poor fit — doses may not scale linearly, or replicates disagree."))
    })

    output$freq_csv <- downloadHandler("LDA_frequencies.csv",
      function(f) utils::write.csv(lda()$est, f, row.names = FALSE))
    output$comp_csv <- downloadHandler("LDA_comparisons.csv",
      function(f) utils::write.csv(comps(), f, row.names = FALSE))

    dl <- function(name, plot_fun, ext)
      downloadHandler(paste0(name, ".", ext),
        function(f) save_plot(f, plot_fun(), input$w, input$h, input$dpi))
    output$bar_pdf <- dl("LDA_frequency", bar_plot,   "pdf")
    output$bar_svg <- dl("LDA_frequency", bar_plot,   "svg")
    output$bar_png <- dl("LDA_frequency", bar_plot,   "png")
    output$crv_pdf <- dl("LDA_curve",     curve_plot, "pdf")
    output$crv_svg <- dl("LDA_curve",     curve_plot, "svg")
    output$crv_png <- dl("LDA_curve",     curve_plot, "png")

    helper_files <- reactive(
      file.path(helper_dir, c("00_theme.R", "02_lda_stats.R", "03_lda_plots.R")))

    export_opts <- reactive(list(
      conf = input$conf, adjust = input$adjust, pairs = pair_list(),
      ref = input$ref %||% "", comp = input$comp %||% "ref",
      bar_title = input$bar_title, bar_sub = input$bar_sub,
      drc_title = input$drc_title, drc_sub = input$drc_sub, ylab = input$ylab,
      dots = input$dots, labels = input$labels,
      x_step = input$x_step, w = input$w, h = input$h))

    output$bundle <- downloadHandler("LDA_analysis.zip", function(f) {
      export_bundle(
        f, cur_data(),
        lda_script(input$mode, fit_groups(), cols(), export_opts()),
        helpers = helper_files())
    })

    output$batch_button <- renderUI({
      lv <- batch_levels()
      if (is.null(lv) || length(lv) < 2) return(NULL)
      downloadButton(ns("batch_zip"),
                     sprintf("Download all %d figures (.zip)", length(lv)),
                     class = "btn-outline-primary")
    })
    outputOptions(output, "batch_button", suspendWhenHidden = FALSE)

    output$batch_note <- renderUI({
      lv <- batch_levels()
      if (is.null(lv) || length(lv) < 2) return(NULL)
      div(class = "form-text mt-2",
          sprintf("Batch ZIP: %s \u2014 PDF and PNG for every level of '%s', ",
                  paste(lv, collapse = ", "), batch_col()),
          "with per-level and pooled tables.")
    })

    output$batch_zip <- downloadHandler("LDA_batch.zip", function(f) {
      withProgress(message = "Building figures", value = 0, {
        lda_batch_bundle(
          f, dat(), batch_col(), input$mode, sel_groups(), cols(), export_opts(),
          helpers = helper_files(),
          script  = lda_script(input$mode, fit_groups(), cols(), export_opts()),
          progress = function(frac, lv)
            setProgress(value = frac, detail = lv))
      })
    })
  })
}
