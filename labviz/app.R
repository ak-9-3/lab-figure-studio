# ═══════════════════════════════════════════════════════════════════════════
#  Lab Figure Studio
#  Upload raw assay data → publication figures → take the code home.
#
#  Run:  open this file in RStudio and click "Run App"
#        or  shiny::runApp("labviz")
#
#  Everything in R/ is auto-loaded by Shiny. Each assay is its own module,
#  so adding a new one means adding one file — nothing else changes.
# ═══════════════════════════════════════════════════════════════════════════

# In webR/shinylive every package is baked into the bundle already, and there
# is no CRAN to install from — so only offer to install when running locally.
if (!(identical(R.version$os, "emscripten") ||
      grepl("emscripten", R.version$platform, fixed = TRUE))) {
  REQUIRED <- c("shiny", "bslib", "bsicons", "ggplot2", "scales", "ggrepel",
                "zip", "readxl")
  missing <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing: ", paste(missing, collapse = ", "))
    utils::install.packages(missing)
  }
}
library(shiny); library(bslib); library(ggplot2)

# Shiny draws every on-screen plot with ragg whenever it is installed, and in
# the browser that costs a 13 MB Unicode download before the first figure.
# R's own png() looks the same here, so use it.
options(shiny.useragg = FALSE)

HELPER_DIR <- normalizePath("R", mustWork = FALSE)

theme_app <- bs_theme(
  version = 5, bootswatch = "flatly",
  primary = "#0072B2",
  # A system font stack rather than Google Fonts. font_google(local = TRUE)
  # cannot download inside webR, and local = FALSE would make every visitor's
  # browser call fonts.googleapis.com — the only third-party request the app
  # would make. Dropping it means the page talks to nobody but its own origin.
  base_font    = font_collection("Inter", "-apple-system", "BlinkMacSystemFont",
                                 "Segoe UI", "Roboto", "Helvetica Neue", "Arial",
                                 "sans-serif"),
  heading_font = font_collection("Inter", "-apple-system", "BlinkMacSystemFont",
                                 "Segoe UI", "Roboto", "Helvetica Neue", "Arial",
                                 "sans-serif"),
  font_scale = 0.95
)

ui <- page_navbar(
  id = "nav",
  title = span(bsicons::bs_icon("bar-chart-line-fill"), " Lab Figure Studio"),
  window_title = "Lab Figure Studio",
  theme = theme_app, fillable = FALSE,
  header = tags$head(includeCSS("www/custom.css"), tags$script(HTML(COLOR_INPUT_JS))),

  nav_panel(
    "Start", icon = bsicons::bs_icon("house"),
    div(class = "container-lg py-3", style = "max-width: 1000px;",
      h1(class = "fw-bold mb-1", "Turn raw assay data into figures you can publish"),
      p(class = "lead text-muted",
        "Bring a spreadsheet. Leave with vector figures, statistics, and the R code behind them."),
      layout_columns(
        col_widths = c(4, 4, 4), class = "my-4",
        card(class = "h-100",
          card_body(
            h5(class = "fw-bold", bsicons::bs_icon("1-circle-fill"), " Load your data"),
            p(class = "small mb-0",
              "Upload a CSV, XLSX, or paste straight out of Excel. ",
              "You never rename a column \u2014 the app shows you its guesses and ",
              "you correct any it got wrong. Every tab has a demo file to copy over."))),
        card(class = "h-100",
          card_body(
            h5(class = "fw-bold", bsicons::bs_icon("2-circle-fill"), " Shape the figure"),
            p(class = "small mb-0",
              "Choose groups and their order, pick colours, set the comparisons ",
              "and confidence level. The figure on screen is exactly what gets ",
              "exported \u2014 no surprises at download time."))),
        card(class = "h-100",
          card_body(
            h5(class = "fw-bold", bsicons::bs_icon("3-circle-fill"), " Take it away"),
            p(class = "small mb-0",
              "PDF and SVG stay editable in Illustrator. Or take the ZIP: your ",
              "data, the analysis code, and a script that rebuilds the figures ",
              "on your own machine.")))
      ),
      h4(class = "fw-bold mt-4 mb-2", "Pick an assay"),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(class = "h-100", card_body(
          h5(class = "fw-bold mb-1", "Limiting dilution"),
          p(class = "small text-muted",
            "Cells seeded and organoids counted, or wells scored +/- per dose. ",
            "Frequency bar chart with confidence intervals and significance ",
            "brackets, plus the dose\u2013response curve."),
          actionButton("go_lda", "Open", class = "btn-sm btn-primary"))),
        card(class = "h-100", card_body(
          h5(class = "fw-bold mb-1", "RT-qPCR"),
          p(class = "small text-muted",
            "A long-format Ct table becomes a \u0394\u0394Ct bar chart faceted by ",
            "gene, with tests run on \u0394Ct where they belong."),
          actionButton("go_qpcr", "Open", class = "btn-sm btn-primary"))),
        card(class = "h-100", card_body(
          h5(class = "fw-bold mb-1", "Volcano"),
          p(class = "small text-muted",
            "Any DESeq2, edgeR, or limma results table. Thresholds, gene labels, ",
            "and a hits table you can download."),
          actionButton("go_volc", "Open", class = "btn-sm btn-primary")))
      ),
      div(class = "alert alert-success border mt-4",
        h6(class = "fw-bold mb-1",
           bsicons::bs_icon("shield-lock-fill"), " Your data never leaves this computer"),
        p(class = "small mb-0",
          "This page runs R itself, inside your browser, compiled to WebAssembly. ",
          "There is no server on the other end to receive anything. Files you open ",
          "are read into this tab's sandboxed memory, and everything is discarded ",
          "when you close it. See ", strong("How it works"),
          " for how to check that claim yourself.")),
      div(class = "alert alert-light border mt-3 small",
        bsicons::bs_icon("collection"), " ",
        strong("Got a sheet with many experiments in it?"),
        " On the limiting-dilution tab, set ", strong("Split figures by"),
        " to a column such as cell line or screen. You get one figure per level, ",
        "all in a single ZIP \u2014 which is usually the actual request behind ",
        "\u201ccan you make me a figure?\u201d")
    )
  ),
  nav_panel(
    "Limiting dilution", icon = bsicons::bs_icon("circle-half"),
    layout_sidebar(sidebar = sidebar(width = 350, open = "open", ldaSidebarUI("lda")),
                   ldaUI("lda"))
  ),
  nav_panel(
    "RT-qPCR", icon = bsicons::bs_icon("thermometer-half"),
    layout_sidebar(sidebar = sidebar(width = 350, open = "open", qpcrSidebarUI("qpcr")),
                   qpcrUI("qpcr"))
  ),
  nav_panel(
    "Volcano", icon = bsicons::bs_icon("triangle"),
    layout_sidebar(sidebar = sidebar(width = 350, open = "open", volcanoSidebarUI("volc")),
                   volcanoUI("volc"))
  ),
  nav_spacer(),
  nav_panel(
    "How it works", icon = bsicons::bs_icon("info-circle"),
    card(card_body(markdown(
"## Getting your data in

Every tab accepts the same three things: the **example dataset**, an
**uploaded** `.csv` / `.tsv` / `.xlsx`, or data **pasted straight out of
Excel**. You never have to rename a column — after you upload, the app shows
a *Match your columns* panel where you point your headers at the fields it
needs. It guesses correctly most of the time.

## Limiting dilution

Two assay layouts, one model.

| Layout | You record | Model |
|---|---|---|
| **Counted organoids** | cells seeded, organoids formed | `cloglog(p) = β` — every cell is a trial |
| **Scored wells** | cells/well, wells tested, wells positive | `cloglog(p) = log(dose) + β` — classic ELDA |

Both are a binomial GLM with a complementary log-log link, so the
organoid-forming frequency is `f = exp(β)`, the confidence interval is Wald
on the log scale, and comparing two groups is a Wald test on `β_A − β_B`.
The second layout is exactly the ELDA method of Hu & Smyth (2009), written
as an offset instead of a separate function.

The **goodness-of-fit** line under the comparison table is the residual
deviance test. A small p-value means one frequency per group does *not*
explain the data — usually non-linear dose scaling, or replicates that
disagree with each other. Worth looking at before you trust the number.

## RT-qPCR

Technical replicates are averaged first. ΔCt is taken against the
housekeeping gene *within the same biological replicate*, ΔΔCt against the
mean ΔCt of the reference group, and the bar height is `2^-ΔΔCt`.

Statistics run on **ΔCt values, not fold changes** — fold changes are
log-normal, so a t-test on them is the wrong test. The p-values in the table
are Welch t-tests on ΔCt, Benjamini–Hochberg adjusted across all
gene × group comparisons.

## Volcano

Point any differential-expression table at it — DESeq2, edgeR, limma, or a
spreadsheet a collaborator emailed you. Cap the y-axis if a couple of genes
with astronomically small p-values are flattening everything else; capped
points are drawn as triangles so nobody is misled.

## Where your data goes

Nowhere. That is the whole design, not a policy promise.

This app is compiled to WebAssembly and served as plain static files. When you
open the page your browser downloads R itself — the interpreter, the graphics
devices, ggplot2, every package — and runs the entire analysis locally. There
is no R process on a server, because there is no server beyond the one that
handed you the files.

So when you upload a CSV, the browser reads it into this tab's sandboxed
memory. The bytes never travel anywhere. The same is true of anything you
paste. Figures and ZIPs are generated in the tab and handed straight to your
downloads folder. Close the tab and it is all gone — nothing is written to
disk, no cookie, no account, no log of what you analysed.

**Check it yourself**, because you should not take this on faith:

- Open DevTools (F12) → **Network**, then load a file and make a figure. You
  will see requests only while the page first loads, all to this site's own
  address, and none once R has started.
- Or simpler: load the page, turn off your wifi, and keep working. Everything
  still works, because there is nothing left to talk to.

Two things this does *not* protect against, so that the boundary is clear.
The **app** is public — anyone with the link can open it, even when the source
repository is private, unless your organisation is on GitHub Enterprise Cloud
with private Pages enabled. And anything committed to that repository is
public too, which is why real data must never be committed to it. The link
being hard to guess is not a security measure.

## Taking it with you

**A note on PDF.** On machines without a working cairo library (macOS without
XQuartz, for instance) the PDF device falls back to R's built-in one, which
cannot draw characters outside Latin-1 — a Δ or a • in a title comes out as a
dot. SVG and PNG are unaffected, so use SVG if you see that. Servers running
Linux, including shinyapps.io, have cairo and are fine.

**PDF** and **SVG** are vector — every point, line and letter stays a
separate editable object in Illustrator, Inkscape or Affinity. Use those for
papers; PNG is for slides and Slack.

The **Download data + R code** button on the limiting-dilution tab gives you
a ZIP with your data as `data.csv`, the fitting and plotting functions in
`R/`, and an `analysis.R` that rebuilds the exact figures on screen. Open it
in RStudio, hit Source, and change whatever you like. Nothing is locked
inside the app.
    ")))
  )
)

server <- function(input, output, session) {
  observeEvent(input$go_lda,  nav_select("nav", "Limiting dilution"))
  observeEvent(input$go_qpcr, nav_select("nav", "RT-qPCR"))
  observeEvent(input$go_volc, nav_select("nav", "Volcano"))

  ldaServer("lda", helper_dir = HELPER_DIR)
  qpcrServer("qpcr")
  volcanoServer("volc")
}

shinyApp(ui, server)
