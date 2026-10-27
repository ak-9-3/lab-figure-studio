# ── Getting data in ───────────────────────────────────────────────────────
# One module handles data entry for every assay. The design rule here is that
# the person uploading should not have to know, or care, what shape their
# sheet is in. They press Upload. We work it out, say what we worked out, and
# let them correct us.

read_any <- function(path, name = path, sheet = 1) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) {
    validate(need(requireNamespace("readxl", quietly = TRUE),
                  "Excel files need the 'readxl' package. Save as .csv instead."))
    as.data.frame(readxl::read_excel(path, sheet = sheet))
  } else if (ext %in% c("tsv", "txt")) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

# Excel copy-paste is tab-separated; a pasted CSV should work too.
read_pasted <- function(txt) {
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  sep <- if (grepl("\t", txt)) "\t" else if (grepl(";", txt)) ";" else ","
  utils::read.delim(text = txt, sep = sep, stringsAsFactors = FALSE,
                    check.names = FALSE)
}

# A column name means nothing on a sheet with no header row, so every picker
# shows what is actually in the column instead.
col_preview <- function(v, n = 4) {
  v <- v[!is.na(v) & trimws(as.character(v)) != ""]
  if (!length(v)) return("empty")
  txt <- paste(utils::head(trimws(as.character(v)), n), collapse = ", ")
  if (length(v) > n) txt <- paste0(txt, ", …")
  if (nchar(txt) > 42) txt <- paste0(substr(txt, 1, 40), "…")
  txt
}

looks_like_protocol <- function(df) {
  auto  <- sum(grepl("^\\.\\.\\.[0-9]+$|^X[0-9]+$", names(df)))
  blank <- sum(apply(df, 1, function(r)
    all(is.na(r) | trimws(as.character(r)) == "")))
  auto >= 2 || blank >= 1
}

# Fuzzy header matching: strip case, spaces, underscores, dots.
norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Aliases are in priority order: for cell counts, "cells counted" beats
# "cells seeded", because the count is what actually went into the well.
guess_col <- function(cols, aliases, exclude = character()) {
  n <- norm_name(cols); a <- norm_name(aliases)
  # An empty string matches everything under grepl(), so a sheet with blank
  # or auto-generated headers would otherwise "match" every field.
  ok <- nzchar(n) & !grepl("^[0-9]+$", n) & !cols %in% exclude
  a  <- unique(a[nzchar(a)])
  if (!any(ok) || !length(a)) return(NA_character_)
  for (al in a) {
    hit <- which(ok & n == al)
    if (length(hit)) return(cols[hit[1]])
  }
  # Partial matches only for aliases long enough to mean something — "n"
  # would otherwise match "condition".
  for (al in a[nchar(a) >= 3]) {
    hit <- which(ok & (grepl(al, n, fixed = TRUE) |
                       vapply(n, grepl, logical(1), x = al, fixed = TRUE)))
    if (length(hit)) return(cols[hit[1]])
  }
  NA_character_
}

# Guess every field at once, never giving one column to two fields.
guess_all <- function(cols, spec) {
  taken <- character()
  out <- list()
  for (k in names(spec)) {
    g <- guess_col(cols, c(k, spec[[k]]$aliases), exclude = taken)
    out[[k]] <- if (is.na(g)) "" else g
    if (!is.na(g)) taken <- c(taken, g)
  }
  out
}

required_fields <- function(spec)
  names(spec)[!vapply(spec, function(s) isTRUE(s$optional), logical(1))]

# ---- the module -----------------------------------------------------------
dataInputUI <- function(id, spec, example_label = "Example dataset",
                        protocol = FALSE) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("src"), NULL, inline = TRUE,
                 choiceNames  = c(example_label, "Upload", "Paste"),
                 choiceValues = c("example", "file", "paste")),
    conditionalPanel(
      sprintf("input['%s'] == 'file'", ns("src")),
      fileInput(ns("file"), NULL, buttonLabel = "Choose file…",
                placeholder = "csv / tsv / xlsx",
                accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
      uiOutput(ns("sheet_ui"))
    ),
    conditionalPanel(
      sprintf("input['%s'] == 'paste'", ns("src")),
      textAreaInput(ns("paste"), NULL, rows = 7,
                    placeholder = "Select the cells in Excel, copy, paste here."),
      div(class = "form-text mb-2",
          "Headers or no headers, any column order: the app works out which",
          " column is which from the numbers and tells you what it decided.")
    ),
    downloadLink(ns("template"), "Download a blank template / demo CSV"),
    uiOutput(ns("detected")),
    uiOutput(ns("preview")),
    uiOutput(ns("mapping")),
    uiOutput(ns("sanity")),
    uiOutput(ns("status"))
  )
}

dataInputServer <- function(id, spec, example, protocol = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    req_fields <- required_fields(spec)

    output$sheet_ui <- renderUI({
      req(input$file)
      if (!tolower(tools::file_ext(input$file$name)) %in% c("xlsx", "xls"))
        return(NULL)
      sh <- tryCatch(readxl::excel_sheets(input$file$datapath),
                     error = function(e) NULL)
      if (length(sh) < 2) return(NULL)
      selectInput(ns("sheet"), "Sheet", choices = sh, selected = sh[1])
    })

    # ---- read the file both ways, then decide ---------------------------
    # Pasted text goes down exactly the same road as an uploaded file: it is
    # saved to a temporary file and read both ways. So whatever you paste —
    # with headers, without, in blocks — gets the same "work it out" logic.
    pasted <- debounce(reactive(input$paste), 600)
    source_file <- reactive({
      if (identical(input$src, "file")) {
        req(input$file)
        return(list(path = input$file$datapath, name = input$file$name))
      }
      txt <- pasted(); req(txt, nzchar(trimws(txt)))
      ext <- if (grepl("\t", txt)) "tsv" else "csv"
      if (!grepl("\t", txt) && grepl(";", txt)) { txt <- gsub(";", "\t", txt); ext <- "tsv" }
      f <- tempfile(fileext = paste0(".", ext))
      writeLines(txt, f)
      list(path = f, name = basename(f))
    })

    attempts <- reactive({
      src <- source_file()
      sh <- if (identical(input$src, "file")) input$sheet %||% 1 else 1
      list(
        plain = tryCatch(read_any(src$path, src$name, sh),
                         error = function(e) NULL),
        prot  = if (protocol)
          tryCatch(parse_protocol_file(src$path, src$name, sh),
                   error = function(e) NULL) else NULL)
    })

    # How many required fields can we name from the headers alone?
    n_guessed <- function(df) {
      if (is.null(df) || !ncol(df)) return(0L)
      g <- guess_all(names(df), spec)
      sum(vapply(req_fields, function(k) nzchar(g[[k]]), logical(1)))
    }

    auto_mode <- reactive({
      a <- attempts()
      if (!is.null(a$prot) && n_guessed(a$plain) < length(req_fields))
        "protocol" else "plain"
    })
    # "Use my header row" / "Work it out" is a choice about THIS data. New
    # data, new decision: forget the override whenever the input changes.
    override <- reactiveVal(NULL)
    observeEvent(source_file(), override(NULL))
    observeEvent(input$read_as, {
      a <- tryCatch(auto_mode(), error = function(e) NULL)
      override(if (identical(input$read_as, a)) NULL else input$read_as)
    }, ignoreInit = TRUE)
    read_mode <- reactive({
      if (identical(input$src, "example")) return("plain")
      override() %||% auto_mode()
    })

    # ---- tell them what we did, and let them overrule it ----------------
    output$detected <- renderUI({
      if (identical(input$src, "example")) return(NULL)
      a <- tryCatch(attempts(), error = function(e) NULL)
      if (is.null(a)) return(NULL)
      m <- read_mode()
      both <- !is.null(a$prot) && !is.null(a$plain)
      tagList(
        hr(),
        if (m == "protocol") div(
          class = "alert alert-success py-2 px-2 small mb-2",
          bsicons::bs_icon("magic"), strong(" Read as a lab sheet."),
          " No usable header row, so the columns were worked out from the",
          " numbers themselves — no naming needed.")
        else div(
          class = "alert alert-info py-2 px-2 small mb-2",
          bsicons::bs_icon("table"), strong(" Read as a plain table."),
          " The first row was used as the header."),
        if (both) radioButtons(
          ns("read_as"), NULL, inline = TRUE, selected = m,
          choiceNames  = c("Work it out from the numbers", "Use my header row"),
          choiceValues = c("protocol", "plain"))
      )
    })

    # ---- always show them their own data --------------------------------
    output$preview <- renderUI({
      d <- tryCatch(raw(), error = function(e) NULL)
      if (is.null(d) || !nrow(d)) return(NULL)
      div(class = "mb-2",
        div(class = "fw-bold small mb-1",
            sprintf("Your data (%d rows, showing first %d)", nrow(d), min(5, nrow(d)))),
        div(style = "overflow-x:auto; font-size:.75rem;",
            tableOutput(ns("preview_tbl"))))
    })
    output$preview_tbl <- renderTable({
      d <- tryCatch(raw(), error = function(e) NULL)
      req(d, nrow(d) > 0)
      show <- utils::head(as.data.frame(d), 5)
      show[] <- lapply(show, function(x) {
        x <- as.character(x); x[is.na(x)] <- ""; substr(x, 1, 16) })
      show
    }, striped = TRUE, spacing = "xs", width = "100%")

    # ---- raw table, exactly as supplied ---------------------------------
    raw <- reactive({
      switch(input$src,
        example = example,
        { a <- attempts()
          if (read_mode() == "protocol") a$prot else a$plain })
    })

    # Output that is already in our own shape needs no mapping step.
    canonical <- reactive(
      identical(input$src, "example") || identical(read_mode(), "protocol"))

    # ---- the mapping step, only when it is actually needed --------------
    output$mapping <- renderUI({
      df <- tryCatch(raw(), error = function(e) NULL)
      if (is.null(df) || !ncol(df) || canonical()) return(NULL)
      cols <- names(df)
      gs <- guess_all(cols, spec)
      tagList(
        div(class = "fw-bold mb-1", "Which column is which?"),
        div(class = "form-text mb-2",
            "Pick by the values, not the name — each option shows what is",
            " actually in that column."),
        lapply(names(spec), function(k) {
          s <- spec[[k]]
          guess <- gs[[k]]
          blank <- stats::setNames(
            "", if (isTRUE(s$optional)) "— none —" else "— choose —")
          shown <- stats::setNames(cols, vapply(cols, function(cn)
            sprintf("%s  —  %s", cn, col_preview(df[[cn]])), character(1)))
          selectInput(ns(paste0("map_", k)),
                      paste0(s$label, if (isTRUE(s$optional)) " (optional)" else ""),
                      choices  = c(blank, shown),
                      selected = guess)
        })
      )
    })

    # ---- catch an impossible mapping before it reaches the model --------
    # Organoids cannot exceed cells. If they do, offer to swap rather than
    # letting a model error surface three panels away.
    swap_pair <- reactive({
      if (canonical()) return(NULL)
      pairs <- list(c("Organoids", "CellCount"), c("Positive", "Tested"))
      for (pr in pairs) if (all(pr %in% names(spec))) return(pr)
      NULL
    })

    output$sanity <- renderUI({
      pr <- swap_pair(); if (is.null(pr)) return(NULL)
      df <- tryCatch(raw(), error = function(e) NULL); if (is.null(df)) return(NULL)
      pick <- function(k) {
        v <- input[[paste0("map_", k)]]
        if (!is.null(v)) return(v)
        guess_all(names(df), spec)[[k]]
      }
      a <- pick(pr[1]); b <- pick(pr[2])
      if (!nzchar(a) || !nzchar(b) || !a %in% names(df) || !b %in% names(df))
        return(NULL)
      ev <- suppressWarnings(as.numeric(df[[a]]))
      tr <- suppressWarnings(as.numeric(df[[b]]))
      bad <- sum(is.finite(ev) & is.finite(tr) & ev > tr)
      if (!bad) return(NULL)
      div(class = "alert alert-danger py-2 px-2 small",
          bsicons::bs_icon("exclamation-octagon-fill"),
          strong(sprintf(" '%s' is bigger than '%s' in %d rows.",
                         spec[[pr[1]]]$label, spec[[pr[2]]]$label, bad)),
          " That cannot happen, so these two are probably the wrong way round.",
          div(class = "mt-2",
              actionButton(ns("do_swap"), "Swap them", class = "btn-sm btn-danger")))
    })
    outputOptions(output, "sanity", suspendWhenHidden = FALSE)

    observeEvent(input$do_swap, {
      pr <- swap_pair(); req(pr)
      a <- input[[paste0("map_", pr[1])]]; b <- input[[paste0("map_", pr[2])]]
      updateSelectInput(session, paste0("map_", pr[1]), selected = b)
      updateSelectInput(session, paste0("map_", pr[2]), selected = a)
    })

    output$template <- downloadHandler(
      function() paste0(id, "_template.csv"),
      function(f) utils::write.csv(example, f, row.names = FALSE))

    # ---- canonical, type-coerced, validated -----------------------------
    tidy <- reactive({
      df <- raw()
      validate(need(!is.null(df) && nrow(df) > 0,
                    "No data yet — upload, paste, or use the example."))
      if (canonical()) {
        out <- df
      } else {
        out <- data.frame(row.names = seq_len(nrow(df)))
        for (k in names(spec)) {
          sel <- input[[paste0("map_", k)]]
          # Before the picker has rendered its input, use the guess directly,
          # so a well-formed file needs no round-trip and no empty flash.
          if (is.null(sel)) sel <- guess_all(names(df), spec)[[k]]
          if (!nzchar(sel)) {
            validate(need(isTRUE(spec[[k]]$optional),
                          sprintf("Still need to know which column holds '%s'.",
                                  spec[[k]]$label)))
            next
          }
          out[[k]] <- df[[sel]]
        }
      }
      for (k in intersect(names(spec), names(out))) {
        if (identical(spec[[k]]$type, "numeric")) {
          v <- suppressWarnings(as.numeric(gsub("[, ]", "", as.character(out[[k]]))))
          validate(need(!all(is.na(v)),
                        sprintf("'%s' has no numbers in it — wrong column?",
                                spec[[k]]$label)))
          out[[k]] <- v
        } else {
          out[[k]] <- trimws(as.character(out[[k]]))
        }
      }
      out <- out[stats::complete.cases(
        out[, intersect(names(spec), names(out)), drop = FALSE]), , drop = FALSE]
      validate(need(nrow(out) > 0, "Every row had a missing value after mapping."))
      out
    })

    output$status <- renderUI({
      df <- tryCatch(tidy(), error = function(e) NULL)
      if (is.null(df)) return(NULL)
      div(class = "alert alert-success py-1 px-2 mt-2 mb-0 small d-flex justify-content-between align-items-center",
          span(sprintf("%d rows ready.", nrow(df))),
          downloadLink(ns("dl_data"), tagList(bsicons::bs_icon("download"), " this data as CSV")))
    })
    outputOptions(output, "status", suspendWhenHidden = FALSE)

    output$dl_data <- downloadHandler(
      function() paste0(id, "_data.csv"),
      function(f) utils::write.csv(tidy(), f, row.names = FALSE))

    tidy
  })
}
