# ── Reading a protocol sheet ──────────────────────────────────────────────
# The sheets people actually keep are not tidy tables. They are blocks of
# numbers separated by blank rows, with the cell line and media written as a
# floating title above each block, no column headers, and — crucially — two
# different cell counts side by side:
#
#     Cdx2-ApcKO   R-media                  <- floating metadata, not headers
#     cells seeded                          <- a header for one column only
#     100    87    9    DMSO                <- nominal | actual | organoids | group
#     200   147   13    DMSO
#
# 100 is what the protocol said to seed. 87 is what was counted into the well.
# Fitting the model to 100 instead of 87 biases the frequency by ~20%, so the
# importer has to tell them apart rather than take the first number it sees.

# Parse a single cell to a number, tolerating "1,200", " 87 ", "300 cells".
as_num <- function(x) {
  if (is.na(x)) return(NA_real_)
  v <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", as.character(x))))
  if (length(v) != 1 || !is.finite(v)) NA_real_ else v
}

# Words that mean "this is a column header", not "this is my experiment".
HEADER_WORDS <- paste0(
  "cell|seed|plate|organoid|colony|colonies|sphere|count|number|no\\.|#|",
  "dose|well|group|condition|treatment|sample|rep|total|formed|input")

is_header_text <- function(x) {
  x <- tolower(trimws(as.character(x)))
  grepl(HEADER_WORDS, x)
}

# How "round" a set of numbers is — protocol doses are 100/200/300, counted
# cells are 87/147/245.
roundness <- function(v) {
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) return(0)
  max(mean(v %% 100 == 0), mean(v %% 50 == 0), mean(v %% 25 == 0),
      mean(v %% 10 == 0))
}

# ---- split a raw grid into blocks on blank rows ---------------------------
split_blocks <- function(grid) {
  blank <- apply(grid, 1, function(r) all(is.na(r) | trimws(as.character(r)) == ""))
  run <- cumsum(blank)
  idx <- split(which(!blank), run[!blank])
  Filter(function(i) length(i) >= 2, idx)     # a block needs >= 2 rows
}

# ---- classify one block ---------------------------------------------------
parse_block <- function(grid, rows) {
  b <- grid[rows, , drop = FALSE]
  nums <- suppressWarnings(
    matrix(vapply(as.matrix(b), as_num, numeric(1)), nrow = nrow(b)))

  n_num <- rowSums(is.finite(nums))
  is_data <- n_num >= 2                        # at least cells + organoids

  # Everything above the data that is pure text is the block's own label.
  lab_rows <- which(!is_data)
  label <- character()
  for (i in lab_rows) {
    vals <- trimws(as.character(unlist(b[i, ])))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    vals <- vals[!is_header_text(vals)]
    label <- c(label, vals)
  }
  label <- paste(unique(label), collapse = " ")

  d <- b[is_data, , drop = FALSE]
  nm <- nums[is_data, , drop = FALSE]
  if (!nrow(d)) return(NULL)

  # Numeric columns = finite in every data row. Text columns = the rest.
  num_col <- which(colSums(is.finite(nm)) == nrow(nm))
  txt_col <- setdiff(
    which(colSums(!is.na(d) & trimws(as.matrix(d)) != "") == nrow(d)), num_col)
  if (length(num_col) < 2) return(NULL)

  # The organoid column is the one that never exceeds any other numeric
  # column — you cannot form more organoids than you seeded cells.
  tot <- colSums(nm[, num_col, drop = FALSE])
  org <- num_col[which.min(tot)]
  cell_cols <- setdiff(num_col, org)

  # Of the remaining count columns, the protocol dose is the round, repeating
  # one; the real count is the noisy, near-unique one.
  score <- vapply(cell_cols, function(j) {
    v <- nm[, j]
    roundness(v) - length(unique(v)) / length(v)
  }, numeric(1))
  if (length(cell_cols) >= 2) {
    nominal <- cell_cols[which.max(score)]
    actual  <- cell_cols[which.min(score)]
  } else {
    nominal <- NA_integer_
    actual  <- cell_cols[1]
  }

  grp <- if (length(txt_col)) txt_col[which.max(vapply(txt_col, function(j)
    length(unique(trimws(as.character(d[[j]])))), integer(1)))] else NA_integer_

  out <- data.frame(
    CellCount = nm[, actual],
    Organoids = nm[, org],
    Group     = if (is.na(grp)) "ungrouped" else trimws(as.character(d[[grp]])),
    stringsAsFactors = FALSE)
  out$NominalCells <- if (is.na(nominal)) NA_real_ else nm[, nominal]
  out$Batch <- if (nzchar(label)) label else "block"

  attr(out, "cols") <- list(actual = actual, nominal = nominal,
                            organoids = org, group = grp)
  out
}

# ---- the entry point ------------------------------------------------------
# Returns a tidy data.frame plus an attached report of every decision made,
# so the app can show its working instead of asking to be trusted.
parse_protocol_grid <- function(grid) {
  grid <- grid[, colSums(!is.na(grid) & trimws(as.matrix(grid)) != "") > 0,
               drop = FALSE]
  blocks <- split_blocks(grid)
  if (!length(blocks)) stop("No blocks of data found in that sheet.")

  parsed <- lapply(blocks, function(r) parse_block(grid, r))
  keep <- !vapply(parsed, is.null, logical(1))
  parsed <- parsed[keep]
  if (!length(parsed)) stop("Found blocks, but none had two numeric columns.")

  report <- lapply(seq_along(parsed), function(i) {
    p <- parsed[[i]]; cc <- attr(p, "cols")
    data.frame(
      Block = i, Label = p$Batch[1], Rows = nrow(p),
      Groups = paste(unique(p$Group), collapse = ", "),
      `Cells column` = paste0("col ", cc$actual),
      `Protocol column` = if (is.na(cc$nominal)) "none found"
                          else paste0("col ", cc$nominal),
      `Organoids column` = paste0("col ", cc$organoids),
      check.names = FALSE, stringsAsFactors = FALSE)
  })

  # Blocks often share a title (same line, same media) — make them distinct
  # so they can be used to split figures.
  out <- do.call(rbind, lapply(seq_along(parsed), function(i) {
    p <- parsed[[i]]; attr(p, "cols") <- NULL; p
  }))
  labs <- vapply(parsed, function(p) p$Batch[1], character(1))
  if (anyDuplicated(labs)) {
    for (i in seq_along(parsed))
      parsed[[i]]$Batch <- sprintf("%s #%d", labs[i], i)
    out <- do.call(rbind, lapply(parsed, function(p) {
      attr(p, "cols") <- NULL; p }))
  }
  rownames(out) <- NULL

  attr(out, "report") <- do.call(rbind, report)
  attr(out, "has_nominal") <- any(is.finite(out$NominalCells))
  out
}

parse_protocol_file <- function(path, name = path, sheet = 1) {
  ext <- tolower(tools::file_ext(name))
  grid <- if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Excel files need the 'readxl' package.")
    as.data.frame(readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                                     col_types = "text", .name_repair = "minimal"),
                  stringsAsFactors = FALSE)
  } else {
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    utils::read.delim(path, sep = sep, header = FALSE, colClasses = "character",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      blank.lines.skip = FALSE)
  }
  parse_protocol_grid(grid)
}

protocol_sheets <- function(path, name = path) {
  if (tolower(tools::file_ext(name)) %in% c("xlsx", "xls") &&
      requireNamespace("readxl", quietly = TRUE))
    readxl::excel_sheets(path) else "1"
}

# How much would using the protocol numbers change the answer? Shown as a
# warning, because in practice it is never negligible.
nominal_cost <- function(df, groups, mode = "cells") {
  if (!any(is.finite(df$NominalCells))) return(NULL)
  d <- df[is.finite(df$NominalCells), , drop = FALSE]
  g <- intersect(groups, unique(d$Group))
  if (!length(g)) return(NULL)
  a <- tryCatch(lda_fit(d, g, mode)$est, error = function(e) NULL)
  d2 <- d; d2$CellCount <- d2$NominalCells
  b <- tryCatch(lda_fit(d2, g, mode)$est, error = function(e) NULL)
  if (is.null(a) || is.null(b)) return(NULL)
  data.frame(Group = a$Group,
             `Using counted cells` = signif(a$Frequency, 3),
             `Using protocol numbers` = signif(b$Frequency, 3),
             Difference = sprintf("%+.0f%%", 100 * (b$Frequency / a$Frequency - 1)),
             check.names = FALSE)
}
