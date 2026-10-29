# =============================================================================
# tools/repro/compare.R
#
# Comparison of a reproduced output with its reference file.
#   table : formatted manuscript tables, cells like "5.19 (1.42–16.65)".
#           Every cell is compared as text (exact match at the reported
#           precision); numbers are also parsed out of each cell so the
#           report can give absolute/relative differences for mismatches.
#   grid  : numeric spreadsheets (Supplementary Data 2). Numeric columns are
#           compared with a relative tolerance (default 1e-9); text columns
#           exactly. Rows are matched on the key columns, not row order.
#   png   : figures. Byte-identical files PASS. Otherwise both images are
#           decoded and compared pixel by pixel; any pixel difference is
#           reported as FAIL with the fraction of differing pixels, so that a
#           human can inspect it (rendering is not expected to differ when
#           the same R/package versions are used).
#
# Each comparison returns a one-row data.frame:
#   id, reproduced, reference, type, status (PASS/FAIL/BLOCKED), n_compared,
#   n_mismatch, max_abs_diff, max_rel_diff, tolerance, detail
# =============================================================================

read_xlsx_copy <- function(path, ...) {
  # Read via a temporary copy: OneDrive/Excel can hold locks on the original.
  tmp <- tempfile(fileext = ".xlsx")
  if (!file.copy(path, tmp)) stop("Cannot copy ", path)
  on.exit(unlink(tmp))
  suppressMessages(readxl::read_excel(tmp, ...))
}

parse_numbers <- function(x) {
  # all numbers in a formatted cell, e.g. "5.19 (1.42–16.65)" -> 5.19, 1.42, 16.65
  # A '-' directly after a digit is a range separator ("0.00-0.01"), not a sign.
  m <- regmatches(x, gregexpr("(?<![0-9])-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?", x, perl = TRUE))
  lapply(m, as.numeric)
}

result_row <- function(id, rep, ref, type, status, n = NA, nm = NA, mad = NA, mrd = NA, tol = NA, detail = "") {
  data.frame(id = id, reproduced = rep, reference = ref, type = type, status = status,
             n_compared = n, n_mismatch = nm, max_abs_diff = mad, max_rel_diff = mrd,
             tolerance = tol, detail = detail, stringsAsFactors = FALSE)
}

compare_table <- function(id, rep, ref) {
  a <- as.data.frame(read_xlsx_copy(rep, col_types = "text"))
  b <- as.data.frame(read_xlsx_copy(ref, col_types = "text"))
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b)))
    return(result_row(id, rep, ref, "table", "FAIL",
                      detail = sprintf("shape/header differ: reproduced %dx%d, reference %dx%d",
                                       nrow(a), ncol(a), nrow(b), ncol(b))))
  va <- as.matrix(a); vb <- as.matrix(b)
  va[is.na(va)] <- ""; vb[is.na(vb)] <- ""
  mism <- which(trimws(va) != trimws(vb), arr.ind = TRUE)
  diffs <- numeric(0); rels <- numeric(0)
  if (nrow(mism)) {
    for (k in seq_len(nrow(mism))) {
      na <- parse_numbers(va[mism[k, , drop = FALSE]])[[1]]
      nb <- parse_numbers(vb[mism[k, , drop = FALSE]])[[1]]
      if (length(na) == length(nb) && length(na)) {
        diffs <- c(diffs, abs(na - nb)); rels <- c(rels, abs(na - nb) / pmax(abs(nb), 1e-12))
      }
    }
  }
  detail <- if (nrow(mism)) paste(utils::head(sprintf("[%d,%s] '%s' vs '%s'", mism[, 1], colnames(va)[mism[, 2]],
                                                      va[mism], vb[mism]), 5), collapse = "; ") else ""
  result_row(id, rep, ref, "table", if (nrow(mism)) "FAIL" else "PASS",
             n = length(va), nm = nrow(mism),
             mad = if (length(diffs)) max(diffs) else 0, mrd = if (length(rels)) max(rels) else 0,
             tol = "exact text at reported precision", detail = detail)
}

compare_grid <- function(id, rep, ref, rel_tol = 1e-9,
                         keys = c("Setting", "Age group", "Duration (days)", "Entry day (day of year)")) {
  # guess_max over all rows: columns that start with NA (e.g. Death BRR for
  # 18-64, where vaccine death risk is 0) are otherwise guessed as logical and
  # their numbers silently coerced to TRUE/FALSE.
  a <- as.data.frame(read_xlsx_copy(rep, sheet = 1, guess_max = 1e6))
  b <- as.data.frame(read_xlsx_copy(ref, sheet = 1, guess_max = 1e6))
  if (!setequal(names(a), names(b)))
    return(result_row(id, rep, ref, "grid", "FAIL", detail = "column names differ"))
  if (nrow(a) != nrow(b))
    return(result_row(id, rep, ref, "grid", "FAIL", detail = sprintf("row count %d vs %d", nrow(a), nrow(b))))
  ka <- do.call(paste, c(a[keys], sep = "|")); kb <- do.call(paste, c(b[keys], sep = "|"))
  if (anyDuplicated(ka) || !setequal(ka, kb))
    return(result_row(id, rep, ref, "grid", "FAIL", detail = "key columns do not identify the same rows"))
  b <- b[match(ka, kb), names(a)]
  n <- 0; nm <- 0; mad <- 0; mrd <- 0; worst <- ""
  for (col in names(a)) {
    x <- a[[col]]; y <- b[[col]]
    if (is.numeric(x) && is.numeric(y)) {
      both_na <- is.na(x) & is.na(y)
      one_na <- xor(is.na(x), is.na(y))
      d <- abs(x - y); r <- d / pmax(abs(y), 1e-300)
      bad <- one_na | (!both_na & !one_na & d > rel_tol * pmax(abs(y), 1e-12))
      n <- n + length(x); nm <- nm + sum(bad, na.rm = TRUE)
      if (any(!both_na & !one_na)) {
        mad <- max(mad, d[!both_na & !one_na]); mrd <- max(mrd, r[!both_na & !one_na & y != 0], 0)
      }
      if (any(bad, na.rm = TRUE) && !nzchar(worst)) worst <- col
    } else {
      bad <- !(is.na(x) & is.na(y)) & (is.na(x) != is.na(y) | as.character(x) != as.character(y))
      n <- n + length(x); nm <- nm + sum(bad, na.rm = TRUE)
      if (any(bad, na.rm = TRUE) && !nzchar(worst)) worst <- col
    }
  }
  result_row(id, rep, ref, "grid", if (nm == 0) "PASS" else "FAIL", n = n, nm = nm, mad = mad, mrd = mrd,
             tol = sprintf("relative %.0e (numeric); exact (text)", rel_tol),
             detail = if (nzchar(worst)) paste("first mismatching column:", worst) else "")
}

compare_png <- function(id, rep, ref) {
  if (unname(tools::md5sum(rep)) == unname(tools::md5sum(ref)))
    return(result_row(id, rep, ref, "png", "PASS", n = 1, nm = 0, mad = 0, mrd = 0,
                      tol = "byte-identical", detail = "byte-identical"))
  if (!requireNamespace("png", quietly = TRUE))
    return(result_row(id, rep, ref, "png", "FAIL", detail = "files differ; 'png' package unavailable for pixel comparison"))
  a <- png::readPNG(rep); b <- png::readPNG(ref)
  if (!identical(dim(a), dim(b)))
    return(result_row(id, rep, ref, "png", "FAIL",
                      detail = sprintf("image size differs: %s vs %s", paste(dim(a), collapse = "x"), paste(dim(b), collapse = "x"))))
  d <- abs(a - b)
  px <- if (length(dim(d)) == 3) apply(d, c(1, 2), max) else d
  frac <- mean(px > 0)
  result_row(id, rep, ref, "png", if (frac == 0) "PASS" else "FAIL", n = length(px), nm = sum(px > 0),
             mad = max(d), mrd = NA, tol = "pixel-identical",
             detail = if (frac == 0) "pixel-identical (different file encoding)"
                      else sprintf("%.4f%% of pixels differ; inspect visually", 100 * frac))
}

compare_one <- function(id, rep, ref, type) {
  if (!file.exists(ref)) return(result_row(id, rep, ref, type, "BLOCKED", detail = "reference file missing"))
  if (!file.exists(rep)) return(result_row(id, rep, ref, type, "FAIL", detail = "reproduced file was not produced"))
  tryCatch(switch(type,
                  table = compare_table(id, rep, ref),
                  grid  = compare_grid(id, rep, ref),
                  png   = compare_png(id, rep, ref),
                  stop("unknown comparison type ", type)),
           error = function(e) result_row(id, rep, ref, type, "FAIL", detail = paste("comparison error:", conditionMessage(e))))
}
