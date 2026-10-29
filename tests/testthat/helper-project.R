# Shared helpers for the test suite (loaded automatically by testthat).

find_root <- function() {
  d <- normalizePath(getwd(), winslash = "/")
  repeat {
    if (file.exists(file.path(d, "CHIK_benefit_risk_modelling.Rproj"))) return(d)
    p <- dirname(d)
    if (identical(p, d)) stop("project root not found from ", getwd())
    d <- p
  }
}
ROOT <- find_root()

source(file.path(ROOT, "tools/repro/config.R"), local = TRUE)
source(file.path(ROOT, "tools/repro/inputs.R"), local = TRUE)
source(file.path(ROOT, "tools/repro/compare.R"), local = TRUE)
source(file.path(ROOT, "tools/repro/manifest.R"), local = TRUE)

# Evaluate only the top-level *function definitions* of a script (no data
# loading, no simulation), so tests exercise the real published code.
load_functions_from <- function(script, names) {
  exprs <- parse(file.path(ROOT, script), keep.source = FALSE)
  env <- new.env(parent = globalenv())
  for (e in exprs) {
    if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) &&
        is.name(e[[2]]) && as.character(e[[2]]) %in% names &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
      eval(e, env)
    }
  }
  missing <- setdiff(names, ls(env))
  if (length(missing)) stop("functions not found in ", script, ": ", paste(missing, collapse = ", "))
  env
}

rscript_bin <- function() file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
