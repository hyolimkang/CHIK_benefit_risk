# =============================================================================
# tools/repro/sandbox.R
#
# File-I/O sandbox used when the reproduction runner executes a scientific
# script. The scientific scripts are run UNMODIFIED; instead, inside the
# (throwaway) R process that runs one script, the file-reading and
# file-writing functions are wrapped with base::trace() so that:
#
#   * reads  of the author's hard-coded locations are translated to the
#            configured locations (CHIK_DATA_DIR, CHIK_ORV_ROOT); if the same
#            file was already written earlier in this reproduction run, the
#            reproduced copy is read instead (so caches and chained scripts
#            use reproduced, not archived, intermediates);
#   * writes are redirected into CHIK_OUTPUT_DIR (default reproduced_outputs/),
#            mirroring the original relative path, so archived results in
#            01_Data/, 06_Results/, 07_Final_Results/ and the submission
#            folder are never overwritten;
#   * every redirected read and write is logged to
#            <CHIK_OUTPUT_DIR>/logs/io_log.csv.
#
# Only file *locations* are changed. No data values, arguments or scientific
# code paths are altered. Wrapped functions:
#   reads : base::load, base::readRDS, utils::read.table (read.csv),
#           readxl::read_excel
#   writes: base::save, base::saveRDS, utils::write.table (write.csv),
#           ggplot2::ggsave, writexl::write_xlsx, officer print.rdocx(target),
#           grDevices::png / pdf / jpeg
#   wd    : base::setwd (legacy absolute roots are translated)
# =============================================================================

local({
  here <- if (exists("CHIK_REPRO_TOOLS")) CHIK_REPRO_TOOLS else "tools/repro"
  source(file.path(here, "config.R"), local = TRUE)

  cfg <- chik_config()
  env <- new.env()
  env$cfg <- cfg
  env$script <- Sys.getenv("CHIK_REPRO_SCRIPT", "interactive")
  env$log_file <- file.path(cfg$out_dir, "logs", "io_log.csv")
  dir.create(dirname(env$log_file), recursive = TRUE, showWarnings = FALSE)

  starts_with_ci <- function(x, prefix) {
    nchar(x) >= nchar(prefix) && tolower(substr(x, 1, nchar(prefix))) == tolower(prefix)
  }

  env$canonical <- function(path) {
    p <- gsub("\\\\", "/", path)
    for (nm in names(LEGACY_ROOTS)) {
      lr <- LEGACY_ROOTS[[nm]]
      if (starts_with_ci(p, lr)) {
        target <- if (nm == "project") cfg$root else cfg$orv_root
        return(paste0(target, substr(p, nchar(lr) + 1, nchar(p))))
      }
    }
    if (!grepl("^([A-Za-z]:)?/", p)) p <- file.path(gsub("\\\\", "/", getwd()), p)
    p
  }

  # -> list(area = "project"/"orv"/"out"/NA, rel = path relative to that area)
  env$locate <- function(abs) {
    for (area in c("out", "project", "orv")) {
      base <- switch(area, out = cfg$out_dir, project = cfg$root, orv = cfg$orv_root)
      if (starts_with_ci(abs, paste0(base, "/")))
        return(list(area = area, rel = substr(abs, nchar(base) + 2, nchar(abs))))
    }
    list(area = NA_character_, rel = abs)
  }

  env$sandbox_path <- function(loc) {
    if (identical(loc$area, "project")) file.path(cfg$out_dir, loc$rel)
    else if (identical(loc$area, "orv")) file.path(cfg$out_dir, "_orv", loc$rel)
    else NA_character_
  }

  # Plain cat() logger: it must not go through the wrapped write functions.
  env$log <- function(op, from, to) {
    lf <- env$log_file
    q <- function(x) paste0('"', gsub('"', "'", x), '"')
    if (!file.exists(lf)) cat("time,script,op,requested,resolved\n", file = lf)
    cat(paste(q(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), q(env$script), op, q(from), q(to), sep = ","),
        "\n", file = lf, append = TRUE, sep = "")
  }

  env$remap_read <- function(path) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path)) return(path)
    abs <- env$canonical(path)
    loc <- env$locate(abs)
    if (identical(loc$area, "out")) return(abs)
    sb <- env$sandbox_path(loc)
    out <- abs
    if (!is.na(sb) && file.exists(sb)) {
      out <- sb
    } else if (identical(loc$area, "project") && grepl("^01_Data/", loc$rel)) {
      out <- file.path(cfg$data_dir, sub("^01_Data/", "", loc$rel))
    }
    if (!identical(out, gsub("\\\\", "/", path))) env$log("read", path, out)
    out
  }

  env$remap_write <- function(path) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path)) return(path)
    abs <- env$canonical(path)
    loc <- env$locate(abs)
    out <- if (identical(loc$area, "out")) abs else env$sandbox_path(loc)
    if (is.na(out)) {
      tmp <- gsub("\\\\", "/", normalizePath(tempdir(), winslash = "/"))
      out <- if (starts_with_ci(abs, tmp)) abs else file.path(cfg$out_dir, "_other", basename(abs))
    }
    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    env$log("write", path, out)
    out
  }

  env$remap_dir <- function(dir) {
    if (!is.character(dir) || length(dir) != 1L) return(dir)
    out <- env$canonical(dir)
    if (!dir.exists(out)) stop("[sandbox] setwd() target does not exist: ", out,
                               "\n  (requested: ", dir, "). Check CHIK_PROJECT_ROOT / CHIK_ORV_ROOT.", call. = FALSE)
    out
  }

  attach(env, name = "chik_sandbox", warn.conflicts = FALSE)

  tr <- function(fun, arg, kind, ns) {
    if (!requireNamespace(ns, quietly = TRUE)) return(invisible(FALSE))
    f <- switch(kind, read = "remap_read", write = "remap_write", dir = "remap_dir")
    tracer <- bquote(if (!missing(.(as.name(arg))))
      .(as.name(arg)) <- get(.(f), envir = as.environment("chik_sandbox"))(.(as.name(arg))))
    suppressMessages(trace(fun, tracer = tracer, where = asNamespace(ns), print = FALSE))
    invisible(TRUE)
  }

  tr("setwd",       "dir",      "dir",   "base")
  tr("load",        "file",     "read",  "base")
  tr("readRDS",     "file",     "read",  "base")
  tr("read.table",  "file",     "read",  "utils")
  tr("read_excel",  "path",     "read",  "readxl")
  tr("save",        "file",     "write", "base")
  tr("saveRDS",     "file",     "write", "base")
  tr("write.table", "file",     "write", "utils")
  tr("ggsave",      "filename", "write", "ggplot2")
  tr("write_xlsx",  "path",     "write", "writexl")
  tr("png",         "filename", "write", "grDevices")
  tr("pdf",         "file",     "write", "grDevices")
  tr("jpeg",        "filename", "write", "grDevices")
  # Optional: hide archived summary caches so scripts that compute a summary
  # only "if (!file.exists(cache))" recompute it from the posterior instead.
  # A cache already re-written in this run (sandbox copy) stays visible.
  env$hidden_caches <- c("state_ar_hazard_summary_finite.RData", "state_rho_summary_finite.RData")
  # file.access() (not traced) is used instead of file.exists() to avoid
  # re-entering the traced function.
  exists_raw <- function(p) file.access(p, 0) == 0
  env$file_exists_filter <- function(...) {
    paths <- c(...)
    if (!identical(Sys.getenv("CHIK_RECOMPUTE_CACHES"), "1") || !is.character(paths)) return(NULL)
    hide <- vapply(paths, function(p) {
      if (!basename(p) %in% env$hidden_caches) return(FALSE)
      sb <- env$sandbox_path(env$locate(env$canonical(p)))
      !(length(sb) && !is.na(sb) && exists_raw(sb))
    }, TRUE)
    if (!any(hide)) return(NULL)
    env$log("hide-cache", paste(paths[hide], collapse = ";"), "treated as absent")
    unname(!hide & exists_raw(paths))
  }
  # A trace() tracer cannot make file.exists() return early, so instead a
  # masking file.exists() is placed in the attached sandbox environment. Code
  # sourced into the global environment (the scientific scripts) finds it
  # before base::file.exists; package code is unaffected. attach() copied
  # `env` above, so assign into the attached copy.
  if (identical(Sys.getenv("CHIK_RECOMPUTE_CACHES"), "1")) {
    filter <- env$file_exists_filter
    assign("file.exists", function(...) {
      r <- filter(...)
      if (is.null(r)) base::file.exists(...) else r
    }, envir = as.environment("chik_sandbox"))
  }
  if (tr("print.rdocx", "target", "write", "officer")) {
    # S3 dispatch uses the registered method, so re-register the traced copy.
    registerS3method("print", "rdocx", get("print.rdocx", envir = asNamespace("officer")),
                     envir = asNamespace("officer"))
  }
  message("[sandbox] reads : ", cfg$data_dir, " | ", cfg$orv_root)
  message("[sandbox] writes: ", cfg$out_dir)
})
