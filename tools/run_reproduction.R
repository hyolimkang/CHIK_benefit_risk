# =============================================================================
# tools/run_reproduction.R  --  single entry point for reproduction
#
#   Rscript tools/run_reproduction.R --mode preflight [--deep]
#       Mode A. Checks project layout, R packages and every input listed in
#       tools/repro/manifest.R. --deep also loads each available input and
#       checks the R objects / columns the scripts rely on (slower: some
#       inputs are >1 GB). No simulation is run.
#
#   Rscript tools/run_reproduction.R --mode outputs [--only fig1,tab2,...]
#       Mode B. Regenerates manuscript figures and tables from archived model
#       outputs by running the unmodified scientific scripts inside an I/O
#       sandbox, then compares each reproduced file with the submitted
#       reference. Writes reproduced_outputs/validation_report.{csv,md}.
#
#       Add --recompute-caches to rebuild the per-state attack-rate and
#       reporting-rate summaries from the posterior draws instead of loading
#       their archived caches (Supp. Table 1, Supp. Fig. 2).
#
#   Rscript tools/run_reproduction.R --mode full --allow-full
#       Mode C. Re-runs the upstream stages that live in THIS repository
#       (setting classification, travel and outbreak-response benefit-risk
#       draws) into reproduced_outputs/01_Data, then Mode B on top of them.
#       Hours of computation; the transmission-model fit and outbreak
#       simulations in the companion CHIK_ORV_impact project are NOT re-run.
#
#   Rscript tools/run_reproduction.R --mode smoke
#       Mode D. Fast software tests (tests/testthat): parsing, path handling,
#       unit tests of the benefit-risk functions, input validators. Needs no
#       private data.
#
# Nothing in 01_Data/, 06_Results/, 07_Final_Results/ or 03_Manuscript/ is
# written by any mode. See docs/REVIEWER_REPRODUCTION_GUIDE.md.
# =============================================================================

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
TOOLS <- dirname(normalizePath(if (is.na(self)) "tools/run_reproduction.R" else self, winslash = "/"))
source(file.path(TOOLS, "repro", "config.R"))
source(file.path(TOOLS, "repro", "manifest.R"))
source(file.path(TOOLS, "repro", "compare.R"))
source(file.path(TOOLS, "repro", "inputs.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) default else if (i < length(args) && !startsWith(args[i + 1], "--")) args[i + 1] else TRUE
}
mode <- get_arg("--mode", "preflight")
only <- get_arg("--only", NULL)
deep <- isTRUE(get_arg("--deep", FALSE))
# --recompute-caches: scripts 16/32/33 load archived summary caches
# (state_ar_hazard_summary_finite.RData, state_rho_summary_finite.RData) when
# present; with this flag the sandbox hides them so they are recomputed from
# the posterior draws.
if (isTRUE(get_arg("--recompute-caches", FALSE))) Sys.setenv(CHIK_RECOMPUTE_CACHES = "1")

cfg <- chik_config()
setwd(cfg$root)
dir.create(file.path(cfg$out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

cat("Project root : ", cfg$root, "\n",
    "Data dir     : ", cfg$data_dir, if (!dir.exists(cfg$data_dir)) "  [MISSING]", "\n",
    "ORV root     : ", cfg$orv_root, if (!dir.exists(cfg$orv_root)) "  [MISSING]", "\n",
    "Output dir   : ", cfg$out_dir, "\n",
    "R            : ", R.version.string, "\n\n", sep = "")

items <- PUBLICATION_MANIFEST
if (!is.null(only)) {
  ids <- strsplit(only, ",")[[1]]
  unknown <- setdiff(ids, vapply(items, `[[`, "", "id"))
  if (length(unknown)) stop("Unknown --only id(s): ", paste(unknown, collapse = ", "),
                            ". Valid: ", paste(vapply(items, `[[`, "", "id"), collapse = ", "))
  items <- Filter(function(it) it$id %in% ids, items)
}

# ---- Mode A ---------------------------------------------------------------
run_preflight <- function(items, deep) {
  layout_ok <- all(dir.exists(file.path(cfg$root, c("02_Scripts", "tools"))))
  pk <- check_packages()
  inp <- check_inputs(items, cfg, deep = deep)
  scripts <- unique(vapply(items, `[[`, "", "script"))
  scr <- data.frame(script = scripts, exists = file.exists(scripts),
                    parses = vapply(scripts, function(s) !inherits(try(parse(s), silent = TRUE), "try-error"), TRUE))
  cat("Project layout:", if (layout_ok) "OK" else "MISSING 02_Scripts/ or tools/", "\n\n")
  cat("R packages:\n"); print(pk, row.names = FALSE); cat("\n")
  cat("Scripts:\n"); print(scr, row.names = FALSE); cat("\n")
  cat("Inputs:\n"); print(inp[, c("item", "area", "file", "status", "size_mb", "detail")], row.names = FALSE)
  utils::write.csv(inp, file.path(cfg$out_dir, "preflight_inputs.csv"), row.names = FALSE)
  ready <- vapply(items, function(it) all(inp$status[inp$item == it$id] == "OK"), TRUE)
  cat("\nReady for Mode B:", paste(vapply(items[ready], `[[`, "", "id"), collapse = ", "), "\n")
  if (any(!ready)) cat("BLOCKED (missing inputs):", paste(vapply(items[!ready], `[[`, "", "id"), collapse = ", "),
                       "\n  -> see docs/DATA_REQUIREMENTS.md\n")
  missing_pk <- pk$package[!pk$installed & pk$needed_for != "full pipeline only"]
  if (length(missing_pk)) cat("Missing packages for Mode B:", paste(missing_pk, collapse = ", "), "\n")
  invisible(list(inputs = inp, ready = ready, packages = pk))
}

# ---- Mode B ---------------------------------------------------------------
run_script_sandboxed <- function(script, log_name) {
  log_file <- file.path(cfg$out_dir, "logs", paste0(log_name, ".log"))
  t0 <- Sys.time()
  q <- function(x) shQuote(x, type = if (.Platform$OS.type == "windows") "cmd" else "sh")
  # The child inherits these (system2(env=) is not supported on Windows).
  Sys.setenv(CHIK_PROJECT_ROOT = cfg$root, CHIK_DATA_DIR = cfg$data_dir,
             CHIK_ORV_ROOT = cfg$orv_root, CHIK_OUTPUT_DIR = cfg$out_dir)
  status <- system2(rscript, c(q(file.path(TOOLS, "repro", "run_script.R")), q(script)),
                    stdout = log_file, stderr = log_file)
  list(exit = status, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")), log = log_file)
}

run_outputs <- function(items) {
  pre <- check_inputs(items, cfg, deep = FALSE)
  rows <- list(); runs <- list()
  for (it in items) {
    missing <- pre[pre$item == it$id & pre$status != "OK", ]
    if (nrow(missing)) {
      msg <- paste0("missing input(s): ", paste(missing$file, collapse = ", "), " (see docs/DATA_REQUIREMENTS.md)")
      cat(sprintf("[BLOCKED] %-9s %s\n", it$id, msg))
      for (ck in it$checks) rows[[length(rows) + 1]] <- result_row(it$id, ck$reproduced, ck$reference, ck$type, "BLOCKED", detail = msg)
      runs[[it$id]] <- data.frame(id = it$id, script = it$script, exit = NA, seconds = NA, log = "")
      next
    }
    cat(sprintf("[RUN]     %-9s %s ... ", it$id, it$script)); flush.console()
    r <- run_script_sandboxed(it$script, it$id)
    cat(sprintf("exit %d, %.0f s\n", r$exit, r$secs))
    runs[[it$id]] <- data.frame(id = it$id, script = it$script, exit = r$exit, seconds = round(r$secs, 1), log = r$log)
    for (ck in it$checks) {
      res <- compare_one(it$id, file.path(cfg$out_dir, ck$reproduced), file.path(cfg$root, ck$reference), ck$type)
      if (r$exit != 0) res$detail <- paste0("script exited with status ", r$exit, " (", r$log, "); ", res$detail)
      rows[[length(rows) + 1]] <- res
    }
  }
  report <- do.call(rbind, rows)
  report$reproduced <- sub(paste0("^", cfg$out_dir, "/"), "reproduced_outputs/", report$reproduced)
  report$reference <- sub(paste0("^", cfg$root, "/"), "", report$reference)
  runs <- do.call(rbind, runs)
  write_reports(report, runs)
  report
}

write_reports <- function(report, runs) {
  utils::write.csv(report, file.path(cfg$out_dir, "validation_report.csv"), row.names = FALSE)
  utils::write.csv(runs, file.path(cfg$out_dir, "run_log.csv"), row.names = FALSE)
  fmt <- function(x) ifelse(is.na(x), "", format(signif(x, 3)))
  md <- c("# Validation report (Mode B: reproduction from archived model outputs)", "",
          paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "| R", getRversion(),
                "| platform", R.version$platform), "",
          "| Item | Status | Reproduced | Reference | Compared | Mismatches | Max abs diff | Tolerance | Detail |",
          "|---|---|---|---|---|---|---|---|---|",
          sprintf("| %s | **%s** | `%s` | `%s` | %s | %s | %s | %s | %s |", report$id, report$status, report$reproduced,
                  report$reference, fmt(report$n_compared), fmt(report$n_mismatch), fmt(report$max_abs_diff),
                  report$tolerance, gsub("\\|", "/", report$detail)),
          "", "## Script runs", "", "| Item | Script | Exit | Seconds |", "|---|---|---|---|",
          sprintf("| %s | `%s` | %s | %s |", runs$id, runs$script, runs$exit, runs$seconds))
  writeLines(md, file.path(cfg$out_dir, "validation_report.md"), useBytes = TRUE)
  cat("\nSummary:", paste(names(table(report$status)), table(report$status), sep = "=", collapse = "  "), "\n")
  cat("Report :", file.path(cfg$out_dir, "validation_report.md"), "\n")
}

# ---- Mode C ---------------------------------------------------------------
run_full <- function() {
  if (!isTRUE(get_arg("--allow-full", FALSE)))
    stop("Mode C re-runs the upstream benefit-risk simulations (several hours). ",
         "Re-run with --allow-full to confirm. The transmission-model fit and outbreak ",
         "simulations in CHIK_ORV_impact are not re-run by this repository.", call. = FALSE)
  stages <- c("02_Scripts/08_attack_rate_classification_finite.R",
              "02_Scripts/06_brazil_travel_final_finite.R",
              "02_Scripts/03_brazil_all_draws_ori_v3.R")
  for (s in stages) {
    cat("[FULL] ", s, " ... "); flush.console()
    r <- run_script_sandboxed(s, paste0("full_", tools::file_path_sans_ext(basename(s))))
    cat(sprintf("exit %d, %.0f min\n", r$exit, r$secs / 60))
    if (r$exit != 0) stop("Upstream stage failed: ", s, " (log: ", r$log, ")", call. = FALSE)
  }
  # Downstream scripts now read the reproduced intermediates written to
  # reproduced_outputs/01_Data (the sandbox prefers reproduced copies).
  run_outputs(PUBLICATION_MANIFEST)
}

# ---- Mode D ---------------------------------------------------------------
run_smoke <- function() {
  if (!requireNamespace("testthat", quietly = TRUE)) stop("Package 'testthat' is required for --mode smoke.")
  res <- testthat::test_dir(file.path(cfg$root, "tests", "testthat"), reporter = "summary", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  failed <- sum(df$failed) + sum(df$error)
  cat(sprintf("\nTests: %d expectations, %d failed, %d skipped\n", sum(df$nb), failed, sum(df$skipped)))
  if (failed > 0) quit(status = 1)
}

switch(mode,
       preflight = invisible(run_preflight(items, deep)),
       outputs   = { pf <- run_preflight(items, FALSE); cat("\n"); run_outputs(items) },
       full      = run_full(),
       smoke     = run_smoke(),
       stop("Unknown --mode '", mode, "'. Use preflight, outputs, full or smoke."))
