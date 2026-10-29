# =============================================================================
# tools/repro/run_script.R
#
# Runs ONE unmodified scientific script inside the I/O sandbox
# (tools/repro/sandbox.R). Called by tools/run_reproduction.R in a fresh
# R process per script, so that the sandbox wrappers and any global state
# created by the script do not leak into other scripts.
#
#   Rscript tools/repro/run_script.R 02_Scripts/19_main_table1_daly.R
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/repro/run_script.R <script.R>")

self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
CHIK_REPRO_TOOLS <- dirname(normalizePath(self, winslash = "/"))
source(file.path(CHIK_REPRO_TOOLS, "config.R"))
cfg <- chik_config()
setwd(cfg$root)

script <- args[1]
if (!file.exists(script)) stop("Script not found: ", script)
Sys.setenv(CHIK_REPRO_SCRIPT = basename(script))
source(file.path(CHIK_REPRO_TOOLS, "sandbox.R"))

# Start the script from an empty global environment, as Rscript would.
# Run state lives in the attached sandbox environment, which survives the
# rm(list = ls()) below and any rm(list = ls()) inside the scripts.
sb <- as.environment("chik_sandbox")
sb$run_script_path <- normalizePath(script, winslash = "/")
sb$run_t0 <- Sys.time()
rm(list = ls(envir = globalenv(), all.names = TRUE), envir = globalenv())
# No seed is imposed here; scripts set their own seeds where they use RNG.
with(as.environment("chik_sandbox"), {
  source(run_script_path, local = globalenv(), echo = FALSE, encoding = "UTF-8")
  message(sprintf("[run_script] %s finished in %.1f s", basename(run_script_path),
                  as.numeric(difftime(Sys.time(), run_t0, units = "secs"))))
})
