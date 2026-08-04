# =============================================================================
# 00_run_master.R
#
# Run the CHIK benefit-risk pipeline end-to-end.
#
# Project: C:/Users/user/OneDrive/CHIK_benefit_risk
#          (open CHIK_benefit_risk_modelling.Rproj in RStudio)
#
# Usage (RStudio):
#   source("02_Scripts/00_run_master.R")
#
# Usage (terminal):
#   Rscript 02_Scripts/00_run_master.R
#
# Pipeline order:
#   01_setup.R                        -> packages, data, helper functions
#   02_setup_age_props.R              -> age-specific chronic-stage proportions
#   03_brazil_all_draws_ori_v3.R      -> draw-level BRR analysis (slow)
#   04_brr_cloud.R                    -> cloud plots
#   05_brazil_attack_rate_posterior.R -> attack-rate classification
#   06_brazil_travel_final.R          -> travel BRR analysis
#   07_brazil_travel_final_update.R   -> travel BRR update
# =============================================================================

# ---- Project root -----------------------------------------------------------
find_project_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(20)) {
    if (dir.exists(file.path(d, "01_Data")) && dir.exists(file.path(d, "02_Scripts"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  NA_character_
}

detect_project_root <- function() {
  script_path <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE),
    error = function(e) NA_character_
  )
  if (!is.na(script_path)) {
    return(normalizePath(file.path(dirname(script_path), ".."), winslash = "/"))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(normalizePath(file.path(dirname(script_path), ".."), winslash = "/"))
  }

  root <- find_project_root()
  if (!is.na(root)) return(root)

  normalizePath(getwd(), winslash = "/")
}

PROJECT_ROOT <- detect_project_root()
Sys.setenv(CHIK_PROJECT_ROOT = PROJECT_ROOT)
setwd(PROJECT_ROOT)

message("Project root: ", PROJECT_ROOT)

SCRIPTS_DIR <- file.path(PROJECT_ROOT, "02_Scripts")

# ---- Step toggles (set FALSE to skip) ---------------------------------------
RUN_SETUP              <- TRUE
RUN_AGE_PROPS          <- TRUE
RUN_BRR_DRAWS          <- TRUE   # slow
RUN_CLOUD_PLOTS        <- TRUE
RUN_ATTACK_RATE        <- TRUE
RUN_TRAVEL_FINAL       <- TRUE
RUN_TRAVEL_FINAL_UPDATE<- TRUE

# ---- Runner -----------------------------------------------------------------
run_step <- function(label, script_file, enabled = TRUE) {
  if (!enabled) {
    message("[SKIP] ", label)
    return(invisible(NULL))
  }
  path <- file.path(SCRIPTS_DIR, script_file)
  if (!file.exists(path)) {
    stop("Script not found: ", path, call. = FALSE)
  }
  message("\n", strrep("=", 72))
  message("[RUN] ", label, " (", script_file, ")")
  message(strrep("=", 72))
  t_start <- Sys.time()
  source(path, local = FALSE, encoding = "UTF-8")
  elapsed <- difftime(Sys.time(), t_start, units = "mins")
  message("[DONE] ", label, " (", round(as.numeric(elapsed), 2), " min)")
  invisible(elapsed)
}

run_step("Setup",                   "01_setup.R",                        RUN_SETUP)
run_step("Age proportions",         "02_setup_age_props.R",              RUN_AGE_PROPS)
run_step("BRR draws (v3)",          "03_brazil_all_draws_ori_v3.R",      RUN_BRR_DRAWS)
run_step("BRR cloud plots",         "04_brr_cloud.R",                    RUN_CLOUD_PLOTS)
run_step("Attack rate posterior",   "05_brazil_attack_rate_posterior.R", RUN_ATTACK_RATE)
run_step("Travel BRR final",        "06_brazil_travel_final.R",          RUN_TRAVEL_FINAL)
run_step("Travel BRR update",       "07_brazil_travel_final_update.R",   RUN_TRAVEL_FINAL_UPDATE)

message("\nPipeline complete.")
