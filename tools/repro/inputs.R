# =============================================================================
# tools/repro/inputs.R
#
# Input and package checks used by the preflight (Mode A), by Mode B before
# each script is run, and by the tests. A missing or malformed input is
# reported with its role and a pointer to docs/DATA_REQUIREMENTS.md; nothing
# is ever substituted or generated in its place.
# =============================================================================

# R objects (and, for data frames, columns) that the publication scripts use
# from each input. Checked by `--mode preflight --deep`.
REQUIRED_OBJECTS <- list(
  "setting_key.RData" = list(setting_key = character(0)),
  "pop_by_state.RData" = list(pop_by_state = c("region", "tot_pop")),
  "state_ar_hazard_summary_finite.RData" = list(ar_summary = c("region", "ar_med", "ar_lo", "ar_hi")),
  "psa_grid_bra_travel_finite_setting.RData" = list(psa_grid_travel_setting = c(
    "age_group", "duration", "entry_day", "setting", "AR_median", "AR_lo", "AR_hi",
    "brr_sae_prob_gt1", "brr_death_prob_gt1", "brr_daly_median", "brr_daly_lo", "brr_daly_hi",
    "brr_daly_prob_gt1", "daly_averted_median", "daly_sae_median")),
  "psa_df_bra_travel_finite.RData" = list(psa_df = c(
    "draw", "state", "days", "age_group", "averted_10k_sae", "excess_10k_sae",
    "averted_10k_death", "excess_10k_death", "daly_averted", "daly_sae")),
  "draw_level_xy_serostatus_finite.RData" = list(draw_level_xy_serostatus = c(
    "draw_id", "Region", "Coverage", "outcome", "AgeCat", "RR_seropos", "averted_10k",
    "q_seroneg_vacc", "VE_label", "setting", "x_10k_base", "x_10k_adj", "brr_base", "brr_adj")),
  "observed_2022.RData" = list(observed_all = c("Week", "Observed", "region")),
  "brr_weeksweep_summary_finite.RData" = list(brr_weeksweep_summary = c(
    "outcome", "AgeCat", "VE_label", "setting", "Coverage", "week", "brr_base_prob_gt1", "brr_adj_prob_gt1")),
  "brr_weeksweep_summary_finite_ve0.RData" = list(brr_weeksweep_summary = c(
    "outcome", "AgeCat", "VE_label", "setting", "Coverage", "week", "brr_base_prob_gt1", "brr_adj_prob_gt1")),
  "posterior_finite_all.RData" = list("^posterior_" = NULL)   # regex: >=1 object named posterior_*
)

# Validate one .RData file against REQUIRED_OBJECTS. Returns "" if OK,
# otherwise a description of what is missing.
validate_rdata <- function(path, spec = REQUIRED_OBJECTS[[basename(path)]]) {
  if (!file.exists(path)) return("file not found")
  e <- new.env()
  loaded <- tryCatch(load(path, envir = e), error = function(err) structure(conditionMessage(err), class = "err"))
  if (inherits(loaded, "err")) return(paste("cannot load:", loaded))
  if (is.null(spec)) return("")
  problems <- character(0)
  for (obj in names(spec)) {
    if (startsWith(obj, "^")) {
      if (!any(grepl(obj, loaded))) problems <- c(problems, sprintf("no object matching %s", obj))
      next
    }
    if (!obj %in% loaded) { problems <- c(problems, sprintf("object '%s' missing", obj)); next }
    cols <- spec[[obj]]
    if (length(cols)) {
      miss <- setdiff(cols, names(e[[obj]]))
      if (length(miss)) problems <- c(problems, sprintf("'%s' lacks column(s): %s", obj, paste(miss, collapse = ", ")))
      if (is.data.frame(e[[obj]]) && nrow(e[[obj]]) == 0) problems <- c(problems, sprintf("'%s' has 0 rows", obj))
    }
  }
  paste(problems, collapse = "; ")
}

resolve_input <- function(inp, cfg) {
  switch(inp$area,
         data = file.path(cfg$data_dir, inp$file),
         orv  = file.path(cfg$orv_root, inp$file),
         file.path(cfg$root, inp$file))
}

check_inputs <- function(items, cfg, deep = FALSE) {
  rows <- list()
  for (it in items) for (inp in it$inputs) {
    path <- resolve_input(inp, cfg)
    ok <- file.exists(path)
    detail <- if (!ok) sprintf("MISSING (%s). Needed for %s.", inp$role, it$id) else ""
    if (ok && deep && grepl("\\.RData$", path)) {
      v <- validate_rdata(path)
      if (nzchar(v)) { ok <- FALSE; detail <- paste("MALFORMED:", v) } else detail <- "objects/columns OK"
    }
    rows[[length(rows) + 1]] <- data.frame(
      item = it$id, area = inp$area, file = inp$file, path = path,
      status = if (ok) "OK" else if (file.exists(path)) "MALFORMED" else "MISSING",
      size_mb = if (file.exists(path)) round(file.size(path) / 1e6, 1) else NA,
      detail = detail, stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

# Packages used by the publication scripts (static scan, see
# docs/SCRIPT_INVENTORY.md) plus those the reproduction tooling needs.
REQUIRED_PACKAGES <- data.frame(
  package = c("dplyr", "tidyr", "purrr", "tibble", "ggplot2", "patchwork", "scales", "RColorBrewer",
              "mgcv", "sf", "rnaturalearth", "rnaturalearthhires", "officer", "flextable", "writexl", "readxl",
              "png", "testthat",
              "lhs", "truncnorm", "tictoc", "tidyverse", "cowplot", "ggpubr", "reshape", "knitr",
              "kableExtra", "glue", "ggpattern", "pacman", "conflicted"),
  needed_for = c(rep("figures/tables", 16), "comparison", "tests",
                 rep("full pipeline only", 13)),
  stringsAsFactors = FALSE
)

check_packages <- function(pk = REQUIRED_PACKAGES) {
  ip <- utils::installed.packages()
  pk$installed <- pk$package %in% rownames(ip)
  pk$version <- ifelse(pk$installed, ip[match(pk$package, rownames(ip)), "Version"], NA)
  pk
}
