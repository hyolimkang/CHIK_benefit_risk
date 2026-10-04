# =============================================================================
# SECTION 00. Script map (debugging index)
# =============================================================================
# 00) Script map and section index
# 01) Global setup and shared inputs
# 02) Build baseline draw objects (symptomatic burden)
# 03) Shared helper functions (rho scaling)
# 04) Hospitalisation draws
# 05) Fatal draws
# 06) DALY draws  <-- v3: AGE-SPECIFIC chronic proportions on benefit side
# 07) SAE draws and persistence of core draw objects
# 08) BRR pipeline (draw-level benefit vs risk; no object mixing)
#
# NOTE (v3 vs v2):
# The ONLY difference from v2 is in SECTION 06 (DALY draws / benefit side):
#   * make_daly_draws            -> make_daly_draws_age
#   * calc_total_daly_draws      -> calc_total_daly_draws_age
#   * calc_total_daly_draws_rho  -> calc_total_daly_draws_rho_age
#   * all_draws_daly_true now takes age-specific subac_prop / chr_prop vectors
#     (length = n_age_bins = 20) and applies u40 props to bins 1-10 (ages 0-39)
#     and o40 props to bins 11-20 (ages 40-89). This matches the SAE-side
#     age-specific treatment already in v2 via compute_daly_one_age_specific(),
#     restoring benefit/risk consistency.
#
# Prerequisites: source('02_Scripts/setup.R') and
#                source('02_Scripts/setup_age_props.R') before this script.
# =============================================================================

# This script's own header says: "Prerequisites: source('02_Scripts/01_setup.R')
# and source('02_Scripts/02_setup_age_props.R') before this script." -- true
# in the user's normal long-lived interactive session (both already sourced
# earlier that day) but not in a fresh Rscript process. Needed here purely
# for their FUNCTION definitions (compute_daly_one(), compute_daly_one_age_
# specific()) -- both scripts also rebuild lhs_sample as a side effect, but
# that gets correctly overwritten by this script's own load("01_Data/
# lhs_sample.RData") a few lines below, so the cached (staleness-guard-
# checked) version is still what's actually used for data.
#
# MUST run before anything else in this file is defined: 01_setup.R's own
# first line is rm(list = ls(all = TRUE)), which wipes the entire global
# environment -- defining step_msg()/step_done() (or loading libraries)
# before this would just have them erased again.
message("[03_v3_4chains] Sourcing 01_setup.R + 02_setup_age_props.R for function defs (compute_daly_one*) ...")
source("02_Scripts/01_setup.R")
source("02_Scripts/02_setup_age_props.R")

suppressMessages({
  library(tidyverse)
  library(writexl)
  library(flextable)
  library(officer)
  library(conflicted)
})

# ---- Progress messages (console) -------------------------------------------
.step_timer <- NULL
step_msg <- function(label, detail = NULL) {
  stamp <- format(Sys.time(), "%H:%M:%S")
  line <- if (is.null(detail)) label else paste0(label, " — ", detail)
  message(sprintf("[03_v3 %s] START %s", stamp, line))
  .step_timer <<- Sys.time()
}

step_done <- function(label) {
  elapsed <- if (!is.null(.step_timer)) {
    as.numeric(difftime(Sys.time(), .step_timer, units = "secs"))
  } else {
    NA_real_
  }
  stamp <- format(Sys.time(), "%H:%M:%S")
  if (is.na(elapsed)) {
    message(sprintf("[03_v3 %s] DONE   %s", stamp, label))
  } else {
    message(sprintf("[03_v3 %s] DONE   %s (%.1f sec)", stamp, label, elapsed))
  }
  .step_timer <<- NULL
}

step_msg("03_brazil_all_draws_ori_v3_4chains.R", "script started")

# =============================================================================
# 4-CHAIN COMPARISON VARIANT (chat record 2026-09-06)
# Byte-identical to 03_brazil_all_draws_ori_v3.R except: (1) the "flat" data
# loads below point at the *_4chains inputs (built from the 4-chain/5000-iter
# Stan refit, 4,000-draw subsample -- see age_struc_fitting_region_finite_
# 2022_4chains.R / lhs_samples_sir_finite_4chains.R / preui_flat_run_4chains.R
# / model_0803_finite_history_4chains.R), and (2) MODE_SUFFIX for "flat" is
# "_finite_4chains" instead of "_finite", so every save() in this script
# (already parameterised on MODE_SUFFIX throughout) writes new files instead
# of overwriting the 1-chain production outputs. Purpose: run the SAME BRR
# pipeline on both posteriors so the two can be compared draw-for-draw
# (e.g. brr_table_final_long_finite.RData vs _finite_4chains.RData) --
# not a permanent replacement.
# =============================================================================

# Ensure 01_Data paths resolve when sourced via full path
if (!dir.exists("01_Data") && dir.exists(file.path("..", "01_Data"))) {
  setwd(normalizePath(".."))
  message("[03_v3] setwd to project root: ", getwd())
}

# ---- Data layer toggle -------------------------------------------------------
# DATA_MODE selects which simulation dataset (FOI / immunity assumption) to load:
#   "original" : sim_results_vc_ixchiq_model       (long-term average FOI; R0 = 1 - exp(-FOI * age))
#   "flat"     : sim_results_vc_ixchiq_model_flat  (flat immunity, 12-yr exposure; R0 = 1 - exp(-FOI * 12))
DATA_MODE <- "flat"
#DATA_MODE <- "original"

stopifnot(DATA_MODE %in% c("original", "flat"))

# Suffix applied to all save() output filenames below, so re-running under a
# different DATA_MODE never silently overwrites another mode's results.
MODE_SUFFIX <- switch(DATA_MODE,
                       "original" = "_longterm",
                       "flat"     = "_finite_4chains")

# REGIONS_TO_RUN (applies to "original" and "flat", which hold all 11 states):
#   NULL        = analyse all 11 regions
#   c("Bahia")  = subset tables/loops to Bahia only
REGIONS_TO_RUN <- NULL

subset_named_list <- function(x, keep_names) {
  missing <- setdiff(keep_names, names(x))
  if (length(missing)) {
    stop(
      "REGIONS_TO_RUN not found: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  x[keep_names]
}

apply_regions_filter <- function(regions = REGIONS_TO_RUN) {
  if (is.null(regions)) return(invisible(NULL))
  message(
    "[03_v3] Region filter ON (flat data) — analysing only: ",
    paste(regions, collapse = ", ")
  )
  postsim_vc_ixchiq_model <<- subset_named_list(postsim_vc_ixchiq_model, regions)
  if (exists("sim_results_vc_ixchiq_model", envir = .GlobalEnv)) {
    sim_results_vc_ixchiq_model <<-
      subset_named_list(sim_results_vc_ixchiq_model, regions)
  }
  if (exists("combined_nnv_df_region_coverage_model", envir = .GlobalEnv) &&
      "region" %in% names(combined_nnv_df_region_coverage_model)) {
    combined_nnv_df_region_coverage_model <<-
      dplyr::filter(combined_nnv_df_region_coverage_model, region %in% regions)
  }
  invisible(NULL)
}

apply_regions_filter_post_setup <- function(regions = REGIONS_TO_RUN) {
  if (is.null(regions)) return(invisible(NULL))
  posterior_list <<- subset_named_list(posterior_list, regions)
  preui_all <<- subset_named_list(preui_all, regions)
  invisible(NULL)
}

# =============================================================================
# SECTION 01. Global setup and shared inputs
# =============================================================================
step_msg("SECTION 01", "Global setup and shared inputs")

# =============================================================================
# current wd:  "C:/Users/user/OneDrive/CHIK_benefit_risk"
# for the flat immunity assumption 
# open: 
# preui_flat_all.RData
# posterior_flat_all.RData
# sim results flat
# postsim flat
# hosp
# fatal
# le by age
# =============================================================================
load("01_Data/chikv_fatal_hosp_rate.RData")
le_sample <- readRDS("01_Data/le_sample.RDS")
lhs_sample_young <- readRDS("01_Data/lhs_sample_young.RDS")
lhs_old <- readRDS("01_Data/lhs_old.RDS")
load("01_Data/lhs_orv.RData")
load("01_Data/pop_by_state.RData")
load("01_Data/lhs_sample.RData")

if (DATA_MODE == "flat") {
  message("[03_v3_4chains] Data mode: flat, 4-CHAIN inputs (11 states, flat immunity, R0 = 1 - exp(-FOI * (2022-2014)))")
  load("01_Data/preui_flat_all_4chains.RData")
  load("01_Data/posterior_finite_all_4chains.RData")
  load("01_Data/sim_results_vc_ixchiq_model_finite_4chains.RData")
  load("01_Data/postsim_vc_ixchiq_model_finite_4chains.RData")
  load("01_Data/combined_nnv_df_region_coverage_model_finite_4chains.RData")
  if (!is.null(REGIONS_TO_RUN)) {
    apply_regions_filter(REGIONS_TO_RUN)
  }
} else {   "original"
  message("[03_v3] Data mode: original (11 states, long-term avg FOI, R0 = 1 - exp(-FOI * age))")
  load("01_Data/preui_longterm_all.RData")             # individual preui_* objects (11 regions); combined into preui_all below
  load("01_Data/posterior_longterm_all.RData")    # shared pre-vacc posterior fit (rho/base_beta/I0); renamed from posterior_flat_all.RData -- fit with unbounded long-term FOI sero, matches this branch's R0 = 1-exp(-FOI*age)
  load("01_Data/sim_results_vc_ixchiq_model_longterm.RData")
  load("01_Data/postsim_vc_ixchiq_model_longterm.RData")
  load("01_Data/combined_nnv_df_region_coverage_model_longterm.RData")
  if (!is.null(REGIONS_TO_RUN)) {
    apply_regions_filter(REGIONS_TO_RUN)
  }
}

# lhs_sample needs *_o40 columns from 02_setup_age_props.R (v3 DALY benefit side).
# Do NOT blindly reload lhs_sample.RData if 02 already extended lhs_sample in memory.
AGE_PROPS_O40_COLS <- c("subac_o40", "chr6m_o40", "chr12m_o40", "chr30m_o40")

lhs_has_age_props <- function() {
  exists("lhs_sample", envir = globalenv()) &&
    is.data.frame(lhs_sample) &&
    all(AGE_PROPS_O40_COLS %in% colnames(lhs_sample))
}

# Staleness guard (added 2026-08-31, chat record): 01_Data/lhs_sample.RData is
# a CACHE of 01_setup.R's randomLHS() output, not regenerated automatically
# when 01_setup.R changes. It silently went stale from 2026-06-15 to
# 2026-08-27 (01_setup.R's LHS dimensionality changed from k=62, including 11
# now-removed per-state ar_* columns, to k=51 -- which shifts every column's
# values, not just the removed ones, since lhs::randomLHS()'s design depends
# on total k) -- draw_level_xy_serostatus_finite.RData/Main_Table_1_*/
# headline_summary_outbreak/headline_benefit_risk_plane_* were all built from
# that stale cache. Fail loudly instead of silently loading an out-of-sync
# cache again: regenerate via 02_Scripts/regenerate_lhs_sample_cache.R
# whenever 01_setup.R or 02_setup_age_props.R changes.
STALE_STATE_AR_COLS <- c("ar_ce", "ar_bh", "ar_pa", "ar_pn", "ar_rg", "ar_pi", "ar_tc", "ar_ag", "ar_mg", "ar_se", "ar_go")
validate_lhs_sample_fresh <- function() {
  if (any(STALE_STATE_AR_COLS %in% colnames(lhs_sample))) {
    stop(
      "lhs_sample has the old per-state ar_* columns (k=62 design) -- this is the ",
      "2026-06-15 stale cache, not current 01_setup.R (k=51). Re-run ",
      "02_Scripts/regenerate_lhs_sample_cache.R before this script.",
      call. = FALSE
    )
  }
  if (!all(lhs_sample$p_death_vacc_u65 == 0)) {
    stop(
      "lhs_sample$p_death_vacc_u65 is not all exactly 0 -- this predates the ",
      "18-64 vaccine-attributable-death-is-impossible fix in 01_setup.R. ",
      "Re-run 02_Scripts/regenerate_lhs_sample_cache.R before this script.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

if (!lhs_has_age_props()) {
  load("01_Data/lhs_sample.RData")
}
validate_lhs_sample_fresh()

ensure_lhs_age_props <- function() {
  if (lhs_has_age_props()) return(invisible(TRUE))
  age_props_script <- file.path("02_Scripts", "02_setup_age_props.R")
  if (!file.exists(age_props_script)) {
    stop(
      "lhs_sample is missing: ",
      paste(setdiff(AGE_PROPS_O40_COLS, colnames(lhs_sample)), collapse = ", "),
      ". Run source('02_Scripts/01_setup.R') then ",
      "source('02_Scripts/02_setup_age_props.R') before this script.",
      call. = FALSE
    )
  }
  message("[03_v3] lhs_sample missing *_o40 columns — sourcing 02_setup_age_props.R")
  source(age_props_script, local = FALSE, encoding = "UTF-8")
  if (!lhs_has_age_props()) {
    stop("02_setup_age_props.R did not add required *_o40 columns.", call. = FALSE)
  }
  invisible(TRUE)
}

ensure_lhs_age_props()

if (!exists("setting_key", envir = globalenv()) &&
    file.exists("01_Data/setting_key.RData")) {
  load("01_Data/setting_key.RData")
}

age_groups <- c(mean(0:1),
                mean(1:4),
                mean(5:9),
                mean(10:11),
                mean(12:17),
                mean(18:19),
                mean(20:24),
                mean(25:29),
                mean(30:34),
                mean(35:39),
                mean(40:44),
                mean(45:49),
                mean(50:54),
                mean(55:59),
                mean(60:64),
                mean(65:69),
                mean(70:74),
                mean(75:79),
                mean(80:84),
                mean(85:89)
)

conflicted::conflicts_prefer(dplyr::filter, dplyr::rename)
step_done("SECTION 01")

# =============================================================================
# SECTION 02. Build baseline draw objects (symptomatic burden)
# =============================================================================
step_msg("SECTION 02", "Build baseline draw objects (symptomatic burden)")
# Step 2.1) Collect posterior rho draws by region
# "original" and "flat" share the same pre-vacc posterior fit (posterior_flat_all)
posterior_list <- list(
  "Ceará"               = posterior_ce,
  "Bahia"               = posterior_bh,
  "Paraíba"             = posterior_pa,
  "Pernambuco"          = posterior_pn,
  "Rio Grande do Norte" = posterior_rg,
  "Piauí"               = posterior_pi,
  "Alagoas"             = posterior_ag,
  "Tocantins"           = posterior_tc,
  "Minas Gerais"        = posterior_mg,
  "Sergipe"             = posterior_se,
  "Goiás"               = posterior_go
)

# Both preui_flat_all.RData ("flat") and preui_longterm_all.RData
# ("original") provide individual preui_* objects, not a pre-combined list
# -- build preui_all from them either way.
preui_all <- list(
  preui_ce,
  preui_ag,
  preui_bh,
  preui_go,
  preui_mg,
  preui_pa,
  preui_pi,
  preui_pn,
  preui_rg,
  preui_se,
  preui_tc
)

names(preui_all) <- c(
  "Ceará", "Alagoas", "Bahia", "Goiás",
  "Minas Gerais", "Paraíba", "Piauí",
  "Pernambuco", "Rio Grande do Norte",
  "Sergipe", "Tocantins"
)

if (!is.null(REGIONS_TO_RUN)) {
  apply_regions_filter_post_setup(REGIONS_TO_RUN)
  message("[03_v3] posterior_list / preui_all: filtered to REGIONS_TO_RUN")
}

rho_pool <- purrr::imap_dfr(posterior_list, function(post, region_name){
  tibble(
    Region = region_name,
    rho_id = seq_along(post$rho),
    rho    = as.numeric(post$rho)
  )
})

#save(posterior_list, file = "01_Data/posterior_list.RData")

# Step 2.2) Utility to aggregate draw-level pre/post burden totals
calc_total_impact_draws <- function(pre_list, post_arr) {
  pre_arr <- simplify2array(pre_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

# Step 2.3) Build draw-level symptomatic burden (unscaled)
step_msg("Step 2.3", "building all_draws_ix (slow — per region)")
all_draws_ix <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  message(sprintf("[03_v3 %s]   ... all_draws_ix region: %s", format(Sys.time(), "%H:%M:%S"), region_name))
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        draw_totals <- calc_total_impact_draws(pre_list, post_arr)
        
        mutate(draw_totals,
               Region   = region_name,
               VE       = ve_name,
               Coverage = cov_name,
               Scenario = scen_id)
      })
    })
  })
})
step_done("Step 2.3 all_draws_ix")
step_done("SECTION 02")

# =============================================================================
# SECTION 03. Shared helper functions (rho scaling)
# =============================================================================
step_msg("SECTION 03", "Shared helper functions (rho scaling -- now a no-op)")
# Step 3.1) NOTE (2026-08): sim_results_list_rawsymp / age_array_raw_symp
# already carry the TRUE symptomatic burden -- true_symptomatic =
# p_sym * new_e * (1 - VE), with NO rho baked in -- per the
# age_struc_bra_seir_m3.stan refit (rho = p_sym * rho_sym separation) and the
# matching fix in sim_functions_final.R / nnv_list() / postsim_all_ui(). These
# two helpers used to divide by rho to undo the OLD formula's baked-in rho
# (0.5242478 * new_e * (1-VE) * rho); with that rho no longer present upstream,
# dividing here would re-inflate already-true burden by ~1/rho (~7-10x). Kept
# as passthroughs (rather than removed) so every call site that "reuses the
# shared rho-scaling helpers" -- hosp/fatal/DALY draws, diagnostic
# averted_bin_one -- stays fixed via this one change.
set.seed(1)
scale_symp_by_rho_pre <- function(pre_list, rho_vec) {
  pre_list
}

scale_symp_by_rho_post <- function(post_arr, rho_vec) {
  post_arr
}


# Step 3.2) total_*_true now equals total_* (no rho division) -- rho column
# kept (unused for scaling) in case downstream code inspects it.
all_draws_ix_true <- all_draws_ix %>%
  group_by(Region) %>%
  mutate(
    rho = sample(rho_pool$rho[rho_pool$Region == first(Region)],
                 size = n(), replace = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    total_pre_true  = total_pre,
    total_post_true = total_post
  )
step_done("SECTION 03")

# =============================================================================
# SECTION 04. Hospitalisation draws
# =============================================================================
step_msg("SECTION 04", "Hospitalisation draws")
# Step 4.1) Convert symptomatic draws to hospitalisation draws
make_hosp_draws <- function(symp_list, hosp_rate) {
  lapply(symp_list, function(symp_matrix) {
    sweep(symp_matrix, 1, hosp_rate, `*`)
  })
}

calc_total_hosp_draws <- function(pre_list, post_arr, hosp_rate) {
  # Pre
  pre_hosp_list <- make_hosp_draws(pre_list, hosp_rate)
  pre_arr <- simplify2array(pre_hosp_list)
  if (length(dim(pre_arr)) == 2) pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  
  # Post
  post_list <- lapply(seq_len(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_hosp_list <- make_hosp_draws(post_list, hosp_rate)
  post_arr2 <- simplify2array(post_hosp_list)
  if (length(dim(post_arr2)) == 2) post_arr2 <- array(post_arr2, dim = c(dim(post_arr2), 1))
  
  # Align
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  pre_arr   <- pre_arr[ , , 1:n_draws, drop = FALSE]
  post_arr2 <- post_arr2[ , , 1:n_draws, drop = FALSE]
  
  # Sum
  tibble(
    draw_id    = seq_len(n_draws),
    total_pre  = apply(pre_arr,  3, sum, na.rm = TRUE),
    total_post = apply(post_arr2, 3, sum, na.rm = TRUE)
  )
}

calc_total_hosp_draws_rho <- function(pre_list, post_arr, hosp_rate, rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws, drop = FALSE]
  rho_vec  <- rho_vec[1:n_draws]
  
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  calc_total_hosp_draws(pre_list_true, post_arr_true, hosp_rate = hosp_rate)
}

set.seed(1)

# Step 4.2) Build rho-adjusted hospitalisation draw table
step_msg("Step 4.2", "building all_draws_hosp_true (slow — per region)")
all_draws_hosp_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  message(sprintf("[03_v3 %s]   ... all_draws_hosp_true region: %s", format(Sys.time(), "%H:%M:%S"), region_name))
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # region rho pool
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        n_draws <- min(length(pre_list), dim(post_arr)[3])

        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_hosp_draws_rho(
          pre_list = pre_list,
          post_arr = post_arr,
          hosp_rate = hosp,
          rho_vec = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})
step_done("Step 4.2 all_draws_hosp_true")
step_done("SECTION 04")

# =============================================================================
# SECTION 05. Fatal draws
# =============================================================================
step_msg("SECTION 05", "Fatal draws")
# Step 5.1) Split symptomatic burden into hospitalised/non-hospitalised fatal burden
make_fatal_hosp_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate) {
  lapply(symp_list, function(symp_matrix) {
    hospitalised     <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised <- symp_matrix - hospitalised
    
    fatal <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    
    list(
      hosp  = hospitalised,
      fatal = fatal
    )
  })
}

# NOTE:
# This function calculates the hospitalisation component using the same
# fatal-model decomposition pipeline. It is intentionally named differently
# from calc_total_hosp_draws(...) above to avoid accidental function override.
calc_total_hosp_from_fatal_model <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_hosp_list <- lapply(pre_conv, `[[`, "hosp")
  pre_arr <- simplify2array(pre_hosp_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq_len(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_hosp_list <- lapply(post_conv, `[[`, "hosp")
  post_arr2 <- simplify2array(post_hosp_list)
  if (length(dim(post_arr2)) == 2) {
    post_arr2 <- array(post_arr2, dim = c(dim(post_arr2), 1))
  }
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  pre_arr   <- pre_arr[ , , 1:n_draws, drop = FALSE]
  post_arr2 <- post_arr2[ , , 1:n_draws, drop = FALSE]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr2, 3, sum, na.rm = TRUE)
  
  tibble(
    draw_id    = seq_len(n_draws),
    total_pre  = total_pre,
    total_post = total_post
  )
}

# Reuse shared rho-scaling helpers defined in SECTION 03.

calc_total_fatal_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  # ── Pre ────────────────────────────────
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_fatal_list <- lapply(pre_conv, `[[`, "fatal")
  pre_arr <- simplify2array(pre_fatal_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }
  
  # ── Post ───────────────────────────────
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_fatal_list <- lapply(post_conv, `[[`, "fatal")
  post_arr <- simplify2array(post_fatal_list)
  
  # ── Align draws ────────────────────────
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]
  
  # ── Sum across all ages × weeks ────────
  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

calc_total_fatal_draws_rho <- function(pre_list, post_arr,
                                       hosp_rate, fatal_rate, nh_fatal_rate,
                                       rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  calc_total_fatal_draws(pre_list_true, post_arr_true,
                         hosp_rate, fatal_rate, nh_fatal_rate)
}


set.seed(1)

# Step 5.2) Build rho-adjusted fatal draw table
step_msg("Step 5.2", "building all_draws_fatal_true (slow — per region)")
all_draws_fatal_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  message(sprintf("[03_v3 %s]   ... all_draws_fatal_true region: %s", format(Sys.time(), "%H:%M:%S"), region_name))
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
  
  # Region-specific rho pool (required)
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
  
  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }
  
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
        
        post_arr <- scen$sim_result$age_array_raw_symp
        
        n_draws <- min(length(pre_list), dim(post_arr)[3])
        
        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
        
        draw_totals <- calc_total_fatal_draws_rho(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,
          rho_vec       = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})
step_done("Step 5.2 all_draws_fatal_true")
step_done("SECTION 05")

# =============================================================================
# SECTION 06. DALY draws
# =============================================================================
step_msg("SECTION 06", "DALY draws (age-specific benefit side)")
# Step 6.1) Convert symptomatic burden to DALY components
# v3: age-specific version. subac_prop and chr_prop are now VECTORS of length
# n_age_bins (= length of the first dimension of symp_matrix). u40 vs o40 is
# resolved OUTSIDE this function by the caller (see build_daly_prop_vecs()).
make_daly_draws_age <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate,
                                dw_hosp, dw_nonhosp, dw_chronic,
                                dur_acute, dur_subacute, dur_chronic,
                                dw_subacute,
                                subac_prop_vec, chr_prop_vec,
                                le_left_vec) {
  lapply(symp_list, function(symp_matrix) {
    n_age <- nrow(symp_matrix)
    stopifnot(
      length(subac_prop_vec) == n_age,
      length(chr_prop_vec)   == n_age,
      length(le_left_vec)    == n_age
    )

    hospitalised     <- sweep(symp_matrix, 1, hosp_rate, `*`)
    non_hospitalised <- symp_matrix - hospitalised
    fatal <- sweep(hospitalised, 1, fatal_rate, `*`) +
      sweep(non_hospitalised, 1, nh_fatal_rate, `*`)

    # Acute YLD (age-independent dw / duration, as in v2)
    yld_acute <- (hospitalised     * dw_hosp    * dur_acute) +
      (non_hospitalised * dw_nonhosp * dur_acute)

    # Subacute YLD: subac_prop applied to hospitalised rows, chr_prop to
    # non-hospitalised (preserves v2 convention; only props are now age vectors)
    yld_subacute <-
      sweep(hospitalised,     1, subac_prop_vec, `*`) * dw_subacute * dur_subacute +
      sweep(non_hospitalised, 1, chr_prop_vec,   `*`) * dw_subacute * dur_subacute

    # Chronic YLD: chr_prop applied to both, per-age
    yld_chronic <-
      sweep(hospitalised,     1, chr_prop_vec, `*`) * dw_chronic * dur_chronic +
      sweep(non_hospitalised, 1, chr_prop_vec, `*`) * dw_chronic * dur_chronic

    yld_total <- yld_acute + yld_subacute + yld_chronic
    yll       <- sweep(fatal, 1, le_left_vec, `*`)
    daly_tot  <- yld_total + yll

    list(
      hosp     = hospitalised,
      fatal    = fatal,
      yld_tot  = yld_total,
      yll      = yll,
      daly_tot = daly_tot
    )
  })
}

# Backward-compatible alias: if someone calls the old name, dispatch to the
# age-aware version by broadcasting scalars to length-n_age vectors. This
# preserves behaviour of v2 exactly when scalars are supplied.
make_daly_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate,
                            dw_hosp, dw_nonhosp, dw_chronic,
                            dur_acute, dur_subacute, dur_chronic,
                            dw_subacute, subac_prop, chr_prop,
                            le_left_vec) {
  n_age <- length(le_left_vec)
  if (length(subac_prop) == 1) subac_prop <- rep(subac_prop, n_age)
  if (length(chr_prop)   == 1) chr_prop   <- rep(chr_prop,   n_age)
  make_daly_draws_age(
    symp_list      = symp_list,
    hosp_rate      = hosp_rate,
    fatal_rate     = fatal_rate,
    nh_fatal_rate  = nh_fatal_rate,
    dw_hosp        = dw_hosp,
    dw_nonhosp     = dw_nonhosp,
    dw_chronic     = dw_chronic,
    dur_acute      = dur_acute,
    dur_subacute   = dur_subacute,
    dur_chronic    = dur_chronic,
    dw_subacute    = dw_subacute,
    subac_prop_vec = subac_prop,
    chr_prop_vec   = chr_prop,
    le_left_vec    = le_left_vec
  )
}

# -------------------------------------------------------------------------
# Helper: build subac / chr proportion vectors of length n_age.
#
# CRITICAL (v3 bug fix): v2 benefit used medians from `lhs_sample_young`
# (a SEPARATE LHS object loaded from 01_Data/lhs_sample_young.RData), NOT
# from `lhs_sample`. These two objects have different Beta draws and
# different medians, so mixing them silently would shift BRR. To preserve
# v2 behaviour exactly for u40 bins (0-39y), this helper pulls u40 medians
# from `lhs_sample_young` and only the new o40 medians from `lhs_sample`
# (which setup_age_props.R extended with *_o40 columns).
#
#   Bins 1-10 (ages 0-39)  -> u40 = median(lhs_sample_young$*)     [== v2]
#   Bins 11-20 (ages 40-89) -> o40 = median(lhs_sample$*_o40)       [v3 new]
# -------------------------------------------------------------------------
build_daly_prop_vecs <- function(n_age = 20, u40_cutoff_bin = 10L) {
  stopifnot(
    exists("lhs_sample_young", envir = globalenv()),
    exists("lhs_sample",       envir = globalenv()),
    all(c("subac_o40","chr6m_o40","chr12m_o40","chr30m_o40") %in%
          colnames(lhs_sample))
  )

  # u40 (bins 1-10): match v2 exactly -> use lhs_sample_young medians
  subac_u40 <- stats::median(lhs_sample_young$subac, na.rm = TRUE)
  chr_u40   <- stats::median(lhs_sample_young$chr6m,  na.rm = TRUE) +
    stats::median(lhs_sample_young$chr12m, na.rm = TRUE) +
    stats::median(lhs_sample_young$chr30m, na.rm = TRUE)

  # o40 (bins 11-20): new sampling from setup_age_props.R (user-supplied
  # Beta parameters for 65+/>40 population)
  subac_o40 <- stats::median(lhs_sample$subac_o40, na.rm = TRUE)
  chr_o40   <- stats::median(lhs_sample$chr6m_o40,  na.rm = TRUE) +
    stats::median(lhs_sample$chr12m_o40, na.rm = TRUE) +
    stats::median(lhs_sample$chr30m_o40, na.rm = TRUE)

  # Bin mapping (matches age_groups definition above):
  # 1=0-1, 2=1-4, 3=5-9, 4=10-11, 5=12-17,
  # 6=18-19, 7=20-24, 8=25-29, 9=30-34, 10=35-39,       <-- u40
  # 11=40-44, 12=45-49, 13=50-54, 14=55-59, 15=60-64,   <-- o40 (40-64)
  # 16=65-69, 17=70-74, 18=75-79, 19=80-84, 20=85-89    <-- o40 (65+)
  idx_o40 <- seq.int(u40_cutoff_bin + 1L, n_age)
  subac_vec <- rep(subac_u40, n_age); subac_vec[idx_o40] <- subac_o40
  chr_vec   <- rep(chr_u40,   n_age); chr_vec[idx_o40]   <- chr_o40

  list(
    subac_prop_vec = subac_vec,
    chr_prop_vec   = chr_vec,
    meta = list(
      subac_u40 = subac_u40, subac_o40 = subac_o40,
      chr_u40   = chr_u40,   chr_o40   = chr_o40,
      u40_cutoff_bin = u40_cutoff_bin,
      n_age = n_age,
      u40_source = "lhs_sample_young (matches v2)",
      o40_source = "lhs_sample$*_o40 (setup_age_props.R)"
    )
  )
}
le_by_age <- function(age_numeric) {
  if (age_numeric <= 1) return(quantile(le_sample$le_1, 0.5))
  if (age_numeric < 20) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 30) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 40) return(quantile(le_sample$le_3, 0.5))
  if (age_numeric < 50) return(quantile(le_sample$le_4, 0.5))
  if (age_numeric < 60) return(quantile(le_sample$le_5, 0.5))
  if (age_numeric < 70) return(quantile(le_sample$le_6, 0.5))
  if (age_numeric < 80) return(quantile(le_sample$le_7, 0.5))
  return(quantile(le_sample$le_8, 0.5))
}

calc_total_daly_draws <- function(pre_list, post_arr,
                                  hosp_rate, fatal_rate, nh_fatal_rate,
                                  dw_hosp, dw_nonhosp, dw_chronic,
                                  dur_acute, dur_subacute, dur_chronic,
                                  dw_subacute, subac_prop, chr_prop,
                                  le_left_vec) {
  # v3: subac_prop / chr_prop are now accepted as EITHER scalar (back-compat)
  # or length-n_age vectors (age-specific). make_daly_draws() dispatches.
  pre_conv <- make_daly_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate,
                              dw_hosp, dw_nonhosp, dw_chronic,
                              dur_acute, dur_subacute, dur_chronic,
                              dw_subacute, subac_prop, chr_prop,
                              le_left_vec)
  pre_daly_list <- lapply(pre_conv, `[[`, "daly_tot")
  pre_arr <- simplify2array(pre_daly_list)
  if (length(dim(pre_arr)) == 2) {
    pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  }

  # Post
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[ , , i])
  post_conv <- make_daly_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate,
                               dw_hosp, dw_nonhosp, dw_chronic,
                               dur_acute, dur_subacute, dur_chronic,
                               dw_subacute, subac_prop, chr_prop,
                               le_left_vec)
  post_daly_list <- lapply(post_conv, `[[`, "daly_tot")
  post_arr <- simplify2array(post_daly_list)

  # Align
  n_draws <- min(dim(pre_arr)[3], dim(post_arr)[3])
  pre_arr  <- pre_arr[ , , 1:n_draws]
  post_arr <- post_arr[ , , 1:n_draws]

  total_pre  <- apply(pre_arr,  3, sum, na.rm = TRUE)
  total_post <- apply(post_arr, 3, sum, na.rm = TRUE)
  
  tibble(draw_id = 1:n_draws,
         total_pre = total_pre,
         total_post = total_post)
}

le_left_vec <- sapply(age_groups, le_by_age)

# Reuse shared rho-scaling helpers defined in SECTION 03.

calc_total_daly_draws_rho <- function(pre_list, post_arr,
                                      hosp_rate, fatal_rate, nh_fatal_rate,
                                      dw_hosp, dw_nonhosp, dw_chronic,
                                      dur_acute, dur_subacute, dur_chronic,
                                      dw_subacute, subac_prop, chr_prop,
                                      le_left_vec,
                                      rho_vec) {
  
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]
  post_arr <- post_arr[,,1:n_draws]
  rho_vec  <- rho_vec[1:n_draws]
  
  # 2) rho scaling
  pre_list_true <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_arr_true <- scale_symp_by_rho_post(post_arr, rho_vec)
  
  calc_total_daly_draws(
    pre_list = pre_list_true,
    post_arr = post_arr_true,
    hosp_rate     = hosp_rate,
    fatal_rate    = fatal_rate,
    nh_fatal_rate = nh_fatal_rate,
    dw_hosp       = dw_hosp,
    dw_nonhosp    = dw_nonhosp,
    dw_chronic    = dw_chronic,
    dur_acute     = dur_acute,
    dur_subacute  = dur_subacute,
    dur_chronic   = dur_chronic,
    dw_subacute   = dw_subacute,
    subac_prop    = subac_prop,
    chr_prop      = chr_prop,
    le_left_vec   = le_left_vec
  )
}


set.seed(1)

# v3: build age-specific subac / chr proportion vectors ONCE (uses global
# lhs_sample extended by 02_Scripts/setup_age_props.R).
.daly_prop_vecs <- build_daly_prop_vecs(n_age = length(le_left_vec),
                                        u40_cutoff_bin = 10L)
.subac_prop_vec <- .daly_prop_vecs$subac_prop_vec
.chr_prop_vec   <- .daly_prop_vecs$chr_prop_vec

message("[v3] DALY benefit age-specific props built:")
message(sprintf("  subac u40=%.4f  o40=%.4f", .daly_prop_vecs$meta$subac_u40,
                .daly_prop_vecs$meta$subac_o40))
message(sprintf("  chr   u40=%.4f  o40=%.4f", .daly_prop_vecs$meta$chr_u40,
                .daly_prop_vecs$meta$chr_o40))
message(sprintf("  u40 bins = 1..%d (ages 0-39) [source: %s]",
                .daly_prop_vecs$meta$u40_cutoff_bin,
                .daly_prop_vecs$meta$u40_source))
message(sprintf("  o40 bins = %d..%d (ages 40+) [source: %s]",
                .daly_prop_vecs$meta$u40_cutoff_bin + 1L,
                .daly_prop_vecs$meta$n_age,
                .daly_prop_vecs$meta$o40_source))

# Sanity check: compare with v2 scalar medians (should match u40 values exactly)
.v2_subac_med <- stats::median(lhs_sample_young$subac, na.rm = TRUE)
.v2_chr_med   <- stats::median(lhs_sample_young$chr6m,  na.rm = TRUE) +
  stats::median(lhs_sample_young$chr12m, na.rm = TRUE) +
  stats::median(lhs_sample_young$chr30m, na.rm = TRUE)
message(sprintf("[v3] Sanity vs v2: subac_u40 - v2_median = %.6f  (should be 0)",
                .daly_prop_vecs$meta$subac_u40 - .v2_subac_med))
message(sprintf("[v3] Sanity vs v2: chr_u40   - v2_median = %.6f  (should be 0)",
                .daly_prop_vecs$meta$chr_u40 - .v2_chr_med))

# Step 6.2) Build rho-adjusted DALY draw table (age-specific benefit side)
step_msg("Step 6.2", "building all_draws_daly_true (slow — per region)")
all_draws_daly_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  message(sprintf("[03_v3 %s]   ... all_draws_daly_true region: %s", format(Sys.time(), "%H:%M:%S"), region_name))
  pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp

  # region rho pool
  rho_pool_region <- posterior_list[[region_name]]$rho
  rho_pool_region <- as.numeric(rho_pool_region)
  rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]

  if (length(rho_pool_region) == 0) {
    stop("No valid rho draws for region: ", region_name)
  }

  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {

        post_arr <- scen$sim_result$age_array_raw_symp

        n_draws <- min(length(pre_list), dim(post_arr)[3])

        rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)

        draw_totals <- calc_total_daly_draws_rho(
          pre_list, post_arr,
          hosp_rate     = hosp,
          fatal_rate    = fatal,
          nh_fatal_rate = nh_fatal,

          dw_hosp       = quantile(lhs_sample_young$dw_hosp, 0.5),
          dw_nonhosp    = quantile(lhs_sample_young$dw_nonhosp, 0.5),
          dw_chronic    = quantile(lhs_sample_young$dw_chronic, 0.5),

          dur_acute     = quantile(lhs_sample_young$dur_acute, 0.5),
          dur_subacute  = quantile(lhs_sample_young$dur_subac, 0.5),
          dur_chronic   = quantile(lhs_sample_young$dur_chronic, 0.5),

          dw_subacute   = quantile(lhs_sample_young$dw_subac, 0.5),

          # ---- v3 CHANGE: age-specific vectors (length = n_age) ----
          subac_prop    = .subac_prop_vec,
          chr_prop      = .chr_prop_vec,
          # ---- end v3 change ----

          le_left_vec   = le_left_vec,
          rho_vec       = rho_vec
        )
        
        dplyr::mutate(draw_totals,
                      Region   = region_name,
                      VE       = ve_name,
                      Coverage = cov_name,
                      Scenario = scen_id)
      })
    })
  })
})
step_done("Step 6.2 all_draws_daly_true")
step_done("SECTION 06")

# =============================================================================
# SECTION 07. SAE draws and persistence of core draw objects
# =============================================================================
step_msg("SECTION 07", "SAE draws and saving core draw tables")
# Step 7.1) SAE = hospitalisation + fatal events
all_draws_sae_true <- all_draws_hosp_true %>%
  dplyr::rename(
    total_pre_hosp  = total_pre,
    total_post_hosp = total_post
  ) %>%
  inner_join(
    all_draws_fatal_true %>%
      dplyr::rename(
        total_pre_fatal  = total_pre,
        total_post_fatal = total_post
      ),
    by = c("draw_id", "Region", "VE", "Coverage", "Scenario")
  ) %>%
  mutate(
    total_pre  = total_pre_hosp  + total_pre_fatal,
    total_post = total_post_hosp + total_post_fatal
  ) %>%
  dplyr::select(draw_id, total_pre, total_post, Region, VE, Coverage, Scenario)



# Step 7.2) Save core draw-level tables
step_msg("Step 7.2", "saving all_draws_* to 01_Data/")
save(all_draws_ix_true, file = paste0("01_Data/all_draws_ix_true", MODE_SUFFIX, ".RData"))
save(all_draws_hosp_true, file = paste0("01_Data/all_draws_hosp_true", MODE_SUFFIX, ".RData"))
save(all_draws_daly_true, file = paste0("01_Data/all_draws_daly_true", MODE_SUFFIX, ".RData"))
save(all_draws_fatal_true, file = paste0("01_Data/all_draws_fatal_true", MODE_SUFFIX, ".RData"))
save(all_draws_sae_true, file = paste0("01_Data/all_draws_sae_true", MODE_SUFFIX, ".RData"))
step_done("Step 7.2 save all_draws_*")
step_done("SECTION 07")

# =============================================================================
# SECTION 08. BRR pipeline (draw-level benefit vs risk; no object mixing)
# Everything is created as *_true and never overwritten.
# =============================================================================
step_msg("SECTION 08", "BRR pipeline (benefit vs risk)")

# -----------------------------------------------------------------------------
# SECTION 08A. Helper functions and lookup maps
# -----------------------------------------------------------------------------
step_msg("SECTION 08A", "Helper functions and lookup maps")

# -------------------------------
# Step 8.0) Optional safety cleanup for conflicting objects
# -------------------------------
# rm(list = ls(pattern = "^(benefit_draws|risk_draws|joint|draw_level_xy|summary_long|summary_df|brr_)"))
# rm(list = ls(pattern = "^(tot_vacc_map|tot_vacc_map2)$"))

# -------------------------------
# Step 8.1) Helper: scenario -> AgeCat mapping (integer scenarios)
# -------------------------------
map_scenario_agecat_int <- function(scen_int) {
  dplyr::case_when(
    scen_int == 1L ~ "1-11",
    scen_int == 2L ~ "12-17",
    scen_int == 3L ~ "18-64",
    scen_int == 4L ~ "65+",
    TRUE ~ NA_character_
  )
}

# -------------------------------
# Step 8.2) Build tot_vacc_map_true from national vaccine data (cov50 only)
# Output keys: Scenario (int), VE, AgeCat, tot_vacc_grp
# -------------------------------
tot_vacc_map_true <- combined_nnv_df_region_coverage_model %>%
  filter(VC == "cov50") %>%
  mutate(
    AgeCat = case_when(
      AgeGroup %in% 2:4   ~ "1-11",
      AgeGroup == 5       ~ "12-17",
      AgeGroup %in% 6:15  ~ "18-64",
      AgeGroup %in% 16:20 ~ "65+",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(AgeCat)) %>%
  group_by(scenario, region, AgeCat, VE) %>%
  summarise(tot_vacc_grp = sum(tot_vacc, na.rm = TRUE), .groups = "drop") %>%
  # keep only the targeted AgeCat for each scenario
  mutate(
    target = case_when(
      scenario == "Scenario_1" & AgeCat == "1-11"  ~ 1L,
      scenario == "Scenario_2" & AgeCat == "12-17" ~ 1L,
      scenario == "Scenario_3" & AgeCat == "18-64" ~ 1L,
      scenario == "Scenario_4" & AgeCat == "65+"   ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  dplyr::filter(target == 1L) %>%
  transmute(
    scenario = as.integer(gsub("Scenario_", "", scenario)),
    AgeCat, VE, tot_vacc_grp, region
  )%>%
  dplyr::rename(Region   = region,
                Scenario = scenario)

# -------------------------------
# Step 8.3) Build DALY parameter draws and benefit draws
# -------------------------------
make_averted_draws_true <- function(df_true, outcome_name, coverage_keep = "cov50") {
  df_true %>%
    dplyr::filter(Coverage == coverage_keep) %>%
    dplyr::mutate(
      Scenario = as.integer(Scenario),
      outcome  = outcome_name,
      baseline = total_pre,
      post     = total_post,
      averted  = total_pre - total_post
    ) %>%
    dplyr::select(
      draw_id, Region, VE, Coverage, Scenario, outcome,
      baseline, post, averted
    )
}

# Stack averted burden by outcome (cov50 only); same logic as brazil_all_draws_ori.R
step_msg("Step 8.3", "benefit_draws_true and daly_pars_true")
benefit_draws_true <- dplyr::bind_rows(
  make_averted_draws_true(all_draws_daly_true,  "DALY"),
  make_averted_draws_true(all_draws_sae_true,   "SAE"),
  make_averted_draws_true(all_draws_fatal_true, "Death")
)

# DALY disability weights / durations etc., keyed by draw_id for compute_daly_one()
daly_pars_true <- lhs_sample %>%
  dplyr::mutate(draw_id = dplyr::row_number()) %>%
  dplyr::select(
    draw_id,
    le_lost_1_11, le_lost_12_17, le_lost_18_64, le_lost_65,
    dw_hosp, dw_nonhosp, dw_subac, dw_chronic,
    dur_acute, dur_nonhosp, dur_subac, dur_6m, dur_12m, dur_30m,
    acute, subac, chr6m, chr12m, chr30m
  )
step_done("Step 8.3 benefit_draws_true")

# -------------------------------
# Step 8.3.a) Reporting/plot helper functions
# -------------------------------
make_pr_gt1_wide <- function(ceac_ob, brr_type_filter) {
  ceac_ob %>%
    filter(
      abs(log10(threshold)) < 1e-12,
      brr_type == brr_type_filter,
      RR_seropos == 0
    ) %>%
    transmute(
      Outcome     = outcome,
      Setting     = setting,
      `Age group` = AgeCat,
      VE_col      = VE_label,
      pr_fmt      = sprintf("%.0f%%", 100 * p_accept)
    ) %>%
    pivot_wider(
      names_from  = VE_col,
      values_from = pr_fmt,
      names_glue  = "{VE_col} Pr(BRR>1)"
    )
}

# Format median and uncertainty interval as "med (lo-hi)"
fmt_ci <- function(med, lo, hi) {
  sprintf("%.1f (%.1f–%.1f)", med, lo, hi)
}


## BRR probability curve and data generation function
plot_brr_ceac_outbreak_ve <- function(ceac_df,
                                      target_outcome = "DALY",
                                      setting_levels = c("Low", "Moderate", "High"),
                                      age_levels     = c("18-64", "65+"),
                                      ve_levels      = NULL,
                                      brr_type_labels = c("brr_base" = "Base (no serostatus)",
                                                          "brr_adj"  = "Adjusted (serostatus)")) {
  df <- ceac_df %>%
    filter(
      outcome    == target_outcome,
      RR_seropos == 0
    ) %>%
    mutate(
      setting   = factor(setting,   levels = setting_levels),
      AgeCat    = factor(AgeCat,    levels = age_levels),
      brr_type  = factor(brr_type,  levels = names(brr_type_labels),
                         labels = brr_type_labels)
    )
  
  if (!is.null(ve_levels)) {
    df <- df %>% mutate(VE_label = factor(VE_label, levels = ve_levels))
  }
  
  ggplot(df, aes(x = threshold, y = p_accept,
                 colour   = AgeCat,
                 linetype = brr_type)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.4) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_linetype_manual(values = c("Base (no serostatus)"  = "solid",
                                     "Adjusted (serostatus)" = "dashed")) +
    facet_grid(setting ~ VE_label) +
    labs(
      x        = "Benefit-risk ratio (BRR)",
      y        = "Probability (BRR > t)",
      title    = paste0("BRR acceptability curve: ", target_outcome),
      colour   = "Age group",
      linetype = "BRR type"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank())
}


make_brr_ceac_outbreak <- function(brr_long,
                                   thresholds = 10^seq(-1, 1, by = 0.02),
                                   group_vars = c("setting","VE_label","AgeCat","outcome")) {
  brr_long %>%
    tidyr::crossing(threshold = thresholds) %>%
    group_by(across(all_of(group_vars)), threshold) %>%
    summarise(
      p_accept = mean(brr > threshold),
      .groups = "drop"
    )
}

step_done("SECTION 08A")

# -----------------------------------------------------------------------------
# SECTION 08B. Core BRR computation tables
# -----------------------------------------------------------------------------
step_msg("SECTION 08B", "Core BRR computation tables")

# -------------------------------
# Step 8.6) Build benefit base
# -------------------------------
step_msg("Step 8.6", "building benefit_base_true")
benefit_base_true <- benefit_draws_true %>%
  dplyr::mutate(
    Scenario = as.integer(Scenario),
    AgeCat   = map_scenario_agecat_int(Scenario),
    risk_band = dplyr::if_else(AgeCat == "65+", "65+", "u65")
  ) %>%
  dplyr::left_join(daly_pars_true, by = "draw_id")
step_done("Step 8.6 benefit_base_true")

### Step 8.10) Added on 2026-04-04: serostatus-specific risk decomposition -----
step_msg("Step 8.10", "serostatus-specific risk decomposition")
age_map <- data.frame(
  age_index = 1:20,
  age_gr = c(
    "<1", "1-4", "5-9", "10-11", "12-17", "18-19", "20-24", "25-29",
    "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
    "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
  ),
  AgeCat = c(
    "1-11", "1-11", "1-11", "1-11",
    "12-17",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "65+", "65+", "65+", "65+", "65+"
  )
)

### Step 8.10.1)
### Baseline harms for 65+ reproduce the original implementation after aligning
### the originally assigned risk draw.
### Under-65 baseline harms intentionally differ because vaccine-attributable
### death risk is fixed at zero in this revised implementation.
calc_q_seromix_for_scenarios_agecat <- function(sim_region_ve_cov, age_map) {
  n_scenarios <- length(sim_region_ve_cov)
  out <- vector("list", n_scenarios)
  
  for (sc in seq_len(n_scenarios)) {
    raw_alloc_array <- sim_region_ve_cov[[sc]]$sim_result$raw_allocation_array
    vacc_to_S_array <- sim_region_ve_cov[[sc]]$sim_result$vacc_to_S_array
    
    stopifnot(all(dim(raw_alloc_array) == dim(vacc_to_S_array)))
    
    n_draws <- dim(raw_alloc_array)[3]
    sc_df <- vector("list", n_draws)
    
    for (d in seq_len(n_draws)) {
      df <- data.frame(
        age_index = seq_len(dim(raw_alloc_array)[1]),
        total_vacc_age = rowSums(raw_alloc_array[, , d, drop = FALSE], na.rm = TRUE),
        seroneg_vacc_age = rowSums(vacc_to_S_array[, , d, drop = FALSE], na.rm = TRUE)
      ) %>%
        dplyr::mutate(
          seropos_vacc_age = pmax(0, total_vacc_age - seroneg_vacc_age)
        ) %>%
        dplyr::left_join(age_map, by = "age_index") %>%
        dplyr::group_by(AgeCat) %>%
        dplyr::summarise(
          total_vacc_age = sum(total_vacc_age, na.rm = TRUE),
          seroneg_vacc_age = sum(seroneg_vacc_age, na.rm = TRUE),
          seropos_vacc_age = sum(seropos_vacc_age, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          Scenario = sc,
          draw_id = d,
          q_seroneg_vacc = ifelse(total_vacc_age > 0, seroneg_vacc_age / total_vacc_age, NA_real_),
          q_seropos_vacc = ifelse(total_vacc_age > 0, seropos_vacc_age / total_vacc_age, NA_real_)
        )
      
      sc_df[[d]] <- df
    }
    
    out[[sc]] <- dplyr::bind_rows(sc_df)
  }
  
  dplyr::bind_rows(out)
}

step_msg("Step 8.10.1", "calc_q_seromix — q_all_regions (slow — per region)")
q_all_regions <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
  message(sprintf("[03_v3 %s]   ... q_all_regions region: %s", format(Sys.time(), "%H:%M:%S"), region_name))
  purrr::imap_dfr(region_list, function(ve_list, ve_name) {
    
    purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
      
      q_tmp <- calc_q_seromix_for_scenarios_agecat(
        sim_region_ve_cov = cov_list$scenario_result,
        age_map = age_map
      )
      
      q_tmp %>%
        dplyr::mutate(
          Region = region_name,
          VE = ve_name,
          Coverage = cov_name
        ) %>%
        dplyr::select(
          Region, VE, Coverage, Scenario, draw_id, AgeCat,
          total_vacc_age, seroneg_vacc_age, seropos_vacc_age,
          q_seroneg_vacc, q_seropos_vacc
        )
    })
  })
})
step_done("Step 8.10.1 q_all_regions")

### Step 8.10.2) Generate risk table from LHS sample
step_msg("Step 8.10.2", "risk_draw_df from lhs_sample")
risk_draw_df <- bind_rows(
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "1-11",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_u65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]
  ),
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "12-17",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_u65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]
  ),
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "18-64",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_u65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]
  ),
  tibble(
    draw_id = 1:nrow(lhs_sample),
    AgeCat = "65+",
    p_sae_vacc_base   = lhs_sample[, "p_sae_vacc_65"],
    p_death_vacc_base = lhs_sample[, "p_death_vacc_65"]
  )
)
step_done("Step 8.10.2 risk_draw_df")

## Step 8.10.3) Combine risk and regional composition data
step_msg("Step 8.10.3", "risk_components_all")
rr_vals <- c(1.0, 0.5, 0.1, 0.0)

risk_components_all <- q_all_regions %>%
  tidyr::crossing(RR_seropos = rr_vals) %>%
  dplyr::left_join(
    risk_draw_df,
    by = c("draw_id", "AgeCat")
  )

risk_components_all <- risk_components_all %>%
  dplyr::mutate(
    # conditional risk by serostatus
    p_sae_vacc_seroneg   = p_sae_vacc_base,
    p_death_vacc_seroneg = p_death_vacc_base,
    
    p_sae_vacc_seropos   = p_sae_vacc_base * RR_seropos,
    p_death_vacc_seropos = p_death_vacc_base * RR_seropos,
    
    # contribution to vaccinated cohort
    p_sae_vacc_seroneg_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seroneg_vacc * p_sae_vacc_seroneg,
      0
    ),
    p_death_vacc_seroneg_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seroneg_vacc * p_death_vacc_seroneg,
      0
    ),
    
    p_sae_vacc_seropos_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seropos_vacc * p_sae_vacc_seropos,
      0
    ),
    p_death_vacc_seropos_contrib = dplyr::if_else(
      total_vacc_age > 0,
      q_seropos_vacc * p_death_vacc_seropos,
      0
    ),
    
    # total adjusted
    p_sae_vacc_adj = p_sae_vacc_seroneg_contrib + p_sae_vacc_seropos_contrib,
    p_death_vacc_adj = p_death_vacc_seroneg_contrib + p_death_vacc_seropos_contrib
  )
step_done("Step 8.10.3 risk_components_all")

# Step 8.10.4) Estimate SAE and death per 10,000 vaccinated
risk_components_all <- risk_components_all %>%
  dplyr::mutate(
    # baseline
    sae_10k_base   = 1e4 * p_sae_vacc_base + 1e4 * p_death_vacc_base,
    death_10k_base = 1e4 * p_death_vacc_base,
    
    # seronegative contribution
    sae_10k_seroneg   = 1e4 * p_sae_vacc_seroneg_contrib + 1e4 * p_death_vacc_seroneg_contrib,
    death_10k_seroneg = 1e4 * p_death_vacc_seroneg_contrib,
    
    # seropositive contribution
    sae_10k_seropos   = 1e4 * p_sae_vacc_seropos_contrib + 1e4 * p_death_vacc_seropos_contrib,
    death_10k_seropos = 1e4 * p_death_vacc_seropos_contrib,
    
    # adjusted total
    sae_10k_adj   = sae_10k_seroneg + sae_10k_seropos,
    death_10k_adj = death_10k_seroneg + death_10k_seropos
  )

# Step 8.10.5) Join serostatus risk decomposition to benefit data
risk_join_all <- risk_components_all %>%
  dplyr::select(
    Region, VE, Coverage, Scenario, draw_id, AgeCat, RR_seropos,
    total_vacc_age, q_seroneg_vacc, q_seropos_vacc,
    p_sae_vacc_base, p_death_vacc_base,
    p_sae_vacc_seroneg, p_death_vacc_seroneg,
    p_sae_vacc_seropos, p_death_vacc_seropos,
    p_sae_vacc_seroneg_contrib, p_death_vacc_seroneg_contrib,
    p_sae_vacc_seropos_contrib, p_death_vacc_seropos_contrib,
    p_sae_vacc_adj, p_death_vacc_adj,
    sae_10k_base, death_10k_base,
    sae_10k_seroneg, death_10k_seroneg,
    sae_10k_seropos, death_10k_seropos,
    sae_10k_adj, death_10k_adj
  )

# Step 8.10.6) Build draw-level burden table
step_msg("Step 8.10.6", "benefit_draw_df and draw_level_xy join")
benefit_draw_df <- benefit_base_true %>%
  dplyr::left_join(
    tot_vacc_map_true,
    by = c("Region", "Scenario", "VE", "AgeCat")
  ) %>%
  dplyr::mutate(
    Scenario = as.integer(Scenario),
    baseline_10k = (baseline / tot_vacc_grp) * 1e4,
    post_10k     = (post / tot_vacc_grp) * 1e4,
    averted_10k  = (averted / tot_vacc_grp) * 1e4,
    pct_reduction = dplyr::if_else(
      baseline > 0,
      100 * averted / baseline,
      NA_real_
    )
  )

benefit_draw_df <- tidyr::crossing(
  benefit_draw_df,
  RR_seropos = rr_vals
)

# join benefit + risk decomposition
draw_level_xy_serostatus <- benefit_draw_df %>%
  dplyr::left_join(
    risk_join_all,
    by = c("Region", "VE", "Coverage", "Scenario", "draw_id", "AgeCat", "RR_seropos")
  )

# Step 8.10.7) Compute vaccine-caused DALY components
step_msg("Step 8.10.7", "rowwise DALY computation (very slow)")
draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    daly_10k_base = if (outcome == "DALY") {
      compute_daly_one_age_specific(
        age_group      = AgeCat,
        sae_10k        = sae_10k_base,
        deaths_sae_10k = death_10k_base,
        draw_id        = draw_id,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    daly_10k_seroneg = if (outcome == "DALY") {
      compute_daly_one_age_specific(
        age_group      = AgeCat,
        sae_10k        = sae_10k_seroneg,
        deaths_sae_10k = death_10k_seroneg,
        draw_id        = draw_id,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    daly_10k_seropos = if (outcome == "DALY") {
      compute_daly_one_age_specific(
        age_group      = AgeCat,
        sae_10k        = sae_10k_seropos,
        deaths_sae_10k = death_10k_seropos,
        draw_id        = draw_id,
        draw_pars = list(
          le_lost_1_11  = le_lost_1_11,
          le_lost_12_17 = le_lost_12_17,
          le_lost_18_64 = le_lost_18_64,
          le_lost_65    = le_lost_65,
          dw_hosp       = dw_hosp,
          dw_nonhosp    = dw_nonhosp,
          dw_subac      = dw_subac,
          dw_chronic    = dw_chronic,
          dur_acute     = dur_acute,
          dur_nonhosp   = dur_nonhosp,
          dur_subac     = dur_subac,
          dur_6m        = dur_6m,
          dur_12m       = dur_12m,
          dur_30m       = dur_30m,
          acute         = acute,
          subac         = subac,
          chr6m         = chr6m,
          chr12m        = chr12m,
          chr30m        = chr30m
        )
      )$daly_sae
    } else {
      NA_real_
    },
    
    daly_10k_adj = if (outcome == "DALY") {
      daly_10k_seroneg + daly_10k_seropos
    } else {
      NA_real_
    }
  ) %>%
  dplyr::ungroup()
step_done("Step 8.10.7 rowwise DALY")

# Step 8.10.8) Add reporting labels
step_msg("Step 8.10.8–8.10.12", "BRR summaries and wide tables")
draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
  dplyr::mutate(
    VE_label = factor(
      VE,
      levels = c("VE0", "VE98.9"),
      labels = c("Disease blocking only", "Disease and infection blocking")
    ),
    setting = unname(setting_key[Region])
  )

draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
  dplyr::mutate(
    x_10k_base = dplyr::case_when(
      outcome == "SAE"   ~ sae_10k_base,
      outcome == "Death" ~ death_10k_base,
      outcome == "DALY"  ~ daly_10k_base,
      TRUE ~ NA_real_
    ),
    x_10k_adj = dplyr::case_when(
      outcome == "SAE"   ~ sae_10k_adj,
      outcome == "Death" ~ death_10k_adj,
      outcome == "DALY"  ~ daly_10k_adj,
      TRUE ~ NA_real_
    ),
    brr_base = dplyr::if_else(
      is.na(x_10k_base) | x_10k_base == 0,
      NA_real_,
      averted_10k / x_10k_base
    ),
    brr_adj = dplyr::if_else(
      is.na(x_10k_adj) | x_10k_adj == 0,
      NA_real_,
      averted_10k / x_10k_adj
    )
  )

## Step 8.10.9) Summary across all age groups (including 1-17)
brr_draw_summary_true <- draw_level_xy_serostatus %>% 
  dplyr::mutate(
    brr_base = ifelse(is.infinite(brr_base), NA, brr_base),
    brr_adj  = ifelse(is.infinite(brr_adj),  NA, brr_adj),
    setting  = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  dplyr::group_by(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting) %>%  
  dplyr::summarise(
    # BRR
    brr_base_med = quantile(brr_base, 0.50,  na.rm = TRUE),
    brr_base_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
    brr_base_hi  = quantile(brr_base, 0.975, na.rm = TRUE),
    
    brr_adj_med  = quantile(brr_adj,  0.50,  na.rm = TRUE),
    brr_adj_lo   = quantile(brr_adj,  0.025, na.rm = TRUE),
    brr_adj_hi   = quantile(brr_adj,  0.975, na.rm = TRUE),
    
    # Averted burden (outcome-specific; distinguished by outcome column)
    av_med = quantile(averted_10k, 0.50,  na.rm = TRUE),
    av_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
    av_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),
    
    # Caused — base
    ca_sae_base_med   = quantile(sae_10k_base,   0.50,  na.rm = TRUE),
    ca_sae_base_lo    = quantile(sae_10k_base,   0.025, na.rm = TRUE),
    ca_sae_base_hi    = quantile(sae_10k_base,   0.975, na.rm = TRUE),
    
    ca_death_base_med = quantile(death_10k_base, 0.50,  na.rm = TRUE),
    ca_death_base_lo  = quantile(death_10k_base, 0.025, na.rm = TRUE),
    ca_death_base_hi  = quantile(death_10k_base, 0.975, na.rm = TRUE),
    
    ca_daly_base_med  = quantile(daly_10k_base,  0.50,  na.rm = TRUE),
    ca_daly_base_lo   = quantile(daly_10k_base,  0.025, na.rm = TRUE),
    ca_daly_base_hi   = quantile(daly_10k_base,  0.975, na.rm = TRUE),
    
    # Caused — adj (serostatus-adjusted)
    ca_sae_adj_med    = quantile(sae_10k_adj,    0.50,  na.rm = TRUE),
    ca_sae_adj_lo     = quantile(sae_10k_adj,    0.025, na.rm = TRUE),
    ca_sae_adj_hi     = quantile(sae_10k_adj,    0.975, na.rm = TRUE),
    
    ca_death_adj_med  = quantile(death_10k_adj,  0.50,  na.rm = TRUE),
    ca_death_adj_lo   = quantile(death_10k_adj,  0.025, na.rm = TRUE),
    ca_death_adj_hi   = quantile(death_10k_adj,  0.975, na.rm = TRUE),
    
    ca_daly_adj_med   = quantile(daly_10k_adj,   0.50,  na.rm = TRUE),
    ca_daly_adj_lo    = quantile(daly_10k_adj,   0.025, na.rm = TRUE),
    ca_daly_adj_hi    = quantile(daly_10k_adj,   0.975, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    scenario  = Scenario,
    age_group = AgeCat
  ) %>%
  dplyr::arrange(outcome, Scenario, setting, AgeCat, VE_label, RR_seropos)

## Step 8.10.10) Filter to >=18 years and RR_seropos = 0
brr_draw_summary_true_filtered <- brr_draw_summary_true %>%
  dplyr::filter(
    AgeCat %in% c("18-64", "65+"),
    RR_seropos == 0.0
  ) %>%
  dplyr::select(-scenario, -age_group)  


## Step 8.10.11) Build wide-format reporting tables
# Step 8.10.12) Merge caused-burden columns into wide table
ca_summary <- brr_draw_summary_true_filtered %>%
  filter(outcome == "DALY") %>%  
  select(Scenario, AgeCat, VE_label, RR_seropos, setting,
         starts_with("ca_"))

brr_table_long_serostatus <- brr_draw_summary_true_filtered %>%
  select(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting,
         brr_base_med, brr_base_lo, brr_base_hi,
         brr_adj_med,  brr_adj_lo,  brr_adj_hi,
         av_med, av_lo, av_hi) %>%
  left_join(ca_summary,
            by = c("Scenario", "AgeCat", "VE_label", "RR_seropos", "setting")) %>%
  mutate(
    brr_base_formatted = sprintf("%.2f (%.2f–%.2f)", brr_base_med, brr_base_lo, brr_base_hi),
    brr_adj_formatted  = sprintf("%.2f (%.2f–%.2f)", brr_adj_med,  brr_adj_lo,  brr_adj_hi),
    av_formatted       = sprintf("%.2f (%.2f–%.2f)", av_med,       av_lo,       av_hi),
    ca_sae_base_formatted   = sprintf("%.2f (%.2f–%.2f)", ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
    ca_death_base_formatted = sprintf("%.2f (%.2f–%.2f)", ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
    ca_daly_base_formatted  = sprintf("%.2f (%.2f–%.2f)", ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
    ca_sae_adj_formatted    = sprintf("%.2f (%.2f–%.2f)", ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
    ca_death_adj_formatted  = sprintf("%.2f (%.2f–%.2f)", ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
    ca_daly_adj_formatted   = sprintf("%.2f (%.2f–%.2f)", ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
  )


brr_long_serostatus <- draw_level_xy_serostatus %>%
  mutate(
    outcome = recode(outcome, "sae" = "SAE", "death" = "Death", "daly" = "DALY"),
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    AgeCat  = factor(AgeCat,  levels = c("18-64", "65+"))
  ) %>%
  filter(
    is.finite(brr_base), brr_base > 0,
    is.finite(brr_adj),  brr_adj  > 0
  ) %>%
  transmute(
    Region, setting, VE_label, AgeCat, outcome, RR_seropos,
    brr_base,
    brr_adj
  )

brr_long_serostatus_long <- brr_long_serostatus %>%
  filter(!is.na(AgeCat)) %>%
  pivot_longer(
    cols      = c(brr_base, brr_adj),
    names_to  = "brr_type",
    values_to = "brr"
  ) %>%
  filter(is.finite(brr), brr > 0)

brr_range <- brr_long_serostatus_long %>%
  filter(RR_seropos == 0) %>%
  summarise(
    min_brr = min(brr, na.rm = TRUE),
    max_brr = max(brr, na.rm = TRUE)
  )

brr_min <- brr_range$min_brr
brr_max <- brr_range$max_brr

thresholds_auto <- 10^seq(
  floor(log10(brr_min)),   
  ceiling(log10(brr_max)), 
  by = 0.02
)
ceac_ob <- make_brr_ceac_outbreak(
  brr_long_serostatus_long,
  group_vars = c("setting", "VE_label", "AgeCat", "outcome", "RR_seropos", "brr_type")
)

pr_gt1_wide_setting_base <- make_pr_gt1_wide(ceac_ob, "brr_base")
pr_gt1_wide_setting_adj  <- make_pr_gt1_wide(ceac_ob, "brr_adj")

ca_summary <- brr_draw_summary_true_filtered %>%
  filter(RR_seropos == 0, outcome == "DALY") %>%
  select(Scenario, AgeCat, VE_label, setting,
         ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi,
         ca_death_base_med, ca_death_base_lo, ca_death_base_hi,
         ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi,
         ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi,
         ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi,
         ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)

brr_draw_for_wide <- brr_draw_summary_true_filtered %>%
  filter(RR_seropos == 0) %>%
  select(outcome, Scenario, AgeCat, VE_label, setting,
         av_med,       av_lo,       av_hi,
         brr_base_med, brr_base_lo, brr_base_hi,
         brr_adj_med,  brr_adj_lo,  brr_adj_hi) %>%
  left_join(ca_summary, by = c("Scenario", "AgeCat", "VE_label", "setting"))

brr_table_wide_serostatus <- brr_draw_for_wide %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    
    av_formatted            = fmt_ci(av_med,       av_lo,       av_hi),
    brr_base_formatted      = fmt_ci(brr_base_med, brr_base_lo, brr_base_hi),
    brr_adj_formatted       = fmt_ci(brr_adj_med,  brr_adj_lo,  brr_adj_hi),
    ca_sae_base_formatted   = fmt_ci(ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
    ca_death_base_formatted = fmt_ci(ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
    ca_daly_base_formatted  = fmt_ci(ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
    ca_sae_adj_formatted    = fmt_ci(ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
    ca_death_adj_formatted  = fmt_ci(ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
    ca_daly_adj_formatted   = fmt_ci(ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
  ) %>%
  select(outcome, Scenario, AgeCat, setting, VE_label,
         av_formatted, brr_base_formatted, brr_adj_formatted,
         ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
         ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted) %>%
  pivot_wider(
    names_from  = VE_label,
    values_from = c(av_formatted, brr_base_formatted, brr_adj_formatted,
                    ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
                    ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted),
    names_glue  = "{VE_label}_{.value}"
  ) %>%
  arrange(outcome, Scenario, setting, AgeCat) %>%
  dplyr::rename(
    Outcome     = outcome,
    `Age group` = AgeCat,
    Setting     = setting
  )

brr_table_wide_serostatus <- brr_table_wide_serostatus %>%
  dplyr::rename(
    # Disease blocking only
    `DB (Averted)`           = `Disease blocking only_av_formatted`,
    `DB (BRR base)`          = `Disease blocking only_brr_base_formatted`,
    `DB (BRR adj)`           = `Disease blocking only_brr_adj_formatted`,
    `DB (SAE caused base)`   = `Disease blocking only_ca_sae_base_formatted`,
    `DB (SAE caused adj)`    = `Disease blocking only_ca_sae_adj_formatted`,
    `DB (Death caused base)` = `Disease blocking only_ca_death_base_formatted`,
    `DB (Death caused adj)`  = `Disease blocking only_ca_death_adj_formatted`,
    `DB (DALY caused base)`  = `Disease blocking only_ca_daly_base_formatted`,
    `DB (DALY caused adj)`   = `Disease blocking only_ca_daly_adj_formatted`,
    # Disease and infection blocking
    `DIB (Averted)`           = `Disease and infection blocking_av_formatted`,
    `DIB (BRR base)`          = `Disease and infection blocking_brr_base_formatted`,
    `DIB (BRR adj)`           = `Disease and infection blocking_brr_adj_formatted`,
    `DIB (SAE caused base)`   = `Disease and infection blocking_ca_sae_base_formatted`,
    `DIB (SAE caused adj)`    = `Disease and infection blocking_ca_sae_adj_formatted`,
    `DIB (Death caused base)` = `Disease and infection blocking_ca_death_base_formatted`,
    `DIB (Death caused adj)`  = `Disease and infection blocking_ca_death_adj_formatted`,
    `DIB (DALY caused base)`  = `Disease and infection blocking_ca_daly_base_formatted`,
    `DIB (DALY caused adj)`   = `Disease and infection blocking_ca_daly_adj_formatted`
  )

pr_gt1_wide_serostatus_base <- pr_gt1_wide_setting_base %>%
  dplyr::rename_with(
    ~ gsub("Pr(BRR>1)", "Pr(BRR_base>1)", .x, fixed = TRUE),
    .cols = -c(Outcome, Setting, `Age group`)
  )

pr_gt1_wide_serostatus_adj <- pr_gt1_wide_setting_adj %>%
  dplyr::rename_with(
    ~ gsub("Pr(BRR>1)", "Pr(BRR_adj>1)", .x, fixed = TRUE),
    .cols = -c(Outcome, Setting, `Age group`)
  )

brr_table_wide_serostatus2 <- brr_table_wide_serostatus %>%
  left_join(pr_gt1_wide_serostatus_base, by = c("Outcome", "Setting", "Age group")) %>%
  left_join(pr_gt1_wide_serostatus_adj,  by = c("Outcome", "Setting", "Age group"))

brr_table_wide_serostatus2 <- brr_table_wide_serostatus2 %>%
  relocate(`DB (Averted)`,                                  .after = `Age group`) %>%
  relocate(`DB (BRR base)`,                                 .after = `DB (Averted)`) %>%
  relocate(`DB (BRR adj)`,                                  .after = `DB (BRR base)`) %>%
  relocate(`Disease blocking only Pr(BRR_base>1)`,          .after = `DB (BRR adj)`) %>%
  relocate(`Disease blocking only Pr(BRR_adj>1)`,           .after = `Disease blocking only Pr(BRR_base>1)`) %>%
  relocate(`DB (SAE caused base)`,                          .after = `Disease blocking only Pr(BRR_adj>1)`) %>%
  relocate(`DB (SAE caused adj)`,                           .after = `DB (SAE caused base)`) %>%
  relocate(`DB (Death caused base)`,                        .after = `DB (SAE caused adj)`) %>%
  relocate(`DB (Death caused adj)`,                         .after = `DB (Death caused base)`) %>%
  relocate(`DB (DALY caused base)`,                         .after = `DB (Death caused adj)`) %>%
  relocate(`DB (DALY caused adj)`,                          .after = `DB (DALY caused base)`) %>%
  relocate(`DIB (Averted)`,                                 .after = `DB (DALY caused adj)`) %>%
  relocate(`DIB (BRR base)`,                                .after = `DIB (Averted)`) %>%
  relocate(`DIB (BRR adj)`,                                 .after = `DIB (BRR base)`) %>%
  relocate(`Disease and infection blocking Pr(BRR_base>1)`, .after = `DIB (BRR adj)`) %>%
  relocate(`Disease and infection blocking Pr(BRR_adj>1)`,  .after = `Disease and infection blocking Pr(BRR_base>1)`) %>%
  relocate(`DIB (SAE caused base)`,                         .after = `Disease and infection blocking Pr(BRR_adj>1)`) %>%
  relocate(`DIB (SAE caused adj)`,                          .after = `DIB (SAE caused base)`) %>%
  relocate(`DIB (Death caused base)`,                       .after = `DIB (SAE caused adj)`) %>%
  relocate(`DIB (Death caused adj)`,                        .after = `DIB (Death caused base)`) %>%
  relocate(`DIB (DALY caused base)`,                        .after = `DIB (Death caused adj)`) %>%
  relocate(`DIB (DALY caused adj)`,                         .after = `DIB (DALY caused base)`) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ tidyr::replace_na(as.character(.x), "beneficial")
    )
  )

step_done("Step 8.10.8–8.10.12 BRR summaries")
step_done("SECTION 08B")

# -----------------------------------------------------------------------------
# SECTION 08C. Reporting table assembly and document export
# -----------------------------------------------------------------------------
step_msg("SECTION 08C", "Export BRR table to Word")

# Step 8.10.13) Build publication-ready BRR table and export to Word
# 1) Extract and rename columns for "Disease blocking only" (base + adjusted BRR/risk)
part1 <- brr_table_wide_serostatus2 %>%
  dplyr::mutate(
    Risk_base = dplyr::case_when(
      Outcome == "SAE"   ~ `DB (SAE caused base)`,
      Outcome == "Death" ~ `DB (Death caused base)`,
      Outcome == "DALY"  ~ `DB (DALY caused base)`,
      TRUE ~ NA_character_
    ),
    Risk_adj = dplyr::case_when(
      Outcome == "SAE"   ~ `DB (SAE caused adj)`,
      Outcome == "Death" ~ `DB (Death caused adj)`,
      Outcome == "DALY"  ~ `DB (DALY caused adj)`,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `DB (Averted)`,
    Risk_base,
    Risk_adj,
    BRR_base = `DB (BRR base)`,
    BRR_adj  = `DB (BRR adj)`,
    prob_base = `Disease blocking only Pr(BRR_base>1)`,
    prob_adj  = `Disease blocking only Pr(BRR_adj>1)`
  ) %>%
  dplyr::mutate(mechanism = "Disease blocking only")

# 2) Extract and rename columns for "Disease and infection blocking"
part2 <- brr_table_wide_serostatus2 %>%
  dplyr::mutate(
    Risk_base = dplyr::case_when(
      Outcome == "SAE"   ~ `DIB (SAE caused base)`,
      Outcome == "Death" ~ `DIB (Death caused base)`,
      Outcome == "DALY"  ~ `DIB (DALY caused base)`,
      TRUE ~ NA_character_
    ),
    Risk_adj = dplyr::case_when(
      Outcome == "SAE"   ~ `DIB (SAE caused adj)`,
      Outcome == "Death" ~ `DIB (Death caused adj)`,
      Outcome == "DALY"  ~ `DIB (DALY caused adj)`,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `DIB (Averted)`,
    Risk_base,
    Risk_adj,
    BRR_base = `DIB (BRR base)`,
    BRR_adj  = `DIB (BRR adj)`,
    prob_base = `Disease and infection blocking Pr(BRR_base>1)`,
    prob_adj  = `Disease and infection blocking Pr(BRR_adj>1)`
  ) %>%
  dplyr::mutate(mechanism = "Disease and infection blocking")

# 3) Stack both mechanisms into one long-format reporting table
brr_table_final_long <- dplyr::bind_rows(part1, part2) %>%
  dplyr::mutate(
    Setting = factor(Setting, levels = c("High", "Moderate", "Low"))
  ) %>%
  dplyr::arrange(Outcome, Setting, `Age group`) %>%
  dplyr::select(
    Outcome, Setting, `Age group`, mechanism,
    Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj, prob_base, prob_adj
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj),
      ~ {
        x <- gsub(" \\(", "\n(", .x, fixed = FALSE)
        dplyr::case_when(
          x %in% c(
            "NA (NA-NA)", "NA (NA–NA)",
            "NA\n(NA-NA)", "NA\n(NA–NA)",
            "NA (NA NA)", "NA\n(NA NA)"
          ) ~ "beneficial",
          TRUE ~ x
        )
      }
    )
  )

# ---- Main table: mechanism pivoted to columns, not rows -------------------
# brr_table_final_long (Outcome x Setting x Age group x mechanism, 36 rows for
# 3 outcomes x 3 settings x 2 ages x 2 mechanisms) is kept below as the full
# supplementary table. For the main manuscript table, mechanism (Disease
# blocking only vs Disease and infection blocking) is pivoted into paired
# columns instead -- the same compression already applied to serostatus
# (base vs adjusted, BRR_base/BRR_adj side by side) -- so rows drop to
# Outcome x Setting x Age group (18 rows) without dropping any information.
# Age group is kept as a row per explicit requirement (not pooled/dropped).
# Pr(BRR>1) dropped from the main table (too dense) -- still in the
# supplementary long table (brr_table_final_long) above for full detail.
metric_cols <- c("Benefit", "Risk_base", "Risk_adj", "BRR_base", "BRR_adj")

brr_table_final_wide <- part1 %>%
  dplyr::select(Outcome, Setting, `Age group`, dplyr::all_of(metric_cols)) %>%
  dplyr::rename_with(~ paste0(.x, "_DB"), dplyr::all_of(metric_cols)) %>%
  dplyr::left_join(
    part2 %>%
      dplyr::select(Outcome, Setting, `Age group`, dplyr::all_of(metric_cols)) %>%
      dplyr::rename_with(~ paste0(.x, "_DIB"), dplyr::all_of(metric_cols)),
    by = c("Outcome", "Setting", "Age group")
  ) %>%
  dplyr::mutate(Setting = factor(Setting, levels = c("High", "Moderate", "Low"))) %>%
  dplyr::arrange(Outcome, Setting, `Age group`) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::ends_with(c("_DB", "_DIB")) & dplyr::matches("Benefit|Risk|BRR"),
      ~ {
        x <- gsub(" \\(", "\n(", .x, fixed = FALSE)
        dplyr::case_when(
          x %in% c("NA (NA-NA)", "NA (NA–NA)", "NA\n(NA-NA)", "NA\n(NA–NA)",
                    "NA (NA NA)", "NA\n(NA NA)") ~ "beneficial",
          TRUE ~ x
        )
      }
    )
  )

metric_labels <- c(
  Benefit   = "Benefit\n(per 10,000)",
  Risk_base = "Risk (base)\n(per 10,000)",
  Risk_adj  = "Risk (adj.)\n(per 10,000)",
  BRR_base  = "BRR (base)",
  BRR_adj   = "BRR (adj.)"
)
header_labels_wide <- c(
  Setting = "Setting", `Age group` = "Age group",
  setNames(metric_labels, paste0(names(metric_labels), "_DB")),
  setNames(metric_labels, paste0(names(metric_labels), "_DIB"))
)

# One table PER outcome (Setting x Age group, 6 rows each) instead of one
# combined 18-row table -- DALY is the main-text table; Death and SAE get
# their own files. Outcome column dropped within each file since it's now
# constant (redundant).
build_outcome_table <- function(outcome_filter, out_file, caption_text) {
  df <- brr_table_final_wide %>%
    dplyr::filter(Outcome == outcome_filter) %>%
    dplyr::select(-Outcome)

  ft <- flextable::flextable(df) %>%
    flextable::set_header_labels(values = as.list(header_labels_wide)) %>%
    flextable::add_header_row(
      values = c("", "", "Disease blocking only", "Disease and infection blocking"),
      colwidths = c(1, 1, length(metric_cols), length(metric_cols))
    ) %>%
    flextable::theme_booktabs() %>%
    flextable::bold(part = "header") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(j = 1:2, align = "left", part = "all") %>%
    flextable::merge_v(j = c("Setting", "Age group")) %>%
    flextable::valign(j = c("Setting", "Age group"), valign = "top") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::autofit()

  doc <- officer::read_docx() %>%
    officer::body_add_par(caption_text, style = "heading 2") %>%
    flextable::body_add_flextable(ft)

  print(doc, target = out_file)
  message("Saved: ", out_file)

  xlsx_file <- sub("\\.docx$", ".xlsx", out_file)
  writexl::write_xlsx(df, xlsx_file)
  message("Saved: ", xlsx_file)
}

# Main text table: DALY only.
build_outcome_table(
  "DALY", paste0("06_Results/BRR_table_ori_setting", MODE_SUFFIX, ".docx"),
  "Benefit-Risk Ratio (BRR) by Setting and Age group -- DALY"
)
# Companion single-outcome tables (Death, SAE) -- own files, not in the main table.
build_outcome_table(
  "Death", paste0("06_Results/BRR_table_ori_setting_Death", MODE_SUFFIX, ".docx"),
  "Benefit-Risk Ratio (BRR) by Setting and Age group -- Death"
)
build_outcome_table(
  "SAE", paste0("06_Results/BRR_table_ori_setting_SAE", MODE_SUFFIX, ".docx"),
  "Benefit-Risk Ratio (BRR) by Setting and Age group -- SAE"
)

# ---- Supplementary table: full long format (mechanism as rows, 36 rows) ---
ft_brr_long <- flextable::flextable(brr_table_final_long) %>%
  flextable::set_header_labels(
    Outcome     = "Outcome",
    Setting     = "Setting",
    `Age group` = "Age group",
    mechanism   = "Vaccine protection\nmechanism",
    Benefit     = "Benefit:\nOutcomes averted\n(per 10,000)",
    Risk_base   = "Risk (base):\nOutcomes attributable\n(per 10,000)",
    Risk_adj    = "Risk (adjusted):\nOutcomes attributable\n(per 10,000)",
    BRR_base    = "BRR (base):\n(Prevented per 1 caused)",
    BRR_adj     = "BRR (adjusted):\n(Prevented per 1 caused)",
    prob_base   = "Probability\n(BRR_base > 1)\n(%)",
    prob_adj    = "Probability\n(BRR_adj > 1)\n(%)"
  ) %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:4, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Outcome", "Setting", "Age group")) %>%
  flextable::valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc_supp <- officer::read_docx() %>%
  officer::body_add_par(
    "Supplementary Table: Benefit-Risk Ratio (BRR) by Outcome, Setting, Age group, and Mechanism (full breakdown)",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft_brr_long)

print(doc_supp, target = paste0("06_Results/BRR_table_ori_setting_supplementary", MODE_SUFFIX, ".docx"))
writexl::write_xlsx(brr_table_final_long, paste0("06_Results/BRR_table_ori_setting_supplementary", MODE_SUFFIX, ".xlsx"))
step_done("SECTION 08C — Word export")

save(brr_table_final_long, file = paste0("01_Data/brr_table_final_long", MODE_SUFFIX, ".RData"))
save(draw_level_xy_serostatus, file = paste0("01_Data/draw_level_xy_serostatus", MODE_SUFFIX, ".RData"))

# -----------------------------------------------------------------------------
# SECTION 08D. BRR probability curve plots
# -----------------------------------------------------------------------------
step_msg("SECTION 08D", "BRR probability curve plots")

## BRR probability curves 

p_daly_ceac_outbreak <- plot_brr_ceac_outbreak_ve(ceac_ob, "DALY") +
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag          = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )


plot_sae   <- plot_brr_ceac_outbreak_ve(ceac_ob, "SAE")+
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag          = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

plot_death <- plot_brr_ceac_outbreak_ve(ceac_ob, "Death")+
  theme(text = element_text(family = "Calibri")) +
  labs(tag = "E") +
  theme(
    plot.tag          = element_text(face = "bold", size = 16),
    plot.tag.position = c(0, 1)
  )

step_done("SECTION 08D")
step_done("SECTION 08")

# =============================================================================
# SECTION 09. DIAGNOSTIC - age-bin share of averted DALY by Setting
# =============================================================================
step_msg("SECTION 09", "Diagnostic — age-bin share of averted DALY")
# Purpose: verify the hypothesis that the v2 -> v3 benefit change is large in
# LOW settings because a bigger fraction of adult averted cases falls into the
# o40 bins (ages 40-64 + 65+), vs HIGH settings where 18-39 (u40) dominates.
#
# For each Region (-> Setting), VE, Coverage, Scenario we compute
#   averted_by_bin[a] = mean_over_draws( sum_weeks(pre_rho[a,,]) -
#                                        sum_weeks(post_rho[a,,]) )
# and then report the share of total averted cases that falls into u40 bins
# (1-10, ages 0-39) vs o40 bins (11-20, ages 40-89). We also break out within
# each AgeCat (1-11 / 12-17 / 18-64 / 65+) so you can see WHERE inside 18-64
# the cases sit (18-39 vs 40-64).
# -----------------------------------------------------------------------------

age_bin_labels <- c(
  "0-1","1-4","5-9","10-11","12-17",
  "18-19","20-24","25-29","30-34","35-39",
  "40-44","45-49","50-54","55-59","60-64",
  "65-69","70-74","75-79","80-84","85-89"
)

# AgeCat mapping consistent with ORI scenario definitions
agecat_by_bin <- c(
  "0-1","1-11","1-11","1-11","12-17",        # bins 1-5
  "18-64","18-64","18-64","18-64","18-64",   # bins 6-10  (18-39, u40)
  "18-64","18-64","18-64","18-64","18-64",   # bins 11-15 (40-64, o40)
  "65+","65+","65+","65+","65+"              # bins 16-20 (o40)
)
u40_flag_by_bin <- c(rep("u40", 10), rep("o40", 10))

# Scenario_X -> targeted AgeCat (for filtering to the vaccinated age group)
scenario_to_agecat <- c(
  Scenario_1 = "1-11",
  Scenario_2 = "12-17",
  Scenario_3 = "18-64",
  Scenario_4 = "65+"
)

# Helper: averted cases per age bin, averaged over draws, for ONE
# region x VE x Coverage x Scenario combination. rho scaling is applied so
# this is on the "true" (under-reported-adjusted) scale that enters DALY.
averted_bin_one <- function(pre_list, post_arr, rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[seq_len(n_draws)]
  post_arr <- post_arr[,, seq_len(n_draws), drop = FALSE]
  rho_vec  <- rho_vec[seq_len(n_draws)]

  # rho scaling (mirror of scale_symp_by_rho_pre / post used in v3)
  pre_true  <- scale_symp_by_rho_pre(pre_list, rho_vec)
  post_true <- scale_symp_by_rho_post(post_arr, rho_vec)

  n_age <- nrow(pre_true[[1]])
  averted_mat <- matrix(0, nrow = n_age, ncol = n_draws)
  for (i in seq_len(n_draws)) {
    averted_mat[, i] <- rowSums(pre_true[[i]], na.rm = TRUE) -
      rowSums(post_true[,, i], na.rm = TRUE)
  }
  rowMeans(averted_mat, na.rm = TRUE)  # length n_age
}

set.seed(1)  # match rho sampling seed used for all_draws_daly_true

diagnostic_averted_bin <- purrr::imap_dfr(
  postsim_vc_ixchiq_model,
  function(region_list, region_name) {
    pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
    rho_pool_region <- as.numeric(posterior_list[[region_name]]$rho)
    rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) &
                                         rho_pool_region > 0]

    purrr::imap_dfr(region_list, function(ve_list, ve_name) {
      purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
        purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
          post_arr <- scen$sim_result$age_array_raw_symp
          n_draws  <- min(length(pre_list), dim(post_arr)[3])
          rho_vec  <- sample(rho_pool_region, size = n_draws, replace = TRUE)

          averted_vec <- averted_bin_one(pre_list, post_arr, rho_vec)

          tibble(
            Region   = region_name,
            Setting  = unname(setting_key[region_name]),
            VE       = ve_name,
            Coverage = cov_name,
            Scenario = scen_id,
            bin      = seq_along(averted_vec),
            age_bin  = age_bin_labels[seq_along(averted_vec)],
            AgeCat_bin = agecat_by_bin[seq_along(averted_vec)],
            u40_flag = u40_flag_by_bin[seq_along(averted_vec)],
            averted  = averted_vec
          )
        })
      })
    })
  }
)

# ---- Q1: In the targeted AgeCat of each Scenario, what fraction of averted
#          cases falls in u40 vs o40 bins, by Setting? -----------------------
share_u40_o40_within_targeted_agecat <- diagnostic_averted_bin %>%
  dplyr::mutate(target_AgeCat = unname(scenario_to_agecat[Scenario])) %>%
  dplyr::filter(AgeCat_bin == target_AgeCat) %>%   # only bins inside vaccinated AgeCat
  dplyr::group_by(Setting, Scenario, AgeCat_bin, VE, Coverage, u40_flag) %>%
  dplyr::summarise(averted = sum(averted, na.rm = TRUE), .groups = "drop_last") %>%
  dplyr::mutate(share = averted / sum(averted, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Scenario, Setting, VE, Coverage, u40_flag)

# Averaged across VE x Coverage so you see one clean number per Setting x AgeCat
share_u40_o40_summary <- share_u40_o40_within_targeted_agecat %>%
  dplyr::group_by(Setting, Scenario, AgeCat_bin, u40_flag) %>%
  dplyr::summarise(share = mean(share, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = u40_flag, values_from = share,
                     values_fill = 0) %>%
  dplyr::mutate(
    Setting = factor(Setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  dplyr::arrange(Scenario, Setting)

message("\n[v3 DIAG] Share of averted cases inside each vaccinated AgeCat ",
        "(u40 bins vs o40 bins) by Setting:")
print(share_u40_o40_summary, n = Inf)

# ---- Q2: Effective chr_prop for each Setting x AgeCat, implied by bin-level
#          case weighting (this is the 'effective' number v3 uses instead of
#          v2's uniform u40 median) --------------------------------------------
chr_u40_med <- .daly_prop_vecs$meta$chr_u40
chr_o40_med <- .daly_prop_vecs$meta$chr_o40

effective_chr_prop <- share_u40_o40_summary %>%
  dplyr::mutate(
    chr_prop_v2 = chr_u40_med,
    chr_prop_v3 = u40 * chr_u40_med + o40 * chr_o40_med,
    fold_change = chr_prop_v3 / chr_prop_v2
  ) %>%
  dplyr::arrange(Scenario, Setting)

message("\n[v3 DIAG] Effective chr_prop (v2 vs v3) by Setting x AgeCat, ",
        "and v3/v2 fold change:")
print(effective_chr_prop, n = Inf)

# ---- Q3: Full age-bin share by Setting (ignores AgeCat grouping).
#          For each Setting the 20 bin shares sum to 1; the distribution across
#          bins reveals where infections actually land in that Setting.
bin_share_by_setting <- diagnostic_averted_bin %>%
  dplyr::group_by(Setting, bin, age_bin, u40_flag) %>%
  dplyr::summarise(averted = sum(averted, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(Setting) %>%                      # <-- normalise ACROSS bins
  dplyr::mutate(share = averted / sum(averted, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Setting = factor(Setting, levels = c("Low","Moderate","High"))) %>%
  dplyr::arrange(Setting, bin)

message("\n[v3 DIAG] Age-bin share of total averted cases by Setting ",
        "(all scenarios pooled):")
print(
  bin_share_by_setting %>%
    tidyr::pivot_wider(id_cols = c(bin, age_bin, u40_flag),
                       names_from = Setting, values_from = share) %>%
    dplyr::arrange(bin),
  n = Inf
)

# Handy objects kept in workspace for follow-up inspection:
#   - diagnostic_averted_bin     : per Region x VE x Cov x Scen x bin
#   - share_u40_o40_summary      : compact Setting x AgeCat x (u40/o40) shares
#   - effective_chr_prop         : implied v3 chr_prop per Setting x AgeCat
#   - bin_share_by_setting       : overall age-bin share per Setting

step_done("SECTION 09")
step_done("03_brazil_all_draws_ori_v3.R — all sections complete")

