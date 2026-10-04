### -----------------------------------------------------------------------
### Vaccination-timing (week) x coverage sweep, finite-history model,
### FULL draw-level output -- same structure as sim_results_vc_ixchiq_model /
### postsim_vc_ixchiq_model (region -> VE -> coverage -> scenario_result),
### one such pair of objects per week value, so 03_brazil_all_draws_ori_v3.R
### can consume any of them unmodified just by swapping which file it loads.
###
### Grid: week in {1,8,16,24,32,42,52} x coverage in {0.10,0.30,0.60,0.90}
### x 11 regions x 4 scenarios, n_draws configurable below (N_DRAWS).
### VE is fixed at the standard "with efficacy" branch (ve_inf = NA, i.e.
### lhs_sample$ve_ix per draw) -- VE=0 sweep is out of scope here -- but the
### output still nests under a VE tag ("VE98.9") to match the shape
### 03_brazil_all_draws_ori_v3.R expects.
###
### Uses the fixed run_simulation_scenarios_ui_ixchiq() (delay/coverage_override
### args added, and the total_coverage_frac coverage bug fixed -- see
### age_struc_fitting_region_func_updated.R). Coverage values here (30%/60%)
### don't have a pre-built LHS uncertainty column, so this sweep drives
### coverage via coverage_override (fixed fraction per cell, no per-draw
### coverage uncertainty) rather than vc_rand_draw.
###
### run_simulation_scenarios_ui_ixchiq() and postsim_all_ui() are pasted in
### verbatim below (not sourced from age_struc_fitting_region_func_updated.R)
### so each parallel worker gets a clean copy without executing that file's
### unrelated top-level code (which depends on objects -- bra_foi_states,
### regions -- this script never builds). If you change either function in
### age_struc_fitting_region_func_updated.R, re-sync the copies here.
### -----------------------------------------------------------------------

library(parallel)
library(purrr)
library(dplyr)

source("01_Script/1_1_Functions/sim_functions_final.R")  # sirv_sim_coverageSwitch_ixchiq()

## ---- Config ---------------------------------------------------------------

WEEKS      <- c(1, 8, 16, 24, 32, 42, 52)
COVERAGES  <- c(0.10, 0.30, 0.60, 0.90)
N_DRAWS    <- 300   # full pipeline uses 1000; reduced here to keep output size
                     # comparable to the existing single-coverage-point run
                     # (~1.4GB sim_results / ~2.9GB postsim) despite the much
                     # larger week x coverage grid. Bump to 1000 if size isn't
                     # a concern -- runtime scales linearly with N_DRAWS.
N_CORES    <- min(8, max(1, parallel::detectCores() - 2))
# NOTE (2026-08): capped at 8 -- detectCores()-2 (26 on this 28-core machine)
# was hitting an intermittent "invalid connection" / sendData.SOCKnode error
# from parallel::clusterExport() on this Windows PSOCK setup (a minimal
# 26-worker cluster test with a trivial object worked fine, so this is
# specific to exporting the larger objects/pasted-in functions below at high
# worker counts, not cluster creation itself). See chat record, 2026-08.

# VE mechanism for this sweep run -- set once here, everything else (nesting
# tag, output filenames) derives from it so a VE0 run can't silently
# overwrite the already-computed VE98.9 sweep files.
#   "VE98.9" -> VE_INF_FIXED = NA  (disease AND infection blocking; per-draw lhs$ve_ix)
#   "VE0"    -> VE_INF_FIXED = 0   (disease blocking only; no infection-blocking effect)
VE_TAG       <- "VE0"
VE_INF_FIXED <- if (VE_TAG == "VE0") 0 else NA
FILE_SUFFIX  <- if (VE_TAG == "VE0") "_ve0" else ""

## ---- Load shared inputs ----------------------------------------------------

load("00_Data/0_2_Processed/lhs_combined_finite.RData")   # posterior_list_finite, lhs_combined_finite, lhs_idx_list_finite, posterior_idx_list_finite
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")  # N_bahia, N_ceara, N_mg, N_pemam, N_pa, N_rg, N_pi, N_ag, N_tc, N_se, N_go
load("00_Data/0_2_Processed/observed_2022.RData")          # observed_ce, observed_bh, ...
load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData")  # hosp, fatal, nh_fatal
load("00_Data/0_2_Processed/prevacc_sets_finite.RData")    # prevacc_sets_finite (pre-vacc baseline, week/coverage independent)
load("00_Data/0_2_Processed/rho_df_finite.RData")          # rho_df (postsim_all_ui() global dependency, m3 refit)

lhs_sample_young <- readRDS("00_Data/0_2_Processed/lhs_sample_young.RDS")
lhs_old          <- readRDS("00_Data/0_2_Processed/lhs_old.RDS")
le_sample        <- readRDS("00_Data/0_2_Processed/le_sample.RDS")

age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
                mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
                mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
                mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))

age_gr_levels <- gsub("–", "-", c(
  "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
))

scenario_labels <- c("1-11y", "12-17y", "18-64y", "65y+")
target_age_list <- list(
  c(0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1)
)

region_names <- c("Bahia", "Ceará", "Minas Gerais", "Pernambuco", "Paraíba",
                  "Rio Grande do Norte", "Piauí", "Alagoas", "Tocantins",
                  "Sergipe", "Goiás")

N_by_region <- list(
  "Bahia"               = as.numeric(N_bahia[["Bahia"]]),
  "Ceará"               = as.numeric(N_ceara[["Ceará"]]),
  "Minas Gerais"        = as.numeric(N_mg[["Minas Gerais"]]),
  "Pernambuco"          = as.numeric(N_pemam[["Pernambuco"]]),
  "Paraíba"             = as.numeric(N_pa[["Paraíba"]]),
  "Rio Grande do Norte" = as.numeric(N_rg[["Rio Grande do Norte"]]),
  "Piauí"               = as.numeric(N_pi[["Piauí"]]),
  "Alagoas"             = as.numeric(N_ag[["Alagoas"]]),
  "Tocantins"           = as.numeric(N_tc[["Tocantins"]]),
  "Sergipe"             = as.numeric(N_se[["Sergipe"]]),
  "Goiás"               = as.numeric(N_go[["Goiás"]])
)

observed_by_region <- list(
  "Bahia" = observed_bh, "Ceará" = observed_ce, "Minas Gerais" = observed_mg,
  "Pernambuco" = observed_pn, "Paraíba" = observed_pa,
  "Rio Grande do Norte" = observed_rg, "Piauí" = observed_pi,
  "Alagoas" = observed_ag, "Tocantins" = observed_tc,
  "Sergipe" = observed_se, "Goiás" = observed_go
)

# posterior_idx_list_finite/lhs_idx_list_finite are built with size 1000
# (see lhs_samples_sir_finite.R) and are already a random sample/permutation,
# so truncating to the first N_DRAWS is equivalent to a random subsample.
posterior_idx_list <- lapply(posterior_idx_list_finite, function(x) x[1:N_DRAWS])
lhs_idx_list        <- lapply(lhs_idx_list_finite,        function(x) x[1:N_DRAWS])

# foi_draws_list[[region_name]] is read inside run_simulation_scenarios_ui_ixchiq()
# but the result (foi_draws) is never used afterward -- dead code carried over
# from the original function. A dummy list with the right names is enough.
foi_draws_list <- setNames(vector("list", length(region_names)), region_names)


## ---- run_simulation_scenarios_ui_ixchiq() -----------------------------
## Verbatim copy of the fixed version in age_struc_fitting_region_func_updated.R
## (L2780 as of this writing). See that file for the authoritative source.

run_simulation_scenarios_ui_ixchiq <- function(target_age_list,
                                               observed,
                                               N, bra_foi_state_summ, age_groups, region_name,
                                               hosp, fatal, nh_fatal,
                                               lhs_sample_young, lhs_old, le_sample,
                                               age_gr_levels,
                                               prevacc_ui = NULL,
                                               posterior,
                                               ve_inf = 0,
                                               total_coverage = 0,
                                               lhs_sample,
                                               finite_history = TRUE,
                                               delay = 2,
                                               coverage_override = NULL
) {
  n_scenarios <- length(target_age_list)
  scenario_result <- vector("list", n_scenarios)
  posterior_idx <- posterior_idx_list[[region_name]]
  lhs_idx       <- lhs_idx_list[[region_name]]
  n_draws <- length(posterior_idx)

  T <- nrow(observed)

  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  default_age_vector <- gsub("–", "-", default_age_vector)
  age_gr <- rep(default_age_vector, T)

  for (s in seq_along(target_age_list)) {
    target <- target_age_list[[s]]
    draw_results_raw_inf  <- vector("list", n_draws)
    draw_results_raw_symp <- vector("list", n_draws)
    draw_results_raw_alloc <- vector("list", n_draws)
    draw_results_vacc_to_S <- vector("list", n_draws)

    foi_draws <- foi_draws_list[[region_name]]
    set.seed(123)

    for (d in 1:n_draws) {

      idx_post <- posterior_idx[d]
      idx_lhs  <- lhs_idx[d]

      # NOTE (2026-08): gamma/rho/sigma re-synced from the fixed version in
      # age_struc_fitting_region_func_updated.R's run_simulation_scenarios_ui_ixchiq()
      # -- must come from `posterior` at the SAME idx_post as beta_observed/I0,
      # not from lhs_sample at idx_lhs (an unrelated, independently-permuted
      # draw). This pasted-in copy had drifted out of sync (still using
      # idx_lhs), which decoupled the transmission force (beta/I0) from the
      # recovery/reporting rates (gamma/rho/sigma) for every "post" draw --
      # confirmed via draw-level cor(baseline, post) ~ 0.76 and 3x+ outliers
      # that vanish once re-paired on idx_post. See chat record, 2026-08.
      base_beta_draw <- posterior$beta_observed[idx_post, ]
      I0_draw        <- posterior$I0[idx_post, ]
      gamma_draw     <- posterior$gamma[idx_post]
      rho_draw       <- posterior$rho[idx_post]
      sigma_draw     <- posterior$sigma[idx_post]
      FOI_rand_draw  <- lhs_sample$foi[idx_lhs]
      ve_rand_draw   <- lhs_sample$ve_ix[idx_lhs]
      wd_rand_draw   <- lhs_sample$wd[idx_lhs]

      vc_rand_draw <- case_when(
        total_coverage == 0.10 ~ lhs_sample$vc10[idx_lhs],
        total_coverage == 0.50 ~ lhs_sample$vc50[idx_lhs],
        total_coverage == 0.90 ~ lhs_sample$vc90[idx_lhs],
        TRUE ~ NA_real_
      )

      if (!is.null(coverage_override)) {
        total_coverage_frac <- coverage_override
      } else if (!is.na(vc_rand_draw)) {
        total_coverage_frac <- vc_rand_draw
      } else {
        stop(
          "total_coverage = ", total_coverage, " has no matching LHS ",
          "uncertainty column (vc10/vc50/vc90 exist for 0.10/0.50/0.90 only). ",
          "Pass coverage_override for other coverage fractions."
        )
      }

      VE_block_draw <- ve_rand_draw

      VE_inf_draw <- if (is.na(ve_inf)) {
        ve_rand_draw
      } else if (ve_inf == 0) {
        0
      } else {
        ve_inf
      }

      if (finite_history) {
        R0_draw <- 1 - exp(-FOI_rand_draw * pmin(age_groups, 2022 - 2014))
      } else {
        R0_draw <- 1 - exp(-FOI_rand_draw * age_groups)
      }

      sim_out <- sirv_sim_coverageSwitch_ixchiq(
        T = T,
        A = length(age_gr_levels),
        N = N,
        r = rep(0, length(age_gr_levels)),
        base_beta = base_beta_draw,
        I0_draw = I0_draw,
        R0 = R0_draw,
        rho = rho_draw,
        gamma = gamma_draw,
        sigma = sigma_draw,
        delay = delay,
        VE_block = VE_block_draw,
        VE_inf = VE_inf_draw,
        target_age = target,
        total_coverage = total_coverage_frac,
        weekly_delivery_speed = wd_rand_draw
      )

      draw_results_raw_inf[[d]]  <- sim_out$age_stratified_cases_raw
      draw_results_raw_symp[[d]] <- sim_out$true_symptomatic
      draw_results_raw_alloc[[d]] <- sim_out$raw_allocation_age
      draw_results_vacc_to_S[[d]] <- sim_out$vacc_delayed

      if (anyNA(sim_out$age_stratified_cases_raw))
        stop("Draw ", d, ": sim_out produced NA")
    }

    age_array_raw_inf  <- array(unlist(draw_results_raw_inf), dim = c(20, T, n_draws))
    age_array_raw_symp <- array(unlist(draw_results_raw_symp), dim = c(20, T, n_draws))
    raw_allocation_array <- array(unlist(draw_results_raw_alloc), dim = c(20, T, n_draws))
    vacc_to_S_array      <- array(unlist(draw_results_vacc_to_S), dim = c(20, T, n_draws))

    median_by_age_rawinf <- apply(age_array_raw_inf, c(1, 2), median)
    low95_by_age_rawinf  <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawinf   <- apply(age_array_raw_inf, c(1, 2), quantile, probs = 0.975)

    median_by_age_rawsymp <- apply(age_array_raw_symp, c(1, 2), median)
    low95_by_age_rawsymp  <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.025)
    hi95_by_age_rawsymp   <- apply(age_array_raw_symp, c(1, 2), quantile, probs = 0.975)

    weekly_totals_rawinf       <- apply(age_array_raw_inf, c(2, 3), sum)
    weekly_cases_median_rawinf <- apply(weekly_totals_rawinf, 1, median)
    weekly_cases_low95_rawinf  <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawinf   <- apply(weekly_totals_rawinf, 1, quantile, probs = 0.975)

    weekly_totals_rawsymp       <- apply(age_array_raw_symp, c(2, 3), sum)
    weekly_cases_median_rawsymp <- apply(weekly_totals_rawsymp, 1, median)
    weekly_cases_low95_rawsymp  <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.025)
    weekly_cases_hi95_rawsymp   <- apply(weekly_totals_rawsymp, 1, quantile, probs = 0.975)

    infection_median_df <- as.data.frame.table(median_by_age_rawinf, responseName = "infection")
    infection_low_df    <- as.data.frame.table(low95_by_age_rawinf,  responseName = "infection_lo")
    infection_hi_df     <- as.data.frame.table(hi95_by_age_rawinf,   responseName = "infection_hi")

    colnames(infection_median_df) <- c("AgeGroup", "Week", "infection")
    colnames(infection_low_df)    <- c("AgeGroup", "Week", "infection_lo")
    colnames(infection_hi_df)     <- c("AgeGroup", "Week", "infection_hi")

    infection_df <- infection_median_df %>%
      mutate(
        infection_lo = infection_low_df$infection_lo,
        infection_hi = infection_hi_df$infection_hi,
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )

    sim_df <- as.data.frame.table(median_by_age_rawsymp, responseName = "Cases")
    colnames(sim_df) <- c("AgeGroup", "Week", "Cases")
    sim_df <- sim_df %>%
      mutate(
        Scenario = s,
        AgeGroup = as.numeric(AgeGroup),
        Week = as.numeric(Week)
      )
    sim_df <- sim_df %>% mutate(region_full = region_name)
    sim_df <- left_join(sim_df, infection_df, by = c("AgeGroup", "Week"))

    low95_df <- as.data.frame.table(low95_by_age_rawsymp, responseName = "post_low95")
    hi95_df  <- as.data.frame.table(hi95_by_age_rawsymp, responseName = "post_hi95")
    sim_df <- sim_df %>%
      mutate(
        post_low95 = low95_df$post_low95,
        post_hi95  = hi95_df$post_hi95
      )

    weekly_df <- data.frame(
      Week = 1:T,
      weekly_median = weekly_cases_median_rawsymp,
      weekly_low95 = weekly_cases_low95_rawsymp,
      weekly_hi95 = weekly_cases_hi95_rawsymp
    )

    sim_df <- sim_df %>% mutate(age_gr = rep(default_age_vector, T))
    sim_df$age_gr <- factor(sim_df$age_gr, levels = age_gr_levels)

    sim_df <- sim_df %>%
      mutate(
        hosp_rate = rep(hosp, T),
        hospitalised = Cases * hosp_rate,
        hospitalised_lo = post_low95 * hosp_rate,
        hospitalised_hi = post_hi95 * hosp_rate,
        non_hospitalised = Cases - hospitalised,
        non_hospitalised_lo = post_low95 - hospitalised_lo,
        non_hospitalised_hi = post_hi95 - hospitalised_hi,
        fatality = rep(fatal, T),
        nh_fatality = rep(nh_fatal, T),
        fatal = (hospitalised * fatality + non_hospitalised * nh_fatality),
        fatal_lo = (hospitalised_lo * fatality + non_hospitalised_lo * nh_fatality),
        fatal_hi = (hospitalised_hi * fatality + non_hospitalised_hi * nh_fatality)
      ) %>%
      arrange(Week, AgeGroup) %>%
      group_by(Week, AgeGroup) %>%
      mutate(
        cum_fatal = cumsum(fatal),
        cum_hosp = cumsum(hospitalised)
      ) %>%
      ungroup()

    sim_df <- sim_df %>%
      mutate(
        age_numeric = case_when(
          age_gr == "<1"    ~ 0,
          age_gr == "1-4"   ~ 2.5,
          age_gr == "5-9"   ~ 7,
          age_gr == "10-11" ~ 10.5,
          age_gr == "12-17" ~ 14.5,
          age_gr == "18-19" ~ 18.5,
          age_gr == "20-24" ~ 22,
          age_gr == "25-29" ~ 27,
          age_gr == "30-34" ~ 32,
          age_gr == "35-39" ~ 37,
          age_gr == "40-44" ~ 42,
          age_gr == "45-49" ~ 47,
          age_gr == "50-54" ~ 52,
          age_gr == "55-59" ~ 57,
          age_gr == "60-64" ~ 62,
          age_gr == "65-69" ~ 67,
          age_gr == "70-74" ~ 72,
          age_gr == "75-79" ~ 77,
          age_gr == "80-84" ~ 82,
          age_gr == "85+"   ~ 87.5,
          TRUE              ~ NA_real_
        )
      )

    sim_df <- sim_df %>%
      mutate(
        dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5),
        dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
        dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
        dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5),
        dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
        dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
        dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5)
      ) %>%
      mutate(
        subac_prop = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$subac, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$subac, 0.5)
        ),
        chr_6m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr6m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr6m, 0.5)
        ),
        chr_12m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr12m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr12m, 0.5)
        ),
        chr_30m = case_when(
          age_numeric < 40 ~ quantile(lhs_sample_young$chr30m, 0.5),
          age_numeric >= 40 ~ quantile(lhs_old$chr30m, 0.5)
        ),
        chr_prop = chr_6m + chr_12m + chr_30m
      ) %>%
      mutate(
        yld_acute = (hospitalised * dw_hosp * dur_acute) +
          (non_hospitalised * dw_nonhosp * dur_acute),
        yld_acute_lo = (hospitalised_lo * dw_hosp * dur_acute) +
          ((post_low95 - hospitalised_lo) * dw_nonhosp * dur_acute),
        yld_acute_hi = (hospitalised_hi * dw_hosp * dur_acute) +
          ((post_hi95 - hospitalised_hi) * dw_nonhosp * dur_acute),

        yld_subacute = (hospitalised * subac_prop * dw_subacute * dur_subacute) +
          (non_hospitalised * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_lo = (hospitalised_lo * subac_prop * dw_subacute * dur_subacute) +
          ((post_low95 - hospitalised_lo) * chr_prop * dw_subacute * dur_subacute),
        yld_subacute_hi = (hospitalised_hi * subac_prop * dw_subacute * dur_subacute) +
          ((post_hi95 - hospitalised_hi) * chr_prop * dw_subacute * dur_subacute),

        yld_chronic = (hospitalised * chr_prop * dw_chronic * dur_chronic) +
          (non_hospitalised * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_lo = (hospitalised_lo * chr_prop * dw_chronic * dur_chronic) +
          ((post_low95 - hospitalised_lo) * chr_prop * dw_chronic * dur_chronic),
        yld_chronic_hi = (hospitalised_hi * chr_prop * dw_chronic * dur_chronic) +
          ((post_hi95 - hospitalised_hi) * chr_prop * dw_chronic * dur_chronic),

        yld_total = yld_acute + yld_subacute + yld_chronic,
        yld_total_lo = yld_acute_lo + yld_subacute_lo + yld_chronic_lo,
        yld_total_hi = yld_acute_hi + yld_subacute_hi + yld_chronic_hi
      ) %>%
      mutate(
        le_left = case_when(
          age_numeric <= 1 ~ quantile(le_sample$le_1, 0.5),
          age_numeric > 1 & age_numeric < 10 ~ quantile(le_sample$le_1, 0.5),
          age_numeric >= 10 & age_numeric < 20 ~ quantile(le_sample$le_2, 0.5),
          age_numeric >= 20 & age_numeric < 30 ~ quantile(le_sample$le_2, 0.5),
          age_numeric >= 30 & age_numeric < 40 ~ quantile(le_sample$le_3, 0.5),
          age_numeric >= 40 & age_numeric < 50 ~ quantile(le_sample$le_4, 0.5),
          age_numeric >= 50 & age_numeric < 60 ~ quantile(le_sample$le_5, 0.5),
          age_numeric >= 60 & age_numeric < 70 ~ quantile(le_sample$le_6, 0.5),
          age_numeric >= 70 & age_numeric < 80 ~ quantile(le_sample$le_7, 0.5),
          age_numeric >= 80 ~ quantile(le_sample$le_8, 0.5),
          TRUE ~ NA_real_
        )
      ) %>%
      mutate(
        yll = fatal * le_left,
        yll_lo = fatal_lo * le_left,
        yll_hi = fatal_hi * le_left,
        daly_tot = yld_total + yll,
        daly_tot_lo = yld_total_lo + yll_lo,
        daly_tot_hi = yld_total_hi + yll_hi,
        cum_daly = cumsum(daly_tot)
      )

    if (!is.null(prevacc_ui)) {
      sim_df <- left_join(sim_df,
                          dplyr::select(prevacc_ui, AgeGroup, Week,
                                        pre_Median = Median,
                                        pre_low95 = low95,
                                        pre_hi95 = hi95),
                          by = c("AgeGroup", "Week"))
    }

    scenario_result[[s]] <- list(sim_result = list(age_array_raw_symp           = age_array_raw_symp,
                                                   age_array_raw_inf            = age_array_raw_inf,
                                                   raw_allocation_array         = raw_allocation_array,
                                                   vacc_to_S_array              = vacc_to_S_array,
                                                   weekly_cases_median_rawsymp  = weekly_cases_median_rawsymp,
                                                   weekly_cases_median_rawinf   = weekly_cases_median_rawinf,
                                                   weekly_cases_low95_rawsymp   = weekly_cases_low95_rawsymp,
                                                   weekly_cases_low95_rawinf    = weekly_cases_low95_rawinf,
                                                   weekly_cases_hi95_rawsymp    = weekly_cases_hi95_rawsymp),
                                 sim_out                      = sim_out,
                                 sim_df                       = sim_df,
                                 weekly_df                    = weekly_df
    )
  }

  return(scenario_result)
}


## ---- postsim_all_ui() ------------------------------------------------
## Verbatim copy of age_struc_fitting_region_func_updated.R's active
## definition (L4419 as of this writing) -- unmodified.

postsim_all_ui <- function(scenario_result,
                           observed,
                           age_gr_levels,
                           pre_summary_cases_age,
                           pre_summary_cases,
                           pre_summary_cases_all,
                           region) {
  T <- nrow(observed)
  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )
  default_age_vector <- gsub("–", "-", default_age_vector)
  age_gr <- rep(default_age_vector, T)
  n_scenarios <- length(scenario_result)
  scenario_df  <- lapply(scenario_result, function(x) x$sim_df)

  summary_list <- lapply(seq_along(scenario_df), function(idx) {
    df <- scenario_df[[idx]]
    df$pre_infection      <- pre_summary_cases$Infections
    df$pre_infection_lo   <- pre_summary_cases$Infections_lo
    df$pre_infection_hi   <- pre_summary_cases$Infections_hi
    df$pre_vacc  <- pre_summary_cases$Median
    df$pre_fatal <- pre_summary_cases$fatal
    df$pre_daly  <- pre_summary_cases$daly_tot
    df <- df %>%
      mutate(
        diff   = pre_vacc - Cases,
        impact = diff / pre_vacc * 100
      )
    return(df)
  })

  summary_list_df <- do.call(rbind, summary_list)
  summary_list_df$age_gr <- rep(age_gr, length.out = nrow(summary_list_df))
  summary_list_df$age_gr <- factor(summary_list_df$age_gr, levels = age_gr_levels)

  if (all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      distinct(Week, .keep_all = TRUE)
    summary_list_df <- left_join(summary_list_df, pre_weekly_df, by = "Week")
  }

  final_summ <- lapply(summary_list, function(df) {
    df %>%
      group_by(AgeGroup) %>%
      summarise(
        infection        = sum(infection, na.rm = TRUE),
        infection_lo     = sum(infection_lo, na.rm = TRUE),
        infection_hi     = sum(infection_hi, na.rm = TRUE),
        Median           = sum(Cases, na.rm = TRUE),
        low95            = sum(post_low95, na.rm = TRUE),
        hi95             = sum(post_hi95, na.rm = TRUE),
        hospitalised     = sum(hospitalised, na.rm = TRUE),
        hospitalised_lo  = sum(hospitalised_lo, na.rm = TRUE),
        hospitalised_hi  = sum(hospitalised_hi, na.rm = TRUE),
        yld_acute        = sum(yld_acute, na.rm = TRUE),
        yld_acute_lo     = sum(yld_acute_lo, na.rm = TRUE),
        yld_acute_hi     = sum(yld_acute_hi, na.rm = TRUE),
        yld_subacute     = sum(yld_subacute, na.rm = TRUE),
        yld_subacute_lo  = sum(yld_subacute_lo, na.rm = TRUE),
        yld_subacute_hi  = sum(yld_subacute_hi, na.rm = TRUE),
        yld_chronic      = sum(yld_chronic, na.rm = TRUE),
        yld_chronic_lo   = sum(yld_chronic_lo, na.rm = TRUE),
        yld_chronic_hi   = sum(yld_chronic_hi, na.rm = TRUE),
        yld_total        = sum(yld_total, na.rm = TRUE),
        yld_total_lo     = sum(yld_total_lo, na.rm = TRUE),
        yld_total_hi     = sum(yld_total_hi, na.rm = TRUE),
        yll              = sum(yll, na.rm = TRUE),
        yll_lo           = sum(yll_lo, na.rm = TRUE),
        yll_hi           = sum(yll_hi, na.rm = TRUE),
        daly_tot         = sum(daly_tot, na.rm = TRUE),
        daly_tot_lo      = sum(daly_tot_lo, na.rm = TRUE),
        daly_tot_hi      = sum(daly_tot_hi, na.rm = TRUE),
        fatal            = sum(fatal, na.rm = TRUE),
        fatal_lo         = sum(fatal_lo, na.rm = TRUE),
        fatal_hi         = sum(fatal_hi, na.rm = TRUE)
      ) %>%
      mutate(
        pre_infection         = pre_summary_cases_age$Infections,
        pre_infection_lo      = pre_summary_cases_age$Infections_lo,
        pre_infection_hi      = pre_summary_cases_age$Infections_hi,

        pre_vacc              = pre_summary_cases_age$Median,
        pre_vacc_low95        = pre_summary_cases_age$low95,
        pre_vacc_hi95         = pre_summary_cases_age$hi95,

        pre_yld_acute         = pre_summary_cases_age$yld_acute,
        pre_yld_acute_low95   = pre_summary_cases_age$yld_acute_lo,
        pre_yld_acute_hi      = pre_summary_cases_age$yld_acute_hi,

        pre_yld_subacute      = pre_summary_cases_age$yld_subacute,
        pre_yld_subacute_low95= pre_summary_cases_age$yld_subacute_lo,
        pre_yld_subacute_hi   = pre_summary_cases_age$yld_subacute_hi,

        pre_yld_chronic       = pre_summary_cases_age$yld_chronic,
        pre_yld_chronic_low95 = pre_summary_cases_age$yld_chronic_lo,
        pre_yld_chronic_hi    = pre_summary_cases_age$yld_chronic_hi,

        pre_yld_total         = pre_summary_cases_age$yld_total,
        pre_yld_total_low95   = pre_summary_cases_age$yld_total_lo,
        pre_yld_total_hi      = pre_summary_cases_age$yld_total_hi,

        pre_yll               = pre_summary_cases_age$yll,
        pre_yll_low95         = pre_summary_cases_age$yll_lo,
        pre_yll_hi            = pre_summary_cases_age$yll_hi,

        pre_daly              = pre_summary_cases_age$daly_tot,
        pre_daly_low95        = pre_summary_cases_age$daly_tot_lo,
        pre_daly_hi           = pre_summary_cases_age$daly_tot_hi,

        pre_fatal             = pre_summary_cases_age$fatal,
        pre_fatal_low95       = pre_summary_cases_age$fatal_lo,
        pre_fatal_hi          = pre_summary_cases_age$fatal_hi,

        pre_hospitalised       = pre_summary_cases_age$hospitalised,
        pre_hospitalised_low95 = pre_summary_cases_age$hospitalised_lo,
        pre_hospitalised_hi    = pre_summary_cases_age$hospitalised_hi,

        diff_inf    = pre_infection - infection,
        diff_inf_lo = pre_infection_lo - infection_lo,
        diff_inf_hi = pre_infection_hi - infection_hi,
        diff       = pre_vacc - Median,
        diff_low   = pre_vacc_low95 - low95,
        diff_hi    = pre_vacc_hi95 - hi95,
        impact     = diff / pre_vacc * 100,
        impact_low = diff_low / pre_vacc_low95 * 100,
        impact_hi  = diff_hi / pre_vacc_hi95 * 100,

        diff_fatal       = pre_fatal - fatal,
        diff_fatal_low   = pre_fatal_low95 - fatal_lo,
        diff_fatal_hi    = pre_fatal_hi - fatal_hi,
        impact_fatal     = diff_fatal / pre_fatal * 100,
        impact_fatal_low = diff_fatal_low / pre_fatal_low95 * 100,
        impact_fatal_hi  = diff_fatal_hi / pre_fatal_hi * 100,

        diff_hosp       = pre_hospitalised - hospitalised,
        diff_hosp_low   = pre_hospitalised_low95 - hospitalised_lo,
        diff_hosp_hi    = pre_hospitalised_hi - hospitalised_hi,
        impact_hosp     = diff_hosp / pre_hospitalised * 100,
        impact_hosp_low = diff_hosp_low / pre_hospitalised_low95 * 100,
        impact_hosp_hi  = diff_hosp_hi / pre_hospitalised_hi * 100,

        diff_yld_acute       = pre_yld_acute - yld_acute,
        diff_yld_acute_low   = pre_yld_acute_low95 - yld_acute_lo,
        diff_yld_acute_hi    = pre_yld_acute_hi - yld_acute_hi,
        impact_yld_acute     = diff_yld_acute / pre_yld_acute * 100,
        impact_yld_acute_low = diff_yld_acute_low / pre_yld_acute_low95 * 100,
        impact_yld_acute_hi  = diff_yld_acute_hi / pre_yld_acute_hi * 100,

        diff_yld_subacute       = pre_yld_subacute - yld_subacute,
        diff_yld_subacute_low   = pre_yld_subacute_low95 - yld_subacute_lo,
        diff_yld_subacute_hi    = pre_yld_subacute_hi - yld_subacute_hi,
        impact_yld_subacute     = diff_yld_subacute / pre_yld_subacute * 100,
        impact_yld_subacute_low = diff_yld_subacute_low / pre_yld_subacute_low95 * 100,
        impact_yld_subacute_hi  = diff_yld_subacute_hi / pre_yld_subacute_hi * 100,

        diff_yld_chronic       = pre_yld_chronic - yld_chronic,
        diff_yld_chronic_low   = pre_yld_chronic_low95 - yld_chronic_lo,
        diff_yld_chronic_hi    = pre_yld_chronic_hi - yld_chronic_hi,
        impact_yld_chronic     = diff_yld_chronic / pre_yld_chronic * 100,
        impact_yld_chronic_low = diff_yld_chronic_low / pre_yld_chronic_low95 * 100,
        impact_yld_chronic_hi  = diff_yld_chronic_hi / pre_yld_chronic_hi * 100,

        diff_yld_total       = pre_yld_total - yld_total,
        diff_yld_total_low   = pre_yld_total_low95 - yld_total_lo,
        diff_yld_total_hi    = pre_yld_total_hi - yld_total_hi,
        impact_yld_total     = diff_yld_total / pre_yld_total * 100,
        impact_yld_total_low = diff_yld_total_low / pre_yld_total_low95 * 100,
        impact_yld_total_hi  = diff_yld_total_hi / pre_yld_total_hi * 100,

        diff_yll       = pre_yll - yll,
        diff_yll_low   = pre_yll_low95 - yll_lo,
        diff_yll_hi    = pre_yll_hi - yll_hi,
        impact_yll     = diff_yll / pre_yll * 100,
        impact_yll_low = diff_yll_low / pre_yll_low95 * 100,
        impact_yll_hi  = diff_yll_hi / pre_yll_hi * 100,

        diff_daly       = pre_daly - daly_tot,
        diff_daly_low   = pre_daly_low95 - daly_tot_lo,
        diff_daly_hi    = pre_daly_hi - daly_tot_hi,
        impact_daly     = diff_daly / pre_daly * 100,
        impact_daly_low = diff_daly_low / pre_daly_low95 * 100,
        impact_daly_hi  = diff_daly_hi / pre_daly_hi * 100
      )
  })

  summary_week <- lapply(seq_along(summary_list), function(i) {
    summary_list[[i]] %>%
      mutate(Scenario = paste0("Scenario_", i),
             pre_vacc = pre_summary_cases$Median,
             pre_vacc_lo = pre_summary_cases$low95,
             pre_vacc_hi = pre_summary_cases$hi95,
             pre_fatal   = pre_summary_cases$fatal,
             pre_fatal_lo   = pre_summary_cases$fatal_lo,
             pre_fatal_hi   = pre_summary_cases$fatal_hi,
             pre_daly       = pre_summary_cases$daly_tot,
             pre_daly_lo    = pre_summary_cases$daly_tot_lo,
             pre_daly_hi    = pre_summary_cases$daly_tot_hi,
             ) %>%
      group_by(Week, Scenario) %>%
      summarise(
        post_inf     = sum(infection),
        post_inf_lo  = sum(infection_lo),
        post_inf_hi  = sum(infection_hi),
        pre_inf      = sum(pre_infection),
        pre_inf_lo   = sum(pre_infection_lo),
        pre_inf_hi   = sum(pre_infection_hi),
        post_cases = sum(Cases, na.rm = TRUE),
        post_cases_lo = sum(post_low95, na.rm = TRUE),
        post_cases_hi = sum(post_hi95, na.rm = TRUE),
        pre_cases  = sum(pre_vacc, na.rm = TRUE),
        pre_cases_lo = sum(pre_vacc_lo, na.rm = TRUE),
        pre_cases_hi  = sum(pre_vacc_hi, na.rm = TRUE),
        diff       = pre_cases - post_cases,
        post_fatal = sum(fatal, na.rm = TRUE),
        post_fatal_lo = sum(fatal_lo, na.rm = TRUE),
        post_fatal_hi = sum(fatal_hi, na.rm = TRUE),

        pre_fatal  = sum(pre_fatal, na.rm = TRUE),
        pre_fatal_lo = sum(pre_fatal_lo, na.rm = TRUE),
        pre_fatal_hi = sum(pre_fatal_hi, na.rm = TRUE),
        diff_fatal = pre_fatal - post_fatal,
        post_daly  = sum(daly_tot, na.rm = TRUE),
        post_daly_lo  = sum(daly_tot_lo, na.rm = TRUE),
        post_daly_hi  = sum(daly_tot_hi, na.rm = TRUE),

        pre_daly   = sum(pre_daly, na.rm = TRUE),
        pre_daly_lo = sum(pre_daly_lo, na.rm = TRUE),
        pre_daly_hi = sum(pre_daly_hi, na.rm = TRUE),
        diff_daly  = pre_daly - post_daly,
        Scenario   = first(Scenario),
        .groups = "drop"
      )
  })

  summary_week_df <- do.call(rbind, summary_week)

  for (i in seq_len(n_scenarios)) {
    scenario_ui <- scenario_result[[i]]$weekly_df %>%
      dplyr::rename(
        post_weekly_median = weekly_median,
        post_weekly_low95  = weekly_low95,
        post_weekly_hi95   = weekly_hi95
      ) %>%
      distinct(Week, .keep_all = TRUE)
    summary_week[[i]] <- left_join(summary_week[[i]], scenario_ui, by = "Week", relationship = "many-to-many")
  }

  summary_week_df <- do.call(rbind, summary_week)

  if (all(c("lo95", "hi95") %in% colnames(pre_summary_cases_all))) {
    pre_weekly_df <- pre_summary_cases_all %>%
      dplyr::select(Week, lo95, hi95) %>%
      distinct(Week, .keep_all = TRUE)
    summary_week_df <- left_join(summary_week_df, pre_weekly_df, by = "Week", relationship = "many-to-many")
  }

  summary_week_df$region <- region

  # NOTE (2026-08): removed /rho_p50 post-scaling -- post_cases/pre_cases/etc.
  # are already TRUE burden (p_sym-only, no rho) out of the m3-refit
  # simulator, so dividing by rho_p50 here would re-inflate by ~1/rho
  # (~7-10x). Same fix as applied to postsim_all_ui() in
  # age_struc_fitting_region_func_updated.R -- see chat record, 2026-08.

  global_impact_week <- summary_week_df %>%
    dplyr::group_by(Scenario) %>%
    dplyr::summarise(
      total_pre  = sum(pre_cases, na.rm=TRUE),
      total_post = sum(post_cases, na.rm=TRUE),
      total_pre_lo = sum(pre_cases_lo, na.rm = TRUE),
      total_pre_hi = sum(pre_cases_hi, na.rm = TRUE),
      total_post_lo = sum(post_cases_lo, na.rm = TRUE),
      total_post_hi = sum(post_cases_hi, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(impact = (total_pre - total_post)/total_pre *100,
                  impact_lo = (total_pre_lo - total_post_lo) / total_pre_lo * 100,
                  impact_hi = (total_pre_hi - total_post_hi) / total_pre_hi * 100
                  )
  annotation_text <- paste0(global_impact_week$Scenario, ": ", round(global_impact_week$impact,1), "%", collapse="\n")
  global_impact_week$region <- region

  vacc_start_week_s1 <- scenario_result[[1]]$sim_out$vacc_start_week[3]
  vacc_end_week_s1   <- scenario_result[[1]]$sim_out$vacc_end_week[3]
  vacc_start_week_s2 <- scenario_result[[2]]$sim_out$vacc_start_week[5]
  vacc_end_week_s2   <- scenario_result[[2]]$sim_out$vacc_end_week[5]
  vacc_start_week_s3 <- scenario_result[[3]]$sim_out$vacc_start_week[13]
  vacc_end_week_s3   <- scenario_result[[3]]$sim_out$vacc_end_week[13]

  return(list(
    scenario_result    = scenario_result,
    scenario_data      = lapply(scenario_result, function(x) x$sim_result),
    scenario_df        = scenario_df,
    summary_list       = summary_list,
    summary_list_df    = summary_list_df,
    final_summ         = final_summ,
    summary_week       = summary_week,
    summary_week_df    = summary_week_df,
    global_impact_week = global_impact_week,
    annotation_text    = annotation_text,
    vacc_weeks = list(
      scenario1 = list(start = vacc_start_week_s1, end = vacc_end_week_s1),
      scenario2 = list(start = vacc_start_week_s2, end = vacc_end_week_s2),
      scenario3 = list(start = vacc_start_week_s3, end = vacc_end_week_s3)
    )
  ))
}


## ---- Build task list (one per region x coverage, for a given week) --------

build_tasks <- function(week_val) {
  purrr::map(region_names, function(region_name) {
    purrr::map(COVERAGES, function(cov_val) {
      list(
        region_name = region_name,
        week        = week_val,
        coverage    = cov_val,
        N           = N_by_region[[region_name]],
        observed    = observed_by_region[[region_name]],
        posterior   = posterior_list_finite[[region_name]],
        lhs_sample  = lhs_combined_finite[[region_name]],
        prevacc     = prevacc_sets_finite[[region_name]]
      )
    })
  }) %>% purrr::flatten()
}

run_one_task <- function(task) {
  scenario_result <- run_simulation_scenarios_ui_ixchiq(
    target_age_list    = target_age_list,
    observed           = task$observed,
    N                  = task$N,
    bra_foi_state_summ = NULL,   # unused when finite_history = TRUE
    age_groups         = age_groups,
    region_name        = task$region_name,
    hosp               = hosp,
    fatal              = fatal,
    nh_fatal           = nh_fatal,
    lhs_sample_young   = lhs_sample_young,
    lhs_old            = lhs_old,
    le_sample          = le_sample,
    age_gr_levels      = age_gr_levels,
    prevacc_ui         = task$prevacc$prevacc_ui,
    posterior          = task$posterior,
    ve_inf             = VE_INF_FIXED,  # VE_TAG / VE_INF_FIXED set at top of script -- 0 = disease-blocking only (VE0), NA = with efficacy (VE98.9)
    total_coverage     = task$coverage,
    lhs_sample         = task$lhs_sample,
    finite_history     = TRUE,
    delay              = task$week,
    coverage_override  = task$coverage
  )

  postsim <- postsim_all_ui(
    scenario_result       = scenario_result,
    observed              = task$observed,
    age_gr_levels         = age_gr_levels,
    pre_summary_cases_age = task$prevacc$pre_summary_age,
    pre_summary_cases     = task$prevacc$prevacc_ui,
    pre_summary_cases_all = task$prevacc$prevacc_ui_all,
    region                = task$region_name
  )

  list(region = task$region_name, coverage = task$coverage,
       sim_results = scenario_result, postsim = postsim)
}


## ---- Run: one week at a time, parallelized across (region x coverage) -----

dir.create("00_Data/0_2_Processed", showWarnings = FALSE, recursive = TRUE)

cl <- parallel::makeCluster(N_CORES)
on.exit(parallel::stopCluster(cl), add = TRUE)

parallel::clusterEvalQ(cl, { library(dplyr); library(purrr) })
parallel::clusterExport(cl, c(
  "run_simulation_scenarios_ui_ixchiq", "postsim_all_ui",
  "sirv_sim_coverageSwitch_ixchiq",
  "posterior_idx_list", "lhs_idx_list", "foi_draws_list",
  "hosp", "fatal", "nh_fatal", "lhs_sample_young", "lhs_old", "le_sample",
  "age_gr_levels", "age_groups", "rho_df", "target_age_list",
  "VE_INF_FIXED",
  "run_one_task"
))

for (week_val in WEEKS) {
  message(sprintf("[finite_weeksweep_fulldraws] week = %d -- starting (%d cores)", week_val, N_CORES))
  t0 <- Sys.time()

  tasks <- build_tasks(week_val)
  results <- parallel::parLapply(cl, tasks, run_one_task)

  cov_tags <- paste0("cov", COVERAGES * 100)

  # region -> VE -> coverage -> scenario_result, matching production nesting.
  sim_results_vc_ixchiq_model <- setNames(
    lapply(region_names, function(region_name) {
      region_results <- Filter(function(r) r$region == region_name, results)
      setNames(list(setNames(lapply(region_results, `[[`, "sim_results"), cov_tags)), VE_TAG)
    }),
    region_names
  )

  postsim_vc_ixchiq_model <- setNames(
    lapply(region_names, function(region_name) {
      region_results <- Filter(function(r) r$region == region_name, results)
      setNames(list(setNames(lapply(region_results, `[[`, "postsim"), cov_tags)), VE_TAG)
    }),
    region_names
  )

  save(sim_results_vc_ixchiq_model,
       file = sprintf("00_Data/0_2_Processed/sim_results_vc_ixchiq_model_finite_week%d%s.RData", week_val, FILE_SUFFIX))
  save(postsim_vc_ixchiq_model,
       file = sprintf("00_Data/0_2_Processed/postsim_vc_ixchiq_model_finite_week%d%s.RData", week_val, FILE_SUFFIX))

  t1 <- Sys.time()
  message(sprintf("[finite_weeksweep_fulldraws] week = %d -- done (%.1f min)",
                  week_val, as.numeric(t1 - t0, units = "mins")))
}

message("[finite_weeksweep_fulldraws] all weeks complete.")
