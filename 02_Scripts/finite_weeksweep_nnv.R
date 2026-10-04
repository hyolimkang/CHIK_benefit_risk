### -----------------------------------------------------------------------
### NNV (number needed to vaccinate) step for the week x coverage sweep --
### the layer that was missing between postsim_vc_ixchiq_model_finite_week*.RData
### (produced by finite_weeksweep_fulldraws.R) and Section 08's BRR calc in
### 03_brazil_all_draws_ori_v3.R (which consumes combined_nnv_df_region_coverage_model,
### not postsim directly).
###
### Mirrors model_0803_finite_history.R's NNV section (vacc_allocation() +
### nnv_list(), then the region/VE/coverage stack into
### combined_nnv_df_region_coverage_model) exactly, just looped over the
### week-specific postsim files instead of a single one.
###
### vacc_allocation() and nnv_list() are pasted in verbatim from
### age_struc_fitting_region_func_updated.R (L5200 and L5685 as of this
### writing) for the same reason finite_weeksweep_fulldraws.R pastes in
### run_simulation_scenarios_ui_ixchiq()/postsim_all_ui() -- avoids sourcing
### that file's unrelated top-level code.
###
### Run this AFTER finite_weeksweep_fulldraws.R has produced
### postsim_vc_ixchiq_model_finite_week*.RData for the weeks you want.
### -----------------------------------------------------------------------

library(purrr)
library(dplyr)
library(tidyr)          # pivot_longer() in vacc_allocation()
library(RColorBrewer)  # vacc_allocation() calls brewer.pal() (unused downstream, but must not error)

WEEKS <- c(1, 8, 16, 24, 32, 42, 52)  # must match finite_weeksweep_fulldraws.R's WEEKS
FILE_SUFFIX <- "_ve0"  # must match the VE_TAG run of finite_weeksweep_fulldraws.R being consumed ("" for VE98.9, "_ve0" for VE0)

load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")  # N_*
load("00_Data/0_2_Processed/observed_2022.RData")          # observed_*
load("00_Data/0_2_Processed/rho_df_finite.RData")          # rho_df (nnv_list() global dependency, m3 refit)

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

age_gr_levels <- gsub("–", "-", c(
  "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
))

# setting_key: verbatim from model_0803_finite_history.R
setting_key <- c(
  "Ceará"               = "High",
  "Bahia"               = "Low",
  "Paraíba"             = "High",
  "Pernambuco"          = "Moderate",
  "Rio Grande do Norte" = "Low",
  "Piauí"               = "High",
  "Tocantins"           = "Moderate",
  "Alagoas"             = "High",
  "Minas Gerais"        = "Low",
  "Sergipe"             = "Low",
  "Goiás"               = "Low"
)

T_weeks <- 52
n_scenarios <- 4
age_gr <- rep(age_gr_levels, T_weeks)  # global read by vacc_allocation()/nnv_list()

target_age_list <- list(  # global read by nnv_list()
  c(0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1)
)


## ---- vacc_allocation() -------------------------------------------------
## Verbatim copy of age_struc_fitting_region_func_updated.R L5200.

vacc_allocation <- function(postsim_all_ui, observed, region) {

  T <- nrow(observed)

  scenario_data <- lapply(seq_along(postsim_all_ui$scenario_result), function(idx){
    list <- postsim_all_ui$scenario_result[[idx]]
    data <- list$sim_out
    return(data)
  })

  raw_allocation_week <- lapply(scenario_data, function(scenario){
    scenario$raw_allocation_age
  })
  scenario_names <- c("Scenario1", "Scenario2", "Scenario3", "Scenario4")

  weekly_allocation_df <- do.call(rbind, lapply(seq_along(raw_allocation_week), function(idx) {
    df <- as.data.frame(raw_allocation_week[[idx]])
    df$Scenario <- scenario_names[idx]
    df$AgeGroup <- 1:nrow(df)
    df
  }))

  colnames(weekly_allocation_df)[1:(ncol(weekly_allocation_df)-2)] <- as.character(1:(ncol(weekly_allocation_df)-2))

  weekly_allocation_df$tot_sum <- rowSums(weekly_allocation_df[,1:T])

  weekly_allocation_long <- weekly_allocation_df %>%
    pivot_longer(
      cols = all_of(as.character(1:T)),
      names_to = "Week",
      values_to = "Vaccinated"
    ) %>%
    mutate(
      Week = as.numeric(Week)
    ) %>%
    mutate(age_gr = rep(age_gr[1:length(age_gr_levels)], each = T, times = n_scenarios)) %>%
    mutate(age_gr = factor(age_gr, levels = unique(age_gr)))

  weekly_allocation_long$region <- region

  total_allocations <- weekly_allocation_long %>% group_by(Scenario) %>% summarise(
    tot_sum = sum(tot_sum)
  )

  return(list(
    scenario_data      = scenario_data,
    raw_allocation_week = raw_allocation_week,
    weekly_allocation_df = weekly_allocation_df,
    weekly_allocation_long = weekly_allocation_long,
    total_allocations = total_allocations
  ))
}


## ---- nnv_list() ----------------------------------------------------
## Verbatim copy of age_struc_fitting_region_func_updated.R L5685 (its
## active/last-defined version).

nnv_list <- function(vacc_allocation,
                     postsim_all_ui,
                     N,
                     region,
                     observed) {

  T <- nrow(observed)
  n_scenarios <- length(target_age_list)

  default_age_vector <- c(
    "<1", "1-4", "5–9", "10-11", "12-17", "18–19", "20–24", "25–29",
    "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
    "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
  )

  default_age_vector <- gsub("[–—]", "-", default_age_vector)
  age_gr_levels <- gsub("[–—]", "-", age_gr_levels)

  age_gr <- rep(default_age_vector, each = T)

  raw_allocation <- lapply(vacc_allocation$scenario_data, function(list){
    raw_allocation <- as.data.frame(list$raw_allocation_age)
    raw_allocation$tot_vacc <- rowSums(raw_allocation)
    list$raw_allocation <- raw_allocation
    tot_df <- setNames(as.data.frame(list$raw_allocation$tot_vacc), "tot_vacc")
  })

  raw_allocation_age <- lapply(seq_along(raw_allocation), function(id){
    df <- raw_allocation[[id]]
    tot_vacc <- sum(df$tot_vacc)
    return(tot_vacc)
  })

  final_summ_df <- do.call(rbind, lapply(seq_along(postsim_all_ui$final_summ), function(idx){

    target <- c("<1 (novacc)", "1-11 years", "1-11 years", "1-11 years", "12-17 years",
                "18-59 years", "18-59 years", "18-59 years", "18-59 years",
                "18-59 years", "18-59 years", "18-59 years", "18-59 years", "18-59 years", "18-59 years",
                "60+ years", "60+ years", "60+ years", "60+ years", "60+ years")

    df <- postsim_all_ui$final_summ[[idx]]
    df <- df %>% mutate(
      tot_vacc = as.numeric(unlist(raw_allocation[[idx]])),
      scenario = paste0("Scenario_", idx),
      target   = target
    )

    summary_by_age <- df %>%
      group_by(scenario, AgeGroup) %>%
      summarise(
        tot_vacc      = sum(tot_vacc,      na.rm = TRUE),
        pre_inf       = sum(pre_infection),
        pre_inf_lo    = sum(pre_infection_lo),
        pre_inf_hi    = sum(pre_infection_hi),
        post_inf      = sum(infection),
        post_inf_lo   = sum(infection_lo),
        post_inf_hi   = sum(infection_hi),
        pre_vacc      = sum(pre_vacc),
        pre_vacc_lo   = sum(pre_vacc_low95),
        pre_vacc_hi   = sum(pre_vacc_hi95),
        post_vacc     = sum(Median),
        post_vacc_lo  = sum(low95),
        post_vacc_hi  = sum(hi95),
        pre_fatal     = sum(pre_fatal),
        pre_fatal_lo  = sum(pre_fatal_low95),
        pre_fatal_hi  = sum(pre_fatal_hi),
        post_fatal    = sum(fatal),
        post_fatal_lo = sum(fatal_lo),
        post_fatal_hi = sum(fatal_hi),
        pre_daly      = sum(pre_daly),
        pre_daly_lo   = sum(pre_daly_low95),
        pre_daly_hi   = sum(pre_daly_hi),
        post_daly     = sum(daly_tot),
        post_daly_lo  = sum(daly_tot_lo),
        post_daly_hi  = sum(daly_tot_hi),
        pre_hosp      = sum(pre_hospitalised),
        pre_hosp_lo   = sum(pre_hospitalised_low95),
        pre_hosp_hi   = sum(pre_hospitalised_hi),
        post_hosp     = sum(hospitalised),
        post_hosp_lo  = sum(hospitalised_lo),
        post_hosp_hi  = sum(hospitalised_hi),
        diff_inf      = sum(diff_inf),
        diff_inf_lo   = sum(diff_inf_lo),
        diff_inf_hi   = sum(diff_inf_hi),
        diff          = sum(diff,          na.rm = TRUE),
        diff_low      = sum(diff_low,      na.rm = TRUE),
        diff_hi       = sum(diff_hi,       na.rm = TRUE),
        diff_fatal    = sum(diff_fatal,    na.rm = TRUE),
        diff_fatal_low= sum(diff_fatal_low,na.rm = TRUE),
        diff_fatal_hi = sum(diff_fatal_hi, na.rm = TRUE),
        diff_daly     = sum(diff_daly,     na.rm = TRUE),
        diff_daly_low = sum(diff_daly_low, na.rm = TRUE),
        diff_daly_hi  = sum(diff_daly_hi,  na.rm = TRUE),
        diff_hosp     = sum(diff_hosp,     na.rm = TRUE),
        diff_hosp_low = sum(diff_hosp_low, na.rm = TRUE),
        diff_hosp_hi  = sum(diff_hosp_hi,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(region = region) %>%
      # NOTE (2026-08): removed /rho_p50 post-scaling + pre_infection/
      # post_infection back-calc -- pre_vacc/post_vacc etc. are already TRUE
      # burden (p_sym-only, no rho) coming out of the m3-refit simulator, so
      # dividing by rho_p50 here would re-inflate by ~1/rho (~7-10x). Same fix
      # as applied to nnv_list() in age_struc_fitting_region_func_updated.R --
      # see chat record, 2026-08.
      mutate(
        scenario_vacc = raw_allocation_age[[idx]]
      ) %>%
      mutate(
        nnv_inf     = scenario_vacc / diff_inf,
        nnv_inf_lo  = scenario_vacc / diff_inf_hi,
        nnv_inf_hi  = scenario_vacc / diff_inf_lo,

        nnv         = scenario_vacc / diff,
        nnv_lo      = scenario_vacc / diff_hi,
        nnv_hi      = scenario_vacc / diff_low,

        nnv_fatal       = scenario_vacc / diff_fatal,
        nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
        nnv_fatal_hi    = scenario_vacc / diff_fatal_low,

        nnv_daly        = scenario_vacc / diff_daly,
        nnv_daly_lo     = scenario_vacc / diff_daly_hi,
        nnv_daly_hi     = scenario_vacc / diff_daly_low,

        nnv_hosp        = scenario_vacc / diff_hosp,
        nnv_hosp_lo     = scenario_vacc / diff_hosp_hi,
        nnv_hosp_hi     = scenario_vacc / diff_hosp_low
      )
    summary_by_age <- summary_by_age %>% mutate(
      target = target
    ) %>%
      relocate(target, .after = scenario)

  }))

  final_summ_df <- final_summ_df %>%
    mutate(
      tot_pop        = rep(N, n_scenarios),
      tot_vacc_prop  = tot_vacc / tot_pop,
      age_gr        = factor(default_age_vector[AgeGroup],
                             levels = gsub("[–—]", "-", default_age_vector))
    )

  final_summ_df$age_gr <- factor(final_summ_df$age_gr, levels = age_gr_levels)
  final_summ_df$region <- region

  per1M_summary <- final_summ_df %>%
    group_by(scenario, target) %>%
    summarise(
      scenario_vacc = first(scenario_vacc),
      tot_pop       = sum(tot_pop),
      vacc_prop     = scenario_vacc / tot_pop,
      diff          = sum(diff,        na.rm = TRUE),
      diff_low      = sum(diff_low,    na.rm = TRUE),
      diff_hi       = sum(diff_hi,     na.rm = TRUE),
      diff_fatal    = sum(diff_fatal,  na.rm = TRUE),
      diff_fatal_low= sum(diff_fatal_low,na.rm = TRUE),
      diff_fatal_hi = sum(diff_fatal_hi, na.rm = TRUE),
      diff_daly     = sum(diff_daly,     na.rm = TRUE),
      diff_daly_low = sum(diff_daly_low, na.rm = TRUE),
      diff_daly_hi  = sum(diff_daly_hi,  na.rm = TRUE),
      diff_hosp     = sum(diff_hosp,     na.rm = TRUE),
      diff_hosp_low = sum(diff_hosp_low, na.rm = TRUE),
      diff_hosp_hi  = sum(diff_hosp_hi,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      per1M_diff       = diff        / scenario_vacc * 1e6,
      per1M_diff_lo    = diff_low    / scenario_vacc * 1e6,
      per1M_diff_hi    = diff_hi     / scenario_vacc * 1e6,

      per1M_fatal      = diff_fatal      / scenario_vacc * 1e6,
      per1M_fatal_lo   = diff_fatal_low  / scenario_vacc * 1e6,
      per1M_fatal_hi   = diff_fatal_hi   / scenario_vacc * 1e6,

      per1M_daly       = diff_daly       / scenario_vacc * 1e6,
      per1M_daly_lo    = diff_daly_low   / scenario_vacc * 1e6,
      per1M_daly_hi    = diff_daly_hi    / scenario_vacc * 1e6,

      per1M_hosp       = diff_hosp       / scenario_vacc * 1e6,
      per1M_hosp_lo    = diff_hosp_low   / scenario_vacc * 1e6,
      per1M_hosp_hi    = diff_hosp_hi    / scenario_vacc * 1e6,

      nnv         = scenario_vacc / diff,
      nnv_lo      = scenario_vacc / diff_hi,
      nnv_hi      = scenario_vacc / diff_low,

      nnv_fatal       = scenario_vacc / diff_fatal,
      nnv_fatal_lo    = scenario_vacc / diff_fatal_hi,
      nnv_fatal_hi    = scenario_vacc / diff_fatal_low,

      nnv_daly        = scenario_vacc / diff_daly,
      nnv_daly_lo     = scenario_vacc / diff_daly_hi,
      nnv_daly_hi     = scenario_vacc / diff_daly_low,

      nnv_hosp        = scenario_vacc / diff_hosp,
      nnv_hosp_lo     = scenario_vacc / diff_hosp_hi,
      nnv_hosp_hi     = scenario_vacc / diff_hosp_low
    )

  list(
    raw_allocation     = raw_allocation,
    raw_allocation_age = raw_allocation_age,
    final_summ_df      = final_summ_df,
    per1M_summary      = per1M_summary
  )
}


## ---- Driver: loop over week files, compute NNV, stack into one table ----

all_weeks_nnv <- purrr::map_dfr(WEEKS, function(week_val) {

  message(sprintf("[finite_weeksweep_nnv] week = %d", week_val))

  load(sprintf("00_Data/0_2_Processed/postsim_vc_ixchiq_model_finite_week%d%s.RData", week_val, FILE_SUFFIX))
  # -> postsim_vc_ixchiq_model (region -> VE -> coverage -> postsim_all_ui output)

  nnv_results_this_week <- purrr::imap(postsim_vc_ixchiq_model, function(ve_cov_list, region_name) {
    purrr::imap(ve_cov_list, function(cov_list, ve_tag) {
      purrr::imap(cov_list, function(postsim_ui, cov_tag) {
        nnv_list(
          vacc_allocation = vacc_allocation(postsim_ui, observed_by_region[[region_name]], region_name),
          postsim_all_ui  = postsim_ui,
          N               = N_by_region[[region_name]],
          region          = region_name,
          observed        = observed_by_region[[region_name]]
        )
      })
    })
  })

  # Save this week's NNV objects too, matching the per-week file pattern
  # already established by finite_weeksweep_fulldraws.R.
  save(nnv_results_this_week,
       file = sprintf("00_Data/0_2_Processed/nnv_results_finite_week%d%s.RData", week_val, FILE_SUFFIX))

  combined_nnv_this_week <- purrr::imap_dfr(nnv_results_this_week, function(ve_cov_list, region_name) {
    purrr::imap_dfr(ve_cov_list, function(cov_list, ve_tag) {
      purrr::imap_dfr(cov_list, function(nnv, cov_tag) {
        nnv$final_summ_df %>%
          mutate(
            region  = region_name,
            VE      = ve_tag,
            VC      = cov_tag,
            setting = setting_key[region_name],
            week    = week_val
          )
      })
    })
  })

  combined_nnv_this_week
})

save(all_weeks_nnv, file = sprintf("00_Data/0_2_Processed/combined_nnv_df_region_coverage_model_finite_weeksweep%s.RData", FILE_SUFFIX))

dir.create("02_Outputs/2_2_Tables", showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(all_weeks_nnv, path = sprintf("02_Outputs/2_2_Tables/combined_nnv_df_region_coverage_model_finite_weeksweep%s.xlsx", FILE_SUFFIX))

message("[finite_weeksweep_nnv] done. Rows: ", nrow(all_weeks_nnv))
message("Columns include: region, VE, VC (coverage), week, scenario, target, diff/diff_fatal/diff_daly/diff_hosp (averted), nnv/nnv_fatal/nnv_daly/nnv_hosp, tot_vacc, ...")
