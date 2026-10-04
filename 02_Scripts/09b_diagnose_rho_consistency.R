# =============================================================================
# 09b_diagnose_rho_consistency.R
#
# For each state x posterior draw, computes:
#   AR_Stan          = sum_t(new_exposed_t) / S0                 (S0 = susceptible
#                      population at start of observed period)
#   rho_post         = the Stan posterior draw's own estimated rho
#   rho_sym_implied  = reported_cases / (sum_t(new_exposed_t) * p_sym)
#                      -- an empirically-backed-out reporting rate, assuming a
#                      fixed external symptomatic proportion (p_sym) and total
#                      TRUE infections (new_exposed, all ages/weeks) from the
#                      same posterior draw, compared against the actual raw
#                      reported case total fed to the Stan fit.
#
# reported_cases is fixed data (same across draws, per state) -- total
# observed cases summed over all 52 weeks, from observed_2022.RData (the same
# raw weekly case data used earlier this session for the setting-key fix).
#
# p_sym = 0.524 (symp_const_fixed, as used in 02b_setup_ar_by_state.R).
#
# Purpose: sanity-check whether the Stan model's own internally-estimated rho
# (which the model docs note combines symptomatic probability AND reporting
# probability) is consistent with an independently-backed-out estimate using
# a fixed external p_sym -- per state, median + 95% UI.
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages(library(dplyr))

load("01_Data/posterior_finite_all.RData")
load("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_2_Processed/observed_2022.RData")

state_list <- list(
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

observed_list <- list(
  "Ceará"               = observed_ce,
  "Bahia"               = observed_bh,
  "Paraíba"             = observed_pa,
  "Pernambuco"          = observed_pn,
  "Rio Grande do Norte" = observed_rg,
  "Piauí"               = observed_pi,
  "Alagoas"             = observed_ag,
  "Tocantins"           = observed_tc,
  "Minas Gerais"        = observed_mg,
  "Sergipe"             = observed_se,
  "Goiás"               = observed_go
)

p_sym <- 0.524  # symp_const_fixed, 02b_setup_ar_by_state.R

diag_one <- function(post, obs_df, state_name) {
  n_draws <- nrow(post$incident_infections)
  S0 <- post$susceptible_fraction[, 1] * post$pop_total
  total_new_exposed <- rowSums(post$incident_infections)
  reported_total <- sum(obs_df$Observed, na.rm = TRUE)  # fixed data, same for every draw

  data.frame(
    state            = state_name,
    draw             = seq_len(n_draws),
    reported_total   = reported_total,
    total_new_exposed = total_new_exposed,
    AR_Stan          = total_new_exposed / S0,
    rho_post         = post$rho,
    rho_sym_implied  = reported_total / (total_new_exposed * p_sym)
  )
}

all_diag <- do.call(rbind, Map(diag_one, state_list, observed_list, names(state_list)))

summ <- all_diag %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    reported_total       = dplyr::first(reported_total),
    AR_Stan_med          = median(AR_Stan),
    AR_Stan_lo           = quantile(AR_Stan, 0.025),
    AR_Stan_hi           = quantile(AR_Stan, 0.975),
    rho_post_med         = median(rho_post),
    rho_post_lo          = quantile(rho_post, 0.025),
    rho_post_hi          = quantile(rho_post, 0.975),
    rho_sym_implied_med  = median(rho_sym_implied),
    rho_sym_implied_lo   = quantile(rho_sym_implied, 0.025),
    rho_sym_implied_hi   = quantile(rho_sym_implied, 0.975),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(AR_Stan_med))

cat("=== AR_Stan, rho_post (Stan-internal), rho_sym_implied (data-backed-out, p_sym=0.524) ===\n\n")
print(as.data.frame(summ), digits = 3)

cat("\n=== rho_post vs rho_sym_implied ratio (median) ===\n")
summ2 <- summ %>% dplyr::mutate(ratio_rho = rho_post_med / rho_sym_implied_med) %>%
  dplyr::select(state, rho_post_med, rho_sym_implied_med, ratio_rho)
print(as.data.frame(summ2), digits = 3)

cat("\n=== Structural check: infections > symptomatic > reported? ===\n")
summ3 <- all_diag %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    total_infections_med  = median(total_new_exposed),
    implied_symptomatic_med = median(total_new_exposed * 0.524),
    reported_total = dplyr::first(reported_total),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    infections_gt_symptomatic = total_infections_med > implied_symptomatic_med,
    symptomatic_gt_reported   = implied_symptomatic_med > reported_total
  ) %>%
  dplyr::arrange(dplyr::desc(total_infections_med))
print(as.data.frame(summ3), digits = 4)
