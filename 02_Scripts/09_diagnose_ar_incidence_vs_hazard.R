# =============================================================================
# 09_diagnose_ar_incidence_vs_hazard.R
#
# Diagnostic requested to test whether the travel-vaccination FOI rescaling
# step (06_brazil_travel_final_finite.R) is a genuine correction or a
# self-calibration artifact.
#
# For each state's Stan posterior (fit_prevacc_<state>_finite, saved fields in
# posterior_finite_all.RData), and for each of the 1000 posterior draws, we
# compare two attack-rate estimates computed WITHIN THE SAME FITTED EPIDEMIC
# TRAJECTORY (no downstream resimulation, no rho-correction, no bootstrap):
#
#   AR_incidence = sum_t(new_exposed_t) / S_start
#                = sum_t(incident_infections[t]) / (susceptible_fraction[1] * pop_total)
#     -- cumulative new infections (S->E transitions) over the observed
#     period, divided by the susceptible population at the start of the
#     observed period.
#
#   AR_hazard = 1 - exp(-sum_t(lambda_t))
#             = 1 - exp(-sum_t(phi_pred[t]))
#     -- the hazard-based attack rate implied by the same model's own FOI
#     trajectory over the observed period.
#
# If these two are already close (same fitted epidemic, so they should be,
# modulo discretisation/susceptible-depletion effects), then the travel
# script's rescaling of FOI to match a SEPARATELY-computed target attack rate
# (tot_inf/total_S0, involving rho-correction and a bootstrapped S0) is not
# correcting a flaw in the fitted model itself -- any gap between raw FOI and
# the travel target must be coming from that separate downstream pipeline,
# not from an AR_incidence vs AR_hazard mismatch within the Stan fit.
#
# Output: printed summary table (median, 95% draw-level range of the ratio
# AR_incidence/AR_hazard, per state) -- no files saved, diagnostic only.
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages(library(dplyr))

load("01_Data/posterior_finite_all.RData")

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

diag_one_state <- function(post, state_name) {
  n_draws <- nrow(post$incident_infections)

  # Susceptible population at the START of the observed period (not the
  # unobserved burn-in start): susceptible_fraction[,1] is already indexed to
  # observed week 1 (model week B+1) in the Stan generated quantities block.
  S_start <- post$susceptible_fraction[, 1] * post$pop_total

  AR_incidence <- rowSums(post$incident_infections) / S_start
  AR_hazard    <- 1 - exp(-rowSums(post$phi_pred))

  data.frame(
    state         = state_name,
    draw          = seq_len(n_draws),
    S_start       = S_start,
    total_new_exp = rowSums(post$incident_infections),
    sum_lambda    = rowSums(post$phi_pred),
    AR_incidence  = AR_incidence,
    AR_hazard     = AR_hazard,
    ratio         = AR_incidence / AR_hazard,
    diff          = AR_incidence - AR_hazard
  )
}

all_diag <- do.call(rbind, Map(diag_one_state, state_list, names(state_list)))

summ <- all_diag %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    AR_incidence_med = median(AR_incidence),
    AR_hazard_med    = median(AR_hazard),
    ratio_med        = median(ratio),
    ratio_lo95       = quantile(ratio, 0.025),
    ratio_hi95       = quantile(ratio, 0.975),
    max_abs_diff     = max(abs(diff)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(abs(ratio_med - 1)))

cat("=== AR_incidence (new_exposed-based) vs AR_hazard (FOI-based), same Stan posterior ===\n")
cat("Ratio near 1.00 across the board => the two are consistent within the fitted\n")
cat("model itself; any FOI/target-AR gap used in travel rescaling would then come\n")
cat("from the SEPARATE downstream pipeline (rho-correction / bootstrapped S0),\n")
cat("not from a flaw in the fitted model.\n\n")
print(as.data.frame(summ), digits = 4)

cat("\n=== Draw-level distribution of ratio, all states pooled ===\n")
print(summary(all_diag$ratio))

# NOTE (2026-08): the AR_lin_S0-vs-AR_hazard comparison that used to live here
# (loading 01_Data/ar_table_by_state_finite.RData, produced by the now-deleted
# 02b_setup_ar_by_state.R) is obsolete -- 06_brazil_travel_final_finite.R has
# since been rewritten to use AR_hazard (phi_pred) directly, so there is no
# more "current target" to compare against. See chat record, 2026-08.
