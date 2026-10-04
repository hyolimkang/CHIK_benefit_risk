# =============================================================================
# 22_export_infection_risk_full_grid.R
#
# Full simulated infection-risk grid, exported as-is: Setting x Age group x
# Duration x Entry timing. P(infection) = 1 - exp(-H_trip), age-independent
# by construction (H_trip/AR_travel computed from state-level FOI alone in
# 06_brazil_travel_final_finite.R, before age-specific SAE/death rates are
# applied -- verified empirically: AR_median is bit-identical across
# age_group for every setting/duration/entry_day cell) -- kept as an explicit
# column anyway per request, not collapsed.
#
# Alongside P(infection), each row also carries the corresponding Benefit,
# Risk, BRR, and Pr(BRR>1) for all three outcomes (SAE, Death, DALY) at that
# exact same setting/age/duration/entry_day cell -- these ARE age-specific
# (unlike infection risk), since they depend on age-specific SAE/death rates
# applied downstream of the shared infection-risk calculation.
#
# 95% UI (2.5%/97.5%) now included for every median column (chat record
# 2026-08-27): 06_brazil_travel_final_finite.R's summarise_grid_by_setting()/
# median_only_df() were patched to compute quantile(0.025/0.975) alongside
# median() -- the raw per-draw array was already fully in memory at that
# collapse step, so this added negligible cost, no new simulation needed.
#
# Duration/entry-day combinations that would run past the observed 2022 FOI
# series (entry_day + duration > 364) are genuinely not computable and are
# left as NA -- not extrapolated or dropped (see figure2a_travel_infection_
# risk_surface.png's investigation of this exact grid).
#
# Output: 06_Results/infection_risk_full_grid.xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(writexl)
})

load("01_Data/psa_grid_bra_travel_finite_setting.RData")  # -> psa_grid_travel_setting

full_grid <- psa_grid_travel_setting %>%
  mutate(
    entry_week = (entry_day - 1) %/% 7 + 1,
    setting    = factor(setting, levels = c("Low", "Moderate", "High")),
    age_group  = factor(age_group, levels = c("18-64", "65+"))
  ) %>%
  arrange(setting, age_group, duration, entry_day) %>%
  transmute(
    Setting       = setting,
    `Age group`   = age_group,
    `Duration (days)` = duration,
    `Entry day (day of year)` = entry_day,
    `Entry week`  = entry_week,
    `P(infection), median` = AR_median,
    `P(infection), 2.5%`   = AR_lo,
    `P(infection), 97.5%`  = AR_hi,

    `SAE: Benefit per 10,000, median` = averted_sae_median,
    `SAE: Benefit per 10,000, 2.5%`   = averted_sae_lo,
    `SAE: Benefit per 10,000, 97.5%`  = averted_sae_hi,
    `SAE: Risk per 10,000, median`    = excess_sae_median,
    `SAE: Risk per 10,000, 2.5%`      = excess_sae_lo,
    `SAE: Risk per 10,000, 97.5%`     = excess_sae_hi,
    `SAE: BRR, median` = brr_sae_median,
    `SAE: BRR, 2.5%`   = brr_sae_lo,
    `SAE: BRR, 97.5%`  = brr_sae_hi,
    `SAE: Pr(BRR>1)`   = brr_sae_prob_gt1,

    `Death: Benefit per 10,000, median` = averted_death_median,
    `Death: Benefit per 10,000, 2.5%`   = averted_death_lo,
    `Death: Benefit per 10,000, 97.5%`  = averted_death_hi,
    `Death: Risk per 10,000, median`    = excess_death_median,
    `Death: Risk per 10,000, 2.5%`      = excess_death_lo,
    `Death: Risk per 10,000, 97.5%`     = excess_death_hi,
    `Death: BRR, median` = brr_death_median,
    `Death: BRR, 2.5%`   = brr_death_lo,
    `Death: BRR, 97.5%`  = brr_death_hi,
    `Death: Pr(BRR>1)`   = brr_death_prob_gt1,

    `DALY: Benefit per 10,000, median` = daly_averted_median,
    `DALY: Benefit per 10,000, 2.5%`   = daly_averted_lo,
    `DALY: Benefit per 10,000, 97.5%`  = daly_averted_hi,
    `DALY: Risk per 10,000, median`    = daly_sae_median,
    `DALY: Risk per 10,000, 2.5%`      = daly_sae_lo,
    `DALY: Risk per 10,000, 97.5%`     = daly_sae_hi,
    `DALY: BRR, median` = brr_daly_median,
    `DALY: BRR, 2.5%`   = brr_daly_lo,
    `DALY: BRR, 97.5%`  = brr_daly_hi,
    `DALY: Pr(BRR>1)`   = brr_daly_prob_gt1
  )

cat("Rows:", nrow(full_grid), "\n")
cat("NA rows (not computable -- entry_day + duration > 364):", sum(is.na(full_grid$`P(infection), median`)), "\n")

writexl::write_xlsx(full_grid, "06_Results/infection_risk_full_grid.xlsx")
message("Saved: 06_Results/infection_risk_full_grid.xlsx")
