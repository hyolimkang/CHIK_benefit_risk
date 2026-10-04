# =============================================================================
# 21_foi_risk_brr_peak_timing.R
#
# Quantifies whether the week of peak daily FOI (per setting) lines up with
# the week of peak travel infection risk and the week of peak Pr(BRR>1) --
# i.e. does the outbreak's actual timing drive both the risk surface and the
# benefit-risk surface peaks in the same direction, at the same time.
#
# FOI: state-level posterior weekly phi_pred, states pooled WITH their draws
# within each setting (same convention as psa_grid_travel_setting's own
# setting-level pooling -- not a simple average of per-state medians), then
# converted to daily via /7 (same as weekly_to_daily_foi() elsewhere).
#
# Infection risk / Pr(BRR>1): from the fine entry_day x duration grid
# (psa_grid_bra_travel_finite_setting.RData), at a representative 30-day
# duration (age-independent for infection risk; DALY outcome for BRR),
# restricted to entry week <=40 -- the range verified elsewhere in this
# session to be free of the boundary extrapolation issue.
#
# Output: 06_Results/foi_risk_brr_peak_timing.xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(writexl)
})

setting_levels <- c("Low", "Moderate", "High")

state_abbrev_map <- c(
  ce = "Ceará", bh = "Bahia", pa = "Paraíba", pn = "Pernambuco",
  rg = "Rio Grande do Norte", pi = "Piauí", ag = "Alagoas",
  tc = "Tocantins", mg = "Minas Gerais", se = "Sergipe", go = "Goiás"
)

load("01_Data/setting_key.RData")  # -> setting_key

# ---- A) Weekly FOI profile, pooled with draws within each setting ----------
foi_cache_file <- "01_Data/setting_weekly_foi_pooled_finite.RData"

if (!file.exists(foi_cache_file)) {
  message("Cache not found -- loading 01_Data/posterior_finite_all.RData (~1.3GB, slow)...")
  load("01_Data/posterior_finite_all.RData")

  weekly_foi_by_setting <- lapply(setting_levels, function(s) {
    states_in_setting <- names(setting_key)[setting_key == s]
    abbrevs <- names(state_abbrev_map)[state_abbrev_map %in% states_in_setting]
    mats <- lapply(abbrevs, function(ab) get(paste0("posterior_", ab))$phi_pred)  # each 1000 x 52
    pooled <- do.call(rbind, mats)  # (1000*n_states) x 52 -- draws AND states pooled together
    apply(pooled, 2, median)  # weekly median FOI, pooled
  })
  names(weekly_foi_by_setting) <- setting_levels

  save(weekly_foi_by_setting, file = foi_cache_file)
  message("Saved: ", foi_cache_file)
} else {
  load(foi_cache_file)  # -> weekly_foi_by_setting
}

foi_peak <- bind_rows(lapply(setting_levels, function(s) {
  wk <- weekly_foi_by_setting[[s]]
  data.frame(
    setting = s,
    foi_peak_week = which.max(wk),
    foi_peak_daily = max(wk) / 7  # weekly -> daily, same convention as weekly_to_daily_foi()
  )
}))

# ---- B) Infection risk + Pr(BRR>1) peak timing, 30-day reference duration --
load("01_Data/psa_grid_bra_travel_finite_setting.RData")  # -> psa_grid_travel_setting

grid_30d <- psa_grid_travel_setting %>%
  mutate(entry_week = (entry_day - 1) %/% 7 + 1) %>%
  filter(duration == 30, entry_week <= 40, age_group == "18-64")

risk_peak <- grid_30d %>%
  group_by(setting) %>%
  slice_max(AR_median, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(setting, risk_peak_week = entry_week, risk_peak_pct = 100 * AR_median)

brr_peak <- grid_30d %>%
  group_by(setting) %>%
  slice_max(brr_daly_prob_gt1, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(setting, brr_peak_week = entry_week, brr_peak_prob_pct = 100 * brr_daly_prob_gt1)

# ---- C) Assemble + alignment check ------------------------------------------
summary_tbl <- foi_peak %>%
  left_join(risk_peak, by = "setting") %>%
  left_join(brr_peak, by = "setting") %>%
  mutate(
    setting = factor(setting, levels = setting_levels),
    weeks_foi_to_risk = risk_peak_week - foi_peak_week,
    weeks_foi_to_brr  = brr_peak_week - foi_peak_week
  ) %>%
  arrange(setting) %>%
  transmute(
    Setting = setting,
    `Peak FOI week` = foi_peak_week,
    `Peak daily FOI` = signif(foi_peak_daily, 3),
    `Peak infection-risk week (30d trip)` = risk_peak_week,
    `Peak infection risk (%)` = round(risk_peak_pct, 3),
    `Peak Pr(BRR>1) week (30d trip)` = brr_peak_week,
    `Peak Pr(BRR>1) (%)` = round(brr_peak_prob_pct, 1),
    `Weeks: FOI peak -> risk peak` = weeks_foi_to_risk,
    `Weeks: FOI peak -> BRR-prob peak` = weeks_foi_to_brr
  )

print(as.data.frame(summary_tbl))

writexl::write_xlsx(summary_tbl, "06_Results/foi_risk_brr_peak_timing.xlsx")
message("Saved: 06_Results/foi_risk_brr_peak_timing.xlsx")
