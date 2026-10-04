# =============================================================================
# 23_travel_brr_range_by_timing.R
#
# For each representative travel duration (7/14/30/90 days), shows how much
# DALY BRR varies purely as a function of entry timing (calendar week of
# travel): a vertical range line from the worst-case (off-season) entry week
# to the best-case (peak-season) entry week, with the pooled/"average traveler"
# median BRR (from the headline discrete-PSA table, which marginalises over a
# uniformly random entry day across the year) overlaid as a reference point.
#
# Motivation (chat record 2026-08-28): headline_summary_travel_*.docx reports
# a single pooled BRR per duration (e.g. DALY, High, 18-64, 7d = 1.49), which
# looks low next to infection_risk_full_grid.xlsx's per-entry-week values
# (same cell can be > 9 at the seasonal peak). Both are correct -- the pooled
# number averages over ~52 weeks of which only ~20 have BRR > 1 -- but the
# headline table alone doesn't show that spread. This figure makes the spread
# explicit.
#
# Inputs:
#   - 01_Data/psa_grid_bra_travel_finite_setting.RData (-> psa_grid_travel_setting)
#     fine entry_day x duration grid, pooled by setting -- gives the min/max
#     across entry weeks for each setting x age x duration.
#   - 01_Data/psa_df_bra_travel_finite.RData (-> psa_df) + 01_Data/setting_key.RData
#     discrete PSA sim (uniform-random entry day) -- gives the pooled/"average
#     traveler" reference BRR, computed identically to
#     11_travel_headline_figures.R's build_travel_synthesis_data().
#
# Output: 06_Results/travel_brr_range_by_timing_DALY.png
#         06_Results/travel_brr_range_by_timing_DALY.xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")

library(dplyr)
library(ggplot2)
library(scales)
library(writexl)

# ---- 1) Min/max BRR across entry weeks, from the fine grid -----------------
# Min = off-season entry (worst timing), max = peak-transmission-season entry
# (best timing) -- keep the actual calendar week each occurs at, since the
# peak week differs by setting (Low ~wk10, Moderate ~wk14-16, High ~wk17; see
# 21_foi_risk_brr_peak_timing.R), so labelling by week ties the range
# directly back to the FOI seasonality rather than an arbitrary "best/worst".
load("01_Data/psa_grid_bra_travel_finite_setting.RData")  # -> psa_grid_travel_setting

grid_prepped <- psa_grid_travel_setting %>%
  filter(duration %in% c(7, 14, 30, 90), is.finite(brr_daly_median)) %>%
  mutate(
    entry_week = (entry_day - 1) %/% 7 + 1,
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    age_group = factor(age_group, levels = c("18-64", "65+")),
    days = factor(paste0(duration, "d"), levels = c("7d", "14d", "30d", "90d"))
  )

min_df <- grid_prepped %>%
  group_by(setting, age_group, days) %>%
  slice_min(brr_daly_median, n = 1, with_ties = FALSE) %>%
  transmute(setting, age_group, days, brr_min = brr_daly_median, week_min = entry_week,
            brr_min_lo = brr_daly_lo, brr_min_hi = brr_daly_hi)

max_df <- grid_prepped %>%
  group_by(setting, age_group, days) %>%
  slice_max(brr_daly_median, n = 1, with_ties = FALSE) %>%
  transmute(setting, age_group, days, brr_max = brr_daly_median, week_max = entry_week,
            brr_max_lo = brr_daly_lo, brr_max_hi = brr_daly_hi)

range_df <- min_df %>% left_join(max_df, by = c("setting", "age_group", "days"))

# ---- 2) Pooled "average traveler" BRR, from the discrete PSA sim -----------
# (identical grouping/summary logic to 11_travel_headline_figures.R's
# build_travel_synthesis_data(), restricted to the DALY outcome.)
load("01_Data/psa_df_bra_travel_finite.RData")  # -> psa_df
load("01_Data/setting_key.RData")               # -> setting_key

pooled_df <- psa_df %>%
  mutate(setting = unname(setting_key[state])) %>%
  filter(!is.na(setting)) %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    days = factor(days, levels = c("7d", "14d", "30d", "90d")),
    age_group = factor(age_group, levels = c("18-64", "65+"))
  ) %>%
  group_by(setting, age_group, days) %>%
  summarise(
    brr_pooled = median(brr_daly, na.rm = TRUE),
    brr_pooled_lo = quantile(brr_daly, 0.025, na.rm = TRUE),
    brr_pooled_hi = quantile(brr_daly, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- range_df %>%
  left_join(pooled_df, by = c("setting", "age_group", "days"))

# ---- 3) Range plot -----------------------------------------------------
theme_hm <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", colour = NA),
      strip.text = element_text(face = "bold", size = rel(1.25)),
      axis.title = element_text(face = "bold", size = rel(1.15)),
      axis.text = element_text(size = rel(1.05), colour = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = rel(1.05)),
      plot.title = element_text(face = "bold", size = rel(1.3))
    )
}

# One flat row per duration: only the x-axis (BRR) carries real information,
# so the connecting line is a straight horizontal segment, not an artificial
# diagonal -- its length (on the log-BRR scale) is the only thing that means
# anything: how far worst-case and best-case entry timing pull apart.
days_levels <- c("7d", "14d", "30d", "90d")
plot_df <- plot_df %>%
  mutate(days = factor(days, levels = days_levels))

p <- ggplot(plot_df, aes(y = days)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = brr_min, xmax = brr_max), colour = "grey40", linewidth = 0.9) +
  geom_point(aes(x = brr_min, shape = "Off-season entry (worst week)"),
             size = 3.6, colour = "#4575B4") +
  geom_point(aes(x = brr_max, shape = "Peak-transmission-season entry (best week)"),
             size = 3.6, colour = "#D73027") +
  geom_point(aes(x = brr_pooled, shape = "Pooled"),
             size = 4.2, colour = "black") +
  scale_shape_manual(
    name = NULL,
    values = c(
      "Off-season entry (worst week)" = 16,
      "Peak-transmission-season entry (best week)" = 16,
      "Pooled" = 18
    )
  ) +
  scale_x_log10(
    name = "DALY BRR (log scale)",
    breaks = scales::breaks_log(n = 6),
    labels = scales::label_number(accuracy = NULL, drop0trailing = TRUE)
  ) +
  scale_y_discrete(name = "Representative travel duration", limits = rev(days_levels)) +
  facet_grid(age_group ~ setting) +
  labs(title = "C. DALY BRR range by entry timing") +
  theme_hm()

ggsave("06_Results/travel_brr_range_by_timing_DALY.png", p, width = 13, height = 8.5, dpi = 300, bg = "white")
message("Saved: 06_Results/travel_brr_range_by_timing_DALY.png")

export_df <- plot_df %>%
  select(setting, age_group, days,
         week_min, brr_min, brr_min_lo, brr_min_hi,
         week_max, brr_max, brr_max_lo, brr_max_hi,
         brr_pooled, brr_pooled_lo, brr_pooled_hi)
write_xlsx(export_df, "06_Results/travel_brr_range_by_timing_DALY.xlsx")
message("Saved: 06_Results/travel_brr_range_by_timing_DALY.xlsx")

print(as.data.frame(export_df))
