# =============================================================================
# 29_seromix_within_campaign_dynamics.R
#
# Within-campaign seromix dynamics: fixes the campaign start at week 2 (the
# reference outbreak-response scenario) and tracks how the CUMULATIVE
# seropositive share among vaccinees evolves over the following weeks of that
# SAME campaign's rollout -- not a comparison across different campaign start
# weeks (that was 28_seromix_dynamics_by_week.R, which conflated "start timing"
# scenarios with "time since start" dynamics; this script fixes that distinction
# per chat record 2026-08-31).
#
# Uses 01_Data/postsim_vc_ixchiq_model_finite.RData -- the reference scenario's
# RAW simulation output (week=2, coverage=50% only; this file predates the
# week-sweep and never collapses its time dimension). raw_allocation_array /
# vacc_to_S_array are [age_index=20, week=52, draw=1000] per region x VE x
# scenario(age-target). Cumulative-sum over the week dimension (not a single
# end-of-campaign total, unlike calc_q_seromix_for_scenarios_agecat() elsewhere)
# gives the running "of everyone vaccinated SO FAR (weeks 2..W), what fraction
# were seronegative at the moment of their dose" -- genuine within-campaign
# dynamics for this one fixed start week.
#
# Output: 06_Results/seromix_within_campaign_dynamics.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

load("01_Data/setting_key.RData")  # -> setting_key (state -> Low/Moderate/High)

cat("Loading reference postsim (~2.9GB, one-time)...\n")
load("01_Data/postsim_vc_ixchiq_model_finite.RData")
cat("Loaded.\n")

regions <- names(postsim_vc_ixchiq_model)
setting_levels <- c("Low", "Moderate", "High")

extract_cumulative <- function(scenario_int, age_label) {
  purrr::map_dfr(regions, function(reg) {
    reg_list <- postsim_vc_ixchiq_model[[reg]]
    purrr::imap_dfr(reg_list, function(ve_list, ve_name) {
      cov_list <- ve_list[["cov50"]]
      sc <- cov_list$scenario_result[[scenario_int]]
      raw_alloc <- sc$sim_result$raw_allocation_array   # [age_index, week, draw]
      vacc_to_S <- sc$sim_result$vacc_to_S_array

      total_by_week_draw   <- apply(raw_alloc, c(2, 3), sum, na.rm = TRUE)  # [week x draw]
      seroneg_by_week_draw <- apply(vacc_to_S, c(2, 3), sum, na.rm = TRUE)

      cum_total   <- apply(total_by_week_draw,   2, cumsum)  # cumulative over week, per draw
      cum_seroneg <- apply(seroneg_by_week_draw, 2, cumsum)

      q_seroneg_cum <- cum_seroneg / cum_total
      med_by_week <- apply(q_seroneg_cum, 1, median, na.rm = TRUE)

      tibble(Region = reg, VE = ve_name, week = seq_len(nrow(cum_total)),
             q_seroneg_cum_med = med_by_week)
    })
  }) %>% mutate(AgeCat = age_label)
}

# Scenario 3 = 18-64 target campaign, Scenario 4 = 65+ target campaign
# (map_scenario_agecat_int convention used throughout this pipeline).
plot_df <- bind_rows(
  extract_cumulative(3, "18-64"),
  extract_cumulative(4, "65+")
) %>%
  mutate(
    setting = unname(setting_key[Region]),
    setting = factor(setting, levels = setting_levels),
    AgeCat  = factor(AgeCat, levels = c("18-64", "65+"), labels = c("18-64 years", "65+ years")),
    Mechanism = factor(VE, levels = c("VE0", "VE98.9"),
                        labels = c("Disease blocking only", "Disease and infection blocking")),
    seropos_pct = 100 * (1 - q_seroneg_cum_med)
  ) %>%
  filter(!is.na(setting), week >= 2) %>%  # week 1 is pre-campaign-start (NA/undefined, 0 doses yet)
  group_by(setting, AgeCat, Mechanism, week) %>%
  summarise(seropos_pct = median(seropos_pct, na.rm = TRUE), .groups = "drop")

setting_colors <- c("Low" = "#1b9e77", "Moderate" = "#d95f02", "High" = "#d7191c")

theme_seromix2 <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", colour = NA),
      strip.text = element_text(face = "bold", size = rel(1.1)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = rel(1.2)),
      plot.subtitle = element_text(colour = "grey30")
    )
}

p <- ggplot(plot_df, aes(x = week, y = seropos_pct, colour = setting, group = setting)) +
  geom_line(linewidth = 1) +
  facet_grid(AgeCat ~ Mechanism) +
  scale_colour_manual(values = setting_colors, name = "Transmission setting") +
  scale_x_continuous(name = "Calendar week (campaign started week 2, 50% coverage target)",
                      breaks = seq(2, 52, by = 10)) +
  scale_y_continuous(name = "Cumulative seropositive among vaccinees so far (%)",
                      labels = label_number(suffix = "%")) +
  labs(
    title = "Within-campaign seromix dynamics (fixed start: week 2)",
    subtitle = "One campaign, tracked forward -- not a comparison across different start weeks"
  ) +
  theme_seromix2()

ggsave("06_Results/seromix_within_campaign_dynamics.png", p, width = 10, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/seromix_within_campaign_dynamics.png")
