# =============================================================================
# 28_seromix_dynamics_by_week.R
#
# Line graph: FINAL (end-of-campaign) seropositive share among vaccinees, for
# each of 7 INDEPENDENT campaign-start-week scenarios (week 1/8/16/24/32/42;
# week 52 excluded -- see below), faceted by Mechanism (DB / D+I) x Age group
# (18-64 / 65+), coloured by Coverage target (10/30/60/90%). Setting fixed to
# High (the three settings track each other closely -- see chat record; not
# worth quadrupling the facet count over a dimension that barely moves this).
#
# IMPORTANT (chat record 2026-08-31): this is a comparison ACROSS separate
# scenarios that each start at a different week, NOT a within-campaign time
# series -- each x-value is its own independent simulation run. For how ONE
# fixed-start campaign's seromix evolves forward in time, see
# 29_seromix_within_campaign_dynamics.R instead; conflating the two was the
# original version of this script's mistake.
#
# week=52 is DROPPED: starting a campaign that late leaves only ~6 of the
# observed-FOI-year's 364 days remaining, so the simulation truncates the
# rollout well short of the target coverage (verified: median total_vacc_age
# collapses to ~1/10th of every other week's value at week=52) -- its "final"
# composition reflects an incomplete, arbitrarily-cut-short cohort, not a
# genuine campaign outcome, so it isn't comparable to the other 6 points.
#
# q_seroneg_vacc isn't saved directly in the week-sweep per-draw export, but
# sae_10k_seroneg = 1e4 * q_seroneg_vacc * p_sae_vacc_base (death=0 for 18-64,
# and the same identity holds for 65+ since RR_seropos=1 defines the seroneg
# contribution as the unweighted base rate) while sae_10k_base = 1e4 *
# p_sae_vacc_base unconditionally -- so their ratio recovers q_seroneg_vacc
# exactly without needing the intermediate q_all_regions object.
#
# Inputs (CHIK_ORV_impact project -- the week-sweep per-draw RData live there,
# not in this repo's 01_Data/):
#   - draw_level_xy_serostatus_finite_weeksweep.RData      (D+I / VE98.9)
#   - draw_level_xy_serostatus_finite_weeksweep_ve0.RData  (DB / VE0)
#
# Output: 06_Results/seromix_dynamics_by_week.png
# =============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

ori_root <- "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"

extract_seromix <- function(rdata_path, mechanism_label) {
  e <- new.env()
  load(rdata_path, envir = e)
  e$all_weeks_brr %>%
    filter(outcome == "DALY", RR_seropos == 0, setting == "High",
           AgeCat %in% c("18-64", "65+"),
           week != 52) %>%   # truncated rollout -- not a comparable endpoint
    mutate(
      q_seroneg = sae_10k_seroneg / sae_10k_base,
      seropos_pct = 100 * (1 - q_seroneg)
    ) %>%
    group_by(week, Coverage, AgeCat) %>%
    summarise(seropos_med = median(seropos_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(Mechanism = mechanism_label)
}

plot_df <- bind_rows(
  extract_seromix(file.path(ori_root, "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep_ve0.RData"), "Disease blocking only"),
  extract_seromix(file.path(ori_root, "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep.RData"), "Disease and infection blocking")
) %>%
  mutate(
    Mechanism = factor(Mechanism, levels = c("Disease blocking only", "Disease and infection blocking")),
    AgeCat    = factor(AgeCat, levels = c("18-64", "65+"), labels = c("18-64 years", "65+ years")),
    coverage_pct = as.integer(sub("cov", "", Coverage)),
    Coverage  = factor(paste0(coverage_pct, "%"), levels = paste0(sort(unique(coverage_pct)), "%"))
  )

# Coverage is an ordered quantity, not an unordered category -- single-hue
# light-to-dark ramp (not a 4-hue categorical set), matching the convention
# already used for BRR heatmaps in this project (continuous/ordered ->
# sequential ramp, never discrete unrelated hues).
coverage_ramp <- c("10%" = "#c6dbef", "30%" = "#6baed6", "60%" = "#2171b5", "90%" = "#08306b")

theme_seromix <- function(base_size = 14) {
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

p <- ggplot(plot_df, aes(x = week, y = seropos_med, colour = Coverage, group = Coverage)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(AgeCat ~ Mechanism) +
  scale_colour_manual(values = coverage_ramp, name = "Coverage target") +
  scale_x_continuous(name = "Campaign start week (each point = a separate scenario)", breaks = c(1, 8, 16, 24, 32, 42)) +
  scale_y_continuous(name = "Final seropositive among vaccinees (%)", labels = label_number(suffix = "%")) +
  labs(
    title = "End-of-campaign serostatus composition, across independent start-week scenarios",
    subtitle = "High-transmission setting; week 52 excluded (rollout truncated by the observed-FOI-year boundary, not comparable)"
  ) +
  theme_seromix()

ggsave("06_Results/seromix_dynamics_by_week.png", p, width = 10, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/seromix_dynamics_by_week.png")

print(as.data.frame(plot_df %>% arrange(Mechanism, AgeCat, Coverage, week)))
