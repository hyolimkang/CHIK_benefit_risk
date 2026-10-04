# =============================================================================
# 30_cumulative_attack_rate_dynamics.R
#
# Population-wide (NOT vaccination-status-split) cumulative attack rate over
# calendar week, for the reference scenario (campaign start week 2, coverage
# 50%). No simulator code changes needed -- age_array_raw_inf (raw new
# infections per age bin per week) is already saved in
# postsim_vc_ixchiq_model_finite.RData, so cumulative-summing it over the week
# dimension directly gives the population's cumulative-infected trajectory
# (approximately the R compartment size over time, since reinfection within
# one year is not modelled here) -- this is a DIFFERENT quantity from
# "seropositive share AMONG VACCINEES" (scripts 28/29): it's the whole
# population's immunity buildup, unsplit by vaccination status. Chat record
# 2026-08-31: getting the vaccination-status-split version (which of the
# vaccinated-while-susceptible cohort later got infected) needs the simulator
# itself modified to track that cohort -- not attempted here.
#
# Output: 06_Results/cumulative_attack_rate_dynamics.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

load("01_Data/setting_key.RData")
setting_levels <- c("Low", "Moderate", "High")

# age_index 1:20 -> age_gr per age_struc_fitting's age_map (bins 1-10 = u40,
# bins 11-20 = o40); 18-64 spans bins 6-15 (20-24 through 60-64), 65+ spans
# bins 16-20 (65-69 through 85+) per the age_map used throughout this pipeline.
age_map <- data.frame(
  age_index = 1:20,
  AgeCat = c(
    "1-11", "1-11", "1-11", "1-11",
    "12-17",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "65+", "65+", "65+", "65+", "65+"
  )
)
# Raw cumulative infection counts (not normalised per-capita) -- shape over
# time is what matters here, not absolute scale.
cat("Loading reference postsim (~2.9GB, one-time)...\n")
load("01_Data/postsim_vc_ixchiq_model_finite.RData")
cat("Loaded.\n")

regions <- names(postsim_vc_ixchiq_model)

extract_cum_attack_rate <- function(scenario_int, age_label) {
  purrr::map_dfr(regions, function(reg) {
    reg_list <- postsim_vc_ixchiq_model[[reg]]
    purrr::imap_dfr(reg_list, function(ve_list, ve_name) {
      sc <- ve_list[["cov50"]]$scenario_result[[scenario_int]]
      raw_inf <- sc$sim_result$age_array_raw_inf   # [age_index, week, draw]

      age_idx <- age_map$age_index[age_map$AgeCat == age_label]
      sub_arr <- raw_inf[age_idx, , , drop = FALSE]
      by_week_draw <- apply(sub_arr, c(2, 3), sum, na.rm = TRUE)  # [week x draw]
      cum_by_week_draw <- apply(by_week_draw, 2, cumsum)

      med_by_week <- apply(cum_by_week_draw, 1, median, na.rm = TRUE)
      tibble(Region = reg, VE = ve_name, week = seq_len(nrow(cum_by_week_draw)),
             cum_infections_med = med_by_week)
    })
  }) %>% mutate(AgeCat = age_label)
}

plot_df <- bind_rows(
  extract_cum_attack_rate(3, "18-64"),
  extract_cum_attack_rate(4, "65+")
) %>%
  mutate(
    setting = factor(unname(setting_key[Region]), levels = setting_levels),
    Mechanism = factor(VE, levels = c("VE0", "VE98.9"),
                        labels = c("Disease blocking only", "Disease and infection blocking")),
    AgeCat = factor(AgeCat, levels = c("18-64", "65+"), labels = c("18-64 years", "65+ years"))
  ) %>%
  filter(!is.na(setting)) %>%
  group_by(setting, AgeCat, Mechanism, week) %>%
  summarise(cum_infections_med = median(cum_infections_med, na.rm = TRUE), .groups = "drop")

setting_colors <- c("Low" = "#1b9e77", "Moderate" = "#d95f02", "High" = "#d7191c")

theme_ar <- function(base_size = 14) {
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

p <- ggplot(plot_df, aes(x = week, y = cum_infections_med, colour = setting, group = setting)) +
  geom_line(linewidth = 1) +
  facet_grid(AgeCat ~ Mechanism) +
  scale_colour_manual(values = setting_colors, name = "Transmission setting") +
  scale_x_continuous(name = "Calendar week", breaks = seq(0, 52, by = 10)) +
  scale_y_continuous(name = "Cumulative new infections (population, raw count)", labels = label_number(scale_cut = cut_short_scale())) +
  labs(
    title = "Population-wide cumulative infection dynamics (R compartment proxy)",
    subtitle = "NOT split by vaccination status -- whole-population attack rate buildup over the year, reference scenario"
  ) +
  theme_ar()

ggsave("06_Results/cumulative_attack_rate_dynamics.png", p, width = 10, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/cumulative_attack_rate_dynamics.png")
