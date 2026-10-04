# =============================================================================
# 31_total_seropositive_dynamics.R
#
# Total population seropositive share over calendar week = cumulative NATURAL
# infections (R, from age_array_raw_inf) + cumulative VACCINE-induced
# seroconversion (people dosed while still susceptible, from vacc_to_S_array
# -- the live-attenuated IXCHIQ vaccine is assumed to cause seroconversion
# regardless of the model's epidemiological VE_inf parameter, which tracks
# transmission-blocking protection, not serology -- so this is added even
# under the DB/VE0 mechanism, where the model's own V compartment would
# otherwise stay at 0). Chat record 2026-08-31.
#
# CAVEAT: this can slightly double-count anyone vaccinated while susceptible
# who later ALSO acquires a natural (breakthrough) infection -- the simulator
# doesn't track that overlap (see 29_seromix_within_campaign_dynamics.R's
# discussion), so under DB especially (no transmission-blocking at all) this
# is a mild upper bound, not an exact figure.
#
# N (population per age bin per region) isn't saved directly; recovered from
# the model's own internal consistency: at this reference scenario's 50%
# coverage target, total doses eventually given ≈ 0.5 * N (confirmed the
# rollout completes well within the 52-week window), so N = final cumulative
# dose count / 0.5.
#
# Output: 06_Results/total_seropositive_dynamics.png
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
COVERAGE_FRAC <- 0.5  # this reference scenario's coverage target

age_map <- data.frame(
  age_index = 1:20,
  AgeCat = c(
    "1-11", "1-11", "1-11", "1-11", "12-17",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "65+", "65+", "65+", "65+", "65+"
  )
)

cat("Loading reference postsim (~2.9GB, one-time)...\n")
load("01_Data/postsim_vc_ixchiq_model_finite.RData")
cat("Loaded.\n")

regions <- names(postsim_vc_ixchiq_model)

extract_total_seropos <- function(scenario_int, age_label) {
  age_idx <- age_map$age_index[age_map$AgeCat == age_label]

  purrr::map_dfr(regions, function(reg) {
    reg_list <- postsim_vc_ixchiq_model[[reg]]
    purrr::imap_dfr(reg_list, function(ve_list, ve_name) {
      sc <- ve_list[["cov50"]]$scenario_result[[scenario_int]]

      raw_inf   <- sc$sim_result$age_array_raw_inf[age_idx, , , drop = FALSE]
      raw_alloc <- sc$sim_result$raw_allocation_array[age_idx, , , drop = FALSE]
      vacc_to_S <- sc$sim_result$vacc_to_S_array[age_idx, , , drop = FALSE]

      inf_by_week_draw   <- apply(raw_inf,   c(2, 3), sum, na.rm = TRUE)
      alloc_by_week_draw <- apply(raw_alloc, c(2, 3), sum, na.rm = TRUE)
      vts_by_week_draw   <- apply(vacc_to_S, c(2, 3), sum, na.rm = TRUE)

      cum_inf   <- apply(inf_by_week_draw,   2, cumsum)  # natural infections (R)
      cum_alloc <- apply(alloc_by_week_draw, 2, cumsum)  # total doses given
      cum_vts   <- apply(vts_by_week_draw,   2, cumsum)  # vaccine-induced seroconversion

      n_weeks <- nrow(cum_inf)
      N_per_draw <- cum_alloc[n_weeks, ] / COVERAGE_FRAC  # recover population size
      N_med <- median(N_per_draw, na.rm = TRUE)

      total_seropos <- cum_inf + cum_vts  # [week x draw]
      seropos_pct <- 100 * total_seropos / matrix(N_per_draw, nrow = n_weeks, ncol = ncol(total_seropos), byrow = TRUE)
      med_by_week <- apply(seropos_pct, 1, median, na.rm = TRUE)

      tibble(Region = reg, VE = ve_name, week = seq_len(n_weeks), seropos_pct_med = med_by_week)
    })
  }) %>% mutate(AgeCat = age_label)
}

plot_df <- bind_rows(
  extract_total_seropos(3, "18-64"),
  extract_total_seropos(4, "65+")
) %>%
  mutate(
    setting = factor(unname(setting_key[Region]), levels = setting_levels),
    Mechanism = factor(VE, levels = c("VE0", "VE98.9"),
                        labels = c("Disease blocking only", "Disease and infection blocking")),
    AgeCat = factor(AgeCat, levels = c("18-64", "65+"), labels = c("18-64 years", "65+ years"))
  ) %>%
  filter(!is.na(setting)) %>%
  group_by(setting, AgeCat, Mechanism, week) %>%
  summarise(seropos_pct_med = median(seropos_pct_med, na.rm = TRUE), .groups = "drop")

setting_colors <- c("Low" = "#1b9e77", "Moderate" = "#d95f02", "High" = "#d7191c")

theme_sp <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", colour = NA),
      strip.text = element_text(face = "bold", size = rel(1.1)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle = element_text(colour = "grey30", size = rel(0.85))
    )
}

p <- ggplot(plot_df, aes(x = week, y = seropos_pct_med, colour = setting, group = setting)) +
  geom_line(linewidth = 1) +
  facet_grid(AgeCat ~ Mechanism) +
  scale_colour_manual(values = setting_colors, name = "Transmission setting") +
  scale_x_continuous(name = "Calendar week", breaks = seq(0, 52, by = 10)) +
  scale_y_continuous(name = "Total seropositive (%)", labels = label_number(suffix = "%")) +
  labs(
    title = "Total population seropositivity: natural infection + vaccine-induced",
    subtitle = "Reference scenario (week 2 start, 50% coverage). Mild upper bound under DB: breakthrough\ninfections among vaccinated-while-susceptible people aren't separately tracked (see chat record)."
  ) +
  theme_sp()

ggsave("06_Results/total_seropositive_dynamics.png", p, width = 10, height = 7.5, dpi = 300, bg = "white")
message("Saved: 06_Results/total_seropositive_dynamics.png")
