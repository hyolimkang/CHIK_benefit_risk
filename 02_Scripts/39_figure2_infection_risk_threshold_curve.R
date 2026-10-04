# =============================================================================
# 39_figure2_infection_risk_threshold_curve.R
#
# Two-panel companion figure:
#   x = cumulative probability of CHIKV infection during travel
#   y = probability that vaccination benefit exceeds risk, Pr(BRR > 1)
#
# Panels are age groups only. The black curve pools the travel grid across
# settings. Reference lines and threshold annotations are intentionally
# omitted; the figure focuses on the pooled and setting-specific relationships.
# =============================================================================

project_root <- getwd()
if (!file.exists(file.path(project_root, "01_Data/psa_grid_bra_travel_finite_setting.RData"))) {
  stop("Run this script from the project root: ", project_root)
}

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

load("01_Data/psa_grid_bra_travel_finite_setting.RData")

d <- psa_grid_travel_setting %>%
  mutate(
    entry_week = (entry_day - 1) %/% 7 + 1,
    age_group = factor(age_group, levels = c("18-64", "65+"),
               labels = c("18-64", "\u226565")),
    setting = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  filter(duration <= 90, entry_week <= 40,
         is.finite(AR_median), is.finite(brr_daly_prob_gt1))

# Bin infection risk within the observed data. The pooled bin summaries define
# the black relationship.
curve <- d %>%
  group_by(age_group) %>%
  mutate(risk_bin = ntile(AR_median, 50)) %>%
  group_by(age_group, risk_bin) %>%
  summarise(
    infection_risk = mean(AR_median),
    probability = mean(brr_daly_prob_gt1),
    .groups = "drop"
  ) %>%
  filter(is.finite(infection_risk), is.finite(probability)) %>%
  arrange(age_group, infection_risk)

curve_setting <- d %>%
  group_by(age_group, setting) %>%
  mutate(risk_bin = ntile(AR_median, 35)) %>%
  group_by(age_group, setting, risk_bin) %>%
  summarise(
    infection_risk = mean(AR_median),
    probability = mean(brr_daly_prob_gt1),
    .groups = "drop"
  ) %>%
  filter(is.finite(infection_risk), is.finite(probability)) %>%
  arrange(age_group, setting, infection_risk)

make_smooth_curve <- function(data, groups) {
  data %>%
  group_by(across(all_of(groups))) %>%
  group_modify(~ {
    o <- order(.x$infection_risk)
    fit <- isoreg(.x$infection_risk[o], .x$probability[o])
    fitted <- tibble(
      infection_risk = .x$infection_risk[o],
      probability = fit$yf
    ) %>%
      group_by(infection_risk) %>%
      summarise(probability = mean(probability), .groups = "drop")
    xout <- seq(min(fitted$infection_risk), max(fitted$infection_risk), length.out = 300)
    tibble(
      infection_risk = xout,
      probability = pmin(pmax(approx(fitted$infection_risk, fitted$probability,
                                     xout, rule = 2)$y, 0), 1)
    )
  }) %>%
  ungroup()
}

smooth_curve <- make_smooth_curve(curve, "age_group")
smooth_curve_setting <- make_smooth_curve(curve_setting, c("age_group", "setting"))

endpoints <- smooth_curve_setting %>%
  group_by(age_group, setting) %>%
  slice_max(infection_risk, n = 1, with_ties = FALSE) %>%
  ungroup()

theme_threshold <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
      panel.border = element_rect(colour = "#c3c2b7", fill = NA, linewidth = 0.3),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 2, colour = "#0b0b0b"),
      axis.title = element_text(face = "bold", colour = "#0b0b0b"),
      axis.text = element_text(colour = "#0b0b0b"),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(face = "bold", colour = "#52514e"),
      plot.caption = element_text(colour = "#52514e", hjust = 0, size = base_size - 2),
      plot.caption.position = "plot",
      panel.spacing = unit(12, "pt"),
      plot.margin = margin(10, 12, 8, 8)
    )
}

p <- ggplot(d, aes(AR_median, brr_daly_prob_gt1)) +
  geom_line(
    data = smooth_curve_setting,
    aes(infection_risk, probability, colour = setting, group = setting),
    linewidth = 1.0, alpha = 0.9, inherit.aes = FALSE
  ) +
  geom_point(
    data = endpoints,
    aes(infection_risk, probability, fill = setting),
    shape = 21, colour = "white", stroke = 0.5, size = 3.2, alpha = 0.9,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = smooth_curve,
    aes(infection_risk, probability),
    colour = "#111111", linewidth = 1.3, inherit.aes = FALSE
  ) +
  facet_wrap(~age_group, nrow = 1) +
  scale_colour_manual(
    values = c(Low = "#2A9D8F", Moderate = "#F4A261", High = "#C44E52"),
    name = "Transmission setting"
  ) +
  scale_fill_manual(
    values = c(Low = "#2A9D8F", Moderate = "#F4A261", High = "#C44E52"),
    guide = "none"
  ) +
  scale_x_continuous(
    name = "Cumulative probability of CHIKV infection during travel (%)",
    labels = percent_format(accuracy = 0.1),
    limits = c(0, NA), expand = expansion(mult = c(0.03, 0.07))
  ) +
  scale_y_continuous(
    name = "Probability that DALY BRR > 1",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)
  ) +
  labs(
    title = "B. Probability that DALY benefit exceeds risk across cumulative infection risk"
  ) +
  theme_threshold()

out_files <- c(
  "06_Results/figure2_infection_risk_threshold_curve.png",
  "07_Final_Results/Main_Figures/figure2_infection_risk_threshold_curve.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p, width = 11, height = 6.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}
