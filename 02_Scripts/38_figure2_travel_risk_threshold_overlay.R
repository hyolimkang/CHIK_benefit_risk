# =============================================================================
# 38_figure2_travel_risk_threshold_overlay.R
#
# Companion to Figure 2: keep the infection-risk surface as the fill and add
# the travel-condition boundary where median DALY benefit equals vaccine risk
# (median DALY BRR = 1). The x/y coordinates remain travel duration x entry
# week; infection risk is not moved to an axis.
#
# Output:
#   06_Results/figure2_travel_risk_surfaces_threshold_overlay.png
#   07_Final_Results/Main_Figures/figure2_travel_risk_surfaces_threshold_overlay.png
# =============================================================================

project_root <- getwd()
if (!file.exists(file.path(project_root, "01_Data/psa_grid_bra_travel_finite_setting.RData"))) {
  stop("Run this script from the project root: ", project_root)
}

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

load("01_Data/psa_grid_bra_travel_finite_setting.RData")

psa_grid_travel_setting <- psa_grid_travel_setting %>%
  mutate(
    entry_week = (entry_day - 1) %/% 7 + 1,
    age_group = factor(age_group, levels = c("18-64", "65+")),
    setting = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  filter(duration <= 90, entry_week <= 40)

ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_gridline <- "#e1e0d9"
ink_baseline <- "#c3c2b7"
threshold_young <- "#111111"
threshold_older <- "#ffffff"

theme_hm <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = ink_baseline, fill = NA, linewidth = 0.3),
      panel.spacing = unit(8, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 2, colour = ink_primary),
      axis.title = element_text(face = "bold", size = base_size + 2, colour = ink_primary),
      axis.text = element_text(size = base_size + 1, colour = ink_primary),
      axis.ticks = element_line(colour = ink_baseline, linewidth = 0.4),
      legend.title = element_text(face = "bold", size = base_size + 1, colour = ink_primary),
      legend.text = element_text(size = base_size, colour = ink_primary),
      legend.key.width = unit(14, "pt"),
      legend.key.height = unit(22, "pt"),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = base_size + 2, colour = ink_primary),
      plot.subtitle = element_text(size = base_size + 1, face = "bold", colour = ink_primary,
                                   margin = margin(b = 5)),
      plot.caption = element_text(size = base_size - 2, colour = ink_secondary,
                                  hjust = 0, lineheight = 1.2, margin = margin(t = 8)),
      plot.caption.position = "plot",
      plot.margin = margin(10, 12, 8, 8)
    )
}

interp_duration <- function(df, yvar, xout) {
  ok <- is.finite(df$duration) & is.finite(df[[yvar]])
  if (sum(ok) < 2) return(tibble::tibble(duration = xout, value = NA_real_))
  fit <- approx(x = df$duration[ok], y = df[[yvar]][ok], xout = xout, rule = 1)
  tibble::tibble(duration = fit$x, value = fit$y)
}

duration_uniform <- seq(min(psa_grid_travel_setting$duration),
                        max(psa_grid_travel_setting$duration), by = 1)

# Infection-risk surface: age-independent, so use one age group only.
risk_A <- psa_grid_travel_setting %>%
  filter(age_group == "18-64") %>%
  group_by(setting, entry_week, duration) %>%
  summarise(AR_median = median(AR_median, na.rm = TRUE), .groups = "drop") %>%
  group_by(setting, entry_week) %>%
  group_modify(~ interp_duration(.x, "AR_median", duration_uniform)) %>%
  ungroup() %>%
  rename(AR_median = value) %>%
  mutate(AR_median = pmax(AR_median, 0))

# Threshold surface: retain age groups because DALY vaccine risk differs by age.
threshold <- psa_grid_travel_setting %>%
  group_by(setting, age_group, entry_week, duration) %>%
  summarise(brr_daly_median = median(brr_daly_median, na.rm = TRUE), .groups = "drop") %>%
  group_by(setting, age_group, entry_week) %>%
  group_modify(~ interp_duration(.x, "brr_daly_median", duration_uniform)) %>%
  ungroup() %>%
  rename(brr_daly_median = value)

p_A <- ggplot(risk_A, aes(duration, entry_week, fill = AR_median)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(
    data = threshold,
    aes(x = duration, y = entry_week, z = brr_daly_median, linetype = age_group),
    breaks = 1,
    colour = ink_primary,
    linewidth = 0.55,
    inherit.aes = FALSE
  ) +
  facet_wrap(~setting, nrow = 1) +
  scale_fill_viridis_c(
    option = "plasma",
    labels = scales::percent_format(accuracy = 0.5),
    name = "P(infection)",
    na.value = ink_gridline
  ) +
  scale_linetype_manual(
    name = "Median DALY BRR = 1",
    values = c("18-64" = "solid", "65+" = "dashed")
  ) +
  scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
  scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 40, 10)) +
  labs(
    title = "A. Infection-risk surface with vaccine benefit-risk threshold",
    subtitle = "Contour: median DALY vaccine benefit equals vaccine risk (BRR = 1)",
    caption = "The heatmap remains P(infection); the contour marks the travel conditions where vaccine benefit changes from below to above vaccine risk."
  ) +
  theme_hm()

# Keep the original Figure 2 Panel B unchanged, for direct visual comparison.
risk_B <- psa_grid_travel_setting %>%
  mutate(prob_gt1 = brr_daly_prob_gt1) %>%
  group_by(setting, age_group, entry_week) %>%
  group_modify(~ interp_duration(.x, "prob_gt1", duration_uniform)) %>%
  ungroup() %>%
  rename(prob_gt1 = value) %>%
  mutate(prob_gt1 = pmin(pmax(prob_gt1, 0), 1))

p_B <- ggplot(risk_B, aes(duration, entry_week, fill = prob_gt1)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = prob_gt1), breaks = 0.5, colour = ink_primary,
               linewidth = 0.35, linetype = "42") +
  facet_grid(age_group ~ setting) +
  scale_fill_gradient2(
    name = "Pr(BRR>1)", low = "#e34948", mid = "#f0efec", high = "#2a78d6",
    midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1), na.value = ink_gridline
  ) +
  scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
  scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 40, 10)) +
  labs(title = "B. DALY benefit-risk surface") +
  theme_hm()

p_fig <- (p_A / p_B) + plot_layout(heights = c(1, 2))

out_files <- c(
  "06_Results/figure2_travel_risk_surfaces_threshold_overlay.png",
  "07_Final_Results/Main_Figures/figure2_travel_risk_surfaces_threshold_overlay.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p_fig, width = 11, height = 11, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}