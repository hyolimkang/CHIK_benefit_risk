# =============================================================================
# 43_figure2B_latest_favourable_start_week.R
#
# Figure 2B decision-boundary summary:
#   x = vaccine coverage (%)
#   y = latest campaign start week with median BRR >= 1
#   columns = Low / Moderate / High setting
#   rows = 18-64 / >=65 years
#   colour = protection mechanism
#   ribbon = range between base and serostatus-adjusted risk assumptions
#
# Boundary values are calculated only at the observed coverage values and are
# connected for display. No week or coverage values are extrapolated.
# =============================================================================

if (!file.exists("06_Results/brr_weeksweep_summary_finite.xlsx") ||
    !file.exists("06_Results/brr_weeksweep_summary_finite_ve0.xlsx")) {
  stop("ORI summary Excel files are missing from 06_Results/")
}

suppressMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

summary_data <- bind_rows(
  read_excel("06_Results/brr_weeksweep_summary_finite.xlsx"),
  read_excel("06_Results/brr_weeksweep_summary_finite_ve0.xlsx")
) %>%
  filter(outcome == "DALY", RR_seropos == 0,
         AgeCat %in% c("18-64", "65+")) %>%
  mutate(
    age_group = factor(AgeCat, levels = c("18-64", "65+"),
                       labels = c("18-64", "≥65")),
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    mechanism = factor(VE_label,
      levels = c("Disease blocking only", "Disease and infection blocking")),
    coverage = as.numeric(sub("cov", "", Coverage)),
    week = as.numeric(week)
  )

latest_week <- function(week, brr) {
  eligible <- week[is.finite(brr) & brr >= 1]
  if (!length(eligible)) return(NA_real_)
  max(eligible)
}

boundary <- summary_data %>%
  group_by(age_group, setting, mechanism, coverage) %>%
  summarise(
    base_week = latest_week(week, brr_base_med),
    adjusted_week = latest_week(week, brr_adj_med),
    .groups = "drop"
  ) %>%
  mutate(
    week_low = pmin(base_week, adjusted_week, na.rm = TRUE),
    week_high = pmax(base_week, adjusted_week, na.rm = TRUE),
    week_mid = rowMeans(cbind(base_week, adjusted_week), na.rm = TRUE),
    week_low = ifelse(is.finite(week_low), week_low, NA_real_),
    week_high = ifelse(is.finite(week_high), week_high, NA_real_),
    week_mid = ifelse(is.finite(week_mid), week_mid, NA_real_)
  )

# Remove only rows with no boundary under either risk assumption. A missing
# edge remains missing rather than being treated as a zero-week boundary.
line_data <- boundary %>%
  filter(is.finite(week_mid)) %>%
  arrange(age_group, setting, mechanism, coverage)

reference <- data.frame(coverage = 50, week = 2)

mechanism_cols <- c(
  "Disease blocking only" = "#2a9d8f",
  "Disease and infection blocking" = "#c44e52"
)

p <- ggplot(line_data, aes(coverage, week_mid, colour = mechanism,
                           group = mechanism)) +
  geom_ribbon(
    aes(ymin = week_low, ymax = week_high, fill = mechanism),
    alpha = 0.16, colour = NA, inherit.aes = TRUE
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  geom_point(
    data = reference, aes(coverage, week), inherit.aes = FALSE,
    shape = 23, size = 2.1, stroke = 0.35,
    fill = "#f5c518", colour = "#6e5700"
  ) +
  facet_grid(age_group ~ setting) +
  scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
  scale_fill_manual(values = mechanism_cols, guide = "none") +
  scale_x_continuous(
    name = "Vaccine coverage (%)", breaks = c(10, 30, 50, 70, 90),
    limits = c(10, 90), expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Latest campaign start week with median BRR ≥ 1",
    breaks = c(1, 8, 16, 24, 32, 40, 52), limits = c(1, 52),
    expand = c(0, 0)
  ) +
  labs(
    title = "B. Latest favourable vaccination start week by coverage",
    subtitle = "Ribbon: range between base-risk and serostatus-adjusted-risk boundaries; diamond: reference scenario (week 2, 50% coverage)"
  ) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#c8c8c8", fill = NA, linewidth = 0.35),
    panel.spacing.x = unit(18, "pt"),
    panel.spacing.y = unit(10, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "#1a1a1a"),
    axis.title = element_text(face = "bold", colour = "#1a1a1a"),
    axis.text = element_text(colour = "#1a1a1a"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(colour = "#5a5a5a", size = 10),
    plot.margin = margin(10, 12, 8, 8)
  )

out_files <- c(
  "06_Results/figure2B_latest_favourable_start_week.png",
  "07_Final_Results/Main_Figures/figure2B_latest_favourable_start_week.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p, width = 11, height = 8.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}
