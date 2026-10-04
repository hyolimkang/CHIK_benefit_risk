# =============================================================================
# 42_figure2B_ori_decision_boundary.R
#
# Figure 2B: six-panel ORI decision-boundary plot.
#   Columns: Low / Moderate / High transmission setting
#   Rows: 18-64 / >=65 years
#   x: vaccination campaign start week
#   y: coverage
#   colour: vaccine protection mechanism
#   linetype: base versus serostatus-adjusted risk assumption
#
# Each line is the interpolated median BRR = 1 crossover contour. The yellow
# diamond marks the reference scenario (week 2, 50% coverage). No heatmap or
# extrapolation beyond the observed week x coverage grid is used.
# =============================================================================

if (!file.exists("06_Results/brr_weeksweep_summary_finite.xlsx") ||
    !file.exists("06_Results/brr_weeksweep_summary_finite_ve0.xlsx")) {
  stop("ORI summary Excel files are missing from 06_Results/")
}

suppressMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(mgcv)
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
    week = as.numeric(week),
    coverage_pct = as.numeric(sub("cov", "", Coverage))
  )

if (!nrow(summary_data)) stop("No ORI DALY rows survived filtering.")

# Fit within the observed 7-week x 4-coverage grid and predict only inside
# that rectangle. This supplies a smooth contour without making new support.
make_surface <- function(data, risk_assumption, value_col, grid_res = 120) {
  data <- data %>% filter(is.finite(.data[[value_col]]))
  week_range <- range(data$week)
  coverage_range <- range(data$coverage_pct)
  prediction_grid <- expand.grid(
    week = seq(week_range[1], week_range[2], length.out = grid_res),
    coverage_pct = seq(coverage_range[1], coverage_range[2], length.out = grid_res)
  )

  fit <- tryCatch(
    mgcv::gam(
      as.formula(paste(value_col, "~ te(week, coverage_pct, k = c(5, 4))")),
      data = data
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- mgcv::gam(
      as.formula(paste(value_col, "~ s(week, k = 5) + s(coverage_pct, k = 4)")),
      data = data
    )
  }

  prediction_grid$BRR <- predict(fit, newdata = prediction_grid)
  prediction_grid$BRR <- pmax(pmin(prediction_grid$BRR, 100), -100)
  prediction_grid$risk_assumption <- risk_assumption
  # Some mechanism/age/setting combinations never cross BRR = 1 within the
  # observed week x coverage rectangle. Return no contour data for those
  # combinations instead of asking ggplot2 to draw a nonexistent contour.
  if (max(prediction_grid$BRR, na.rm = TRUE) < 1 ||
      min(prediction_grid$BRR, na.rm = TRUE) > 1) {
    return(prediction_grid[0, ])
  }
  prediction_grid
}

boundary_data <- bind_rows(
  summary_data %>%
    group_by(age_group, setting, mechanism) %>%
    group_modify(~ make_surface(.x, "Base risk", "brr_base_med")) %>%
    ungroup(),
  summary_data %>%
    group_by(age_group, setting, mechanism) %>%
    group_modify(~ make_surface(.x, "Serostatus-adjusted risk", "brr_adj_med")) %>%
    ungroup()
) %>%
  mutate(risk_assumption = factor(
    risk_assumption,
    levels = c("Base risk", "Serostatus-adjusted risk")
  ))

reference <- data.frame(week = 2, coverage_pct = 50)

p <- ggplot(boundary_data, aes(week, coverage_pct)) +
  geom_contour(
    aes(z = BRR, colour = mechanism, linetype = risk_assumption),
    breaks = 1, linewidth = 0.75, na.rm = TRUE
  ) +
  geom_point(
    data = reference, aes(week, coverage_pct), inherit.aes = FALSE,
    shape = 23, size = 2.4, stroke = 0.35,
    fill = "#f5c518", colour = "#6e5700"
  ) +
  facet_grid(age_group ~ setting) +
  scale_colour_manual(
    values = c(
      "Disease blocking only" = "#2a9d8f",
      "Disease and infection blocking" = "#c44e52"
    ),
    name = "Protection mechanism"
  ) +
  scale_linetype_manual(
    values = c("Base risk" = "solid", "Serostatus-adjusted risk" = "dashed"),
    name = "Risk assumption"
  ) +
  scale_x_continuous(
    name = "Vaccination campaign start week",
    breaks = c(1, 13, 26, 39, 52), expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Coverage",
    breaks = c(10, 30, 60, 90), labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(1, 52), ylim = c(10, 90), expand = FALSE) +
  labs(
    title = "B. Vaccination benefit-risk decision boundaries",
    subtitle = "Lines show the median BRR = 1 crossover contour; diamond = reference scenario (week 2, 50% coverage)"
  ) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "#c8c8c8", fill = NA, linewidth = 0.35),
    panel.spacing = unit(10, "pt"),
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
  "06_Results/figure2B_ori_decision_boundaries.png",
  "07_Final_Results/Main_Figures/figure2B_ori_decision_boundaries.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p, width = 11, height = 8.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}
