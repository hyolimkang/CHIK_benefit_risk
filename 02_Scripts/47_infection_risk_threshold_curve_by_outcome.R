# =============================================================================
# 47_infection_risk_threshold_curve_by_outcome.R
#
# SAE / Death companions to the DALY threshold curve (Figure 1B --
# 39_figure2_infection_risk_threshold_curve.R, stacked as in
# 42_figure2_grid_layout.R), plus one combined DALY | SAE | Death figure.
#
# Same data, filter (duration <= 90 d, entry week <= 40), binning and
# isotonic smoothing as script 39 -- only the Pr(BRR > 1) column changes.
#
# Death, 18-64: vaccine-attributable death risk is fixed at 0 for this age
# group (01_setup.R: p_death_vacc_u65 <- 0), so BRR is undefined and that
# panel is left empty.
#
# Output: 07_Final_Results/Supplementary_Figures/
#           infection_risk_threshold_curve_SAE.png
#           infection_risk_threshold_curve_Death.png
#         07_Final_Results/Main_Figures/infection_risk_threshold_curve_DALY_SAE_Death.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

load("01_Data/psa_grid_bra_travel_finite_setting.RData")

OUTCOME_COLS <- c(DALY = "brr_daly_prob_gt1", SAE = "brr_sae_prob_gt1", Death = "brr_death_prob_gt1")
SETTING_COLS <- c(Low = "#2A9D8F", Moderate = "#F4A261", High = "#C44E52")

d_long <- psa_grid_travel_setting %>%
  mutate(
    entry_week = (entry_day - 1) %/% 7 + 1,
    age_group = factor(age_group, levels = c("18-64", "65+"), labels = c("18-64", "\u226565")),
    setting = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  filter(duration <= 90, entry_week <= 40, is.finite(AR_median)) %>%
  select(age_group, setting, AR_median, all_of(unname(OUTCOME_COLS))) %>%
  pivot_longer(all_of(unname(OUTCOME_COLS)), names_to = "col", values_to = "prob") %>%
  mutate(outcome = factor(names(OUTCOME_COLS)[match(col, OUTCOME_COLS)], levels = names(OUTCOME_COLS))) %>%
  filter(is.finite(prob))

# ---- Binning + isotonic smoothing (verbatim logic from script 39) ----------
bin_curve <- function(data, groups, n_bins) {
  data %>%
    group_by(across(all_of(groups))) %>%
    mutate(risk_bin = ntile(AR_median, n_bins)) %>%
    group_by(across(all_of(c(groups, "risk_bin")))) %>%
    summarise(infection_risk = mean(AR_median), probability = mean(prob), .groups = "drop") %>%
    filter(is.finite(infection_risk), is.finite(probability)) %>%
    arrange(across(all_of(groups)), infection_risk)
}

make_smooth_curve <- function(data, groups) {
  data %>%
    group_by(across(all_of(groups))) %>%
    group_modify(~ {
      o <- order(.x$infection_risk)
      fit <- isoreg(.x$infection_risk[o], .x$probability[o])
      fitted <- tibble(infection_risk = .x$infection_risk[o], probability = fit$yf) %>%
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

smooth_curve <- make_smooth_curve(bin_curve(d_long, c("outcome", "age_group"), 50),
                                  c("outcome", "age_group"))
smooth_curve_setting <- make_smooth_curve(bin_curve(d_long, c("outcome", "age_group", "setting"), 35),
                                          c("outcome", "age_group", "setting"))
endpoints <- smooth_curve_setting %>%
  group_by(outcome, age_group, setting) %>%
  slice_max(infection_risk, n = 1, with_ties = FALSE) %>%
  ungroup()


theme_threshold <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
      panel.border = element_rect(colour = "#c3c2b7", fill = NA, linewidth = 0.3),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 3, colour = "#0b0b0b"),
      axis.title = element_text(face = "bold", colour = "#0b0b0b"),
      axis.text = element_text(colour = "#0b0b0b", size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 2),
      panel.spacing = unit(16, "pt"),
      plot.margin = margin(10, 12, 8, 8)
    )
}

build_plot <- function(outcomes, facet, title, y_name, base_size = 12) {
  sel <- function(df) filter(df, outcome %in% outcomes)
  ggplot() +
    geom_line(data = sel(smooth_curve_setting),
              aes(infection_risk, probability, colour = setting, group = setting),
              linewidth = 1.2, alpha = 0.9) +
    geom_point(data = sel(endpoints),
               aes(infection_risk, probability, fill = setting),
               shape = 21, colour = "white", stroke = 0.5, size = 3.2, alpha = 0.9) +
    geom_line(data = sel(smooth_curve), aes(infection_risk, probability),
              colour = "#111111", linewidth = 1.6) +
    facet +
    scale_colour_manual(values = SETTING_COLS, name = "Transmission setting", drop = FALSE) +
    scale_fill_manual(values = SETTING_COLS, guide = "none") +
    scale_x_continuous(
      name = "Cumulative probability of CHIKV infection during travel (%)",
      # 0.5% steps up to 2.0% only: a 2.5% tick would collide with the next
      # panel's 0% tick once fonts are enlarged.
      breaks = seq(0, 0.02, 0.005),
      labels = function(x) paste0(format(100 * x, drop0trailing = TRUE, trim = TRUE), "%"),
      limits = c(0, NA), expand = expansion(mult = c(0.03, 0.07))
    ) +
    scale_y_continuous(
      name = y_name, labels = percent_format(accuracy = 1),
      limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = c(0, 0)
    ) +
    labs(title = title) +
    theme_threshold(base_size)
}

out_dir <- "07_Final_Results/Supplementary_Figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Individual SAE / Death figures (stacked, like Figure 1B) --------------
for (oc in c("SAE", "Death")) {
  p <- build_plot(
    outcomes = oc,
    facet = facet_wrap(~age_group, ncol = 1, drop = FALSE),
    title = sprintf("Probability that %s benefit exceeds risk\nacross cumulative infection risk", oc),
    y_name = sprintf("Probability that %s BRR > 1", oc)
  )
  out_file <- file.path(out_dir, sprintf("infection_risk_threshold_curve_%s.png", oc))
  ggsave(out_file, p, width = 6.5, height = 9, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}

# ---- Combined DALY | SAE | Death ------------------------------------------
p_all <- build_plot(
  outcomes = names(OUTCOME_COLS),
  facet = facet_grid(age_group ~ outcome, drop = FALSE),
  title = NULL,  # main figure: title lives in the caption
  y_name = "Probability that BRR > 1",
  base_size = 17
)
# Main figure (landscape page in Word) -> Main_Figures, wider canvas.
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
out_file <- "07_Final_Results/Main_Figures/infection_risk_threshold_curve_DALY_SAE_Death.png"
ggsave(out_file, p_all, width = 16, height = 8.5, dpi = 300, bg = "#fcfcfb")
message("Saved: ", out_file)
