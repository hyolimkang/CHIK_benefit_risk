# =============================================================================
# 37_figure2_travel_risk_surface_sae_death.R
#
# Supplementary companion to Figure 2 (17_figure2_travel_risk_surface.R,
# DALY only): a single 3-panel figure --
#   A. Infection-risk surface (age-independent, identical to Figure 2 Panel A)
#   B. SAE benefit-risk surface (Pr(BRR>1), 18-64 x 65+ x setting)
#   C. Death benefit-risk surface (Pr(BRR>1), 65+ only x setting)
#
# Death differs structurally from DALY/SAE: vaccine-attributable death risk
# for 18-64 is fixed at exactly 0 (see 01_setup.R's p_death_vacc_u65 <- 0),
# so brr_death_prob_gt1 is NA for every 18-64 cell (verified: 0/4368 finite
# vs 3345/4368 for 65+) -- Panel C drops the 18-64 row entirely
# (facet_wrap(~setting) instead of facet_grid(age_group~setting)), not just
# leaving it blank/grey (chat record 2026-09-02: previously a separate
# Death-only figure; folded in here as Panel C under the SAE figure instead).
#
# Data source: 01_Data/psa_grid_bra_travel_finite_setting.RData
#   (psa_grid_travel_setting; same object 17_figure2_travel_risk_surface.R
#   uses -- brr_sae_prob_gt1 / brr_death_prob_gt1 columns already present).
#
# Output: 07_Final_Results/Supplementary_Figures/figure2_travel_risk_surfaces_SAE.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

load("01_Data/psa_grid_bra_travel_finite_setting.RData")  # -> psa_grid_travel_setting

psa_grid_travel_setting <- psa_grid_travel_setting %>%
  mutate(entry_week = (entry_day - 1) %/% 7 + 1,
         age_group  = factor(age_group, levels = c("18-64", "65+")),
         setting    = factor(setting, levels = c("Low", "Moderate", "High")))

# ---- Shared styling (matches 17_figure2_travel_risk_surface.R exactly) -----
ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_gridline  <- "#e1e0d9"
ink_baseline  <- "#c3c2b7"
pole_blue     <- "#2a78d6"
pole_red      <- "#e34948"
mid_grey      <- "#f0efec"

theme_hm <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = ink_baseline, fill = NA, linewidth = 0.3),
      panel.spacing    = unit(8, "pt"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size + 2, colour = ink_primary),
      axis.title       = element_text(face = "bold", size = base_size + 2, colour = ink_primary),
      axis.text        = element_text(size = base_size + 1, colour = ink_primary),
      axis.ticks       = element_line(colour = ink_baseline, linewidth = 0.4),
      legend.title      = element_text(face = "bold", size = base_size + 1, colour = ink_primary),
      legend.text       = element_text(size = base_size, colour = ink_primary),
      legend.key.width  = unit(14, "pt"),
      legend.key.height = unit(22, "pt"),
      legend.position   = "right",
      plot.title    = element_text(face = "bold", size = base_size + 2, colour = ink_primary),
      plot.subtitle = element_text(size = base_size + 1, face = "bold", colour = ink_primary, margin = margin(b = 5)),
      plot.caption  = element_text(size = base_size - 2, colour = ink_secondary,
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

psa_grid_travel_setting <- psa_grid_travel_setting %>%
  filter(duration <= 90, entry_week <= 40)

duration_uniform <- seq(min(psa_grid_travel_setting$duration), max(psa_grid_travel_setting$duration), by = 1)

# ---- Panel A -- infection-risk surface (age-independent) -------------------
risk_A <- psa_grid_travel_setting %>%
  filter(age_group == "18-64") %>%
  group_by(setting, entry_week) %>%
  group_modify(~ interp_duration(.x, "AR_median", duration_uniform)) %>%
  ungroup() %>%
  rename(AR_median = value) %>%
  mutate(AR_median = pmax(AR_median, 0))

p_A <- ggplot(risk_A, aes(duration, entry_week, fill = AR_median)) +
  geom_raster(interpolate = TRUE) +
  facet_wrap(~setting, nrow = 1) +
  scale_fill_viridis_c(option = "plasma", labels = scales::percent_format(accuracy = 0.5),
                        name = "P(infection)", na.value = ink_gridline) +
  scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
  scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 40, 10)) +
  labs(title = "A. Infection-risk surface") +
  theme_hm()

# ---- Panel builder -- Pr(<outcome> BRR > 1) surface -------------------------
build_prob_panel <- function(prob_col, age_groups_keep, panel_letter, outcome_label) {
  d <- psa_grid_travel_setting %>%
    filter(age_group %in% age_groups_keep) %>%
    mutate(prob_gt1 = .data[[prob_col]]) %>%
    group_by(setting, age_group, entry_week) %>%
    group_modify(~ interp_duration(.x, "prob_gt1", duration_uniform)) %>%
    ungroup() %>%
    rename(prob_gt1 = value) %>%
    mutate(prob_gt1 = pmin(pmax(prob_gt1, 0), 1),
           age_group = factor(age_group, levels = age_groups_keep))

  p <- ggplot(d, aes(duration, entry_week, fill = prob_gt1)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = prob_gt1), breaks = 0.5, colour = ink_primary, linewidth = 0.35, linetype = "42") +
    scale_fill_gradient2(
      name = "Pr(BRR>1)",
      low = pole_red, mid = mid_grey, high = pole_blue,
      midpoint = 0.5, limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = scales::percent_format(accuracy = 1),
      na.value = ink_gridline
    ) +
    scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
    scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 40, 10)) +
    labs(title = paste0(panel_letter, ". ", outcome_label, " benefit-risk surface")) +
    theme_hm()

  if (length(age_groups_keep) > 1) {
    p <- p + facet_grid(age_group ~ setting)
  } else {
    p <- p + facet_wrap(~setting, nrow = 1)
  }
  p
}

# Panel B: SAE, both age groups (matches DALY's 2-row layout).
p_B <- build_prob_panel("brr_sae_prob_gt1", c("18-64", "65+"), "B", "SAE")

# Panel C: Death, 65+ only -- vaccine-attributable death risk is fixed at 0
# for 18-64, so brr_death_prob_gt1 is undefined (NA) for every 18-64 cell;
# showing an all-grey row would misrepresent "no data" as "no signal".
p_C <- build_prob_panel("brr_death_prob_gt1", "65+", "C", "Death")

dir.create("07_Final_Results/Supplementary_Figures", showWarnings = FALSE, recursive = TRUE)

p_fig <- (p_A / p_B / p_C) +
  plot_layout(heights = c(1, 2, 1)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, colour = ink_primary),
      plot.caption = element_text(size = 8, colour = ink_secondary, hjust = 0)
    )
  )

out_file <- "07_Final_Results/Supplementary_Figures/figure2_travel_risk_surfaces_SAE.png"
ggsave(out_file, p_fig, width = 11, height = 15.5, dpi = 300, bg = "#fcfcfb")
message("Saved: ", out_file)

# Remove the now-superseded standalone Death figure (folded into Panel C above).
old_death_file <- "07_Final_Results/Supplementary_Figures/figure2_travel_risk_surfaces_Death.png"
if (file.exists(old_death_file)) {
  file.remove(old_death_file)
  message("Removed (superseded by Panel C): ", old_death_file)
}
