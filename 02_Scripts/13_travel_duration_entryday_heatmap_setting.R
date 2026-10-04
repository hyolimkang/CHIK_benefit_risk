# =============================================================================
# 13_travel_duration_entryday_heatmap_setting.R
#
# Summarised version of 12_travel_duration_entryday_heatmap.R: instead of one
# panel per state (11 panels), pools states WITHIN each Low/Moderate/High
# transmission setting (setting_key.RData, case-rate-based classification)
# together with draws before summarising -- so 3 settings x 2 age groups = 6
# panels total per outcome, not 22.
#
# Source: 01_Data/psa_grid_bra_travel_finite_setting.RData (psa_grid_travel_setting),
# built by the setting-level pooling block added to 06_brazil_travel_final_finite.R
# (pools states+draws jointly within each setting -- NOT a median-of-state-medians
# approximation, since the raw per-state x per-draw grid arrays were still in
# memory at that point in 06's script).
#
# x = travel duration (7-180 days), y = entry week of year (1-52), fill/contour
# = median BRR, dashed contour at BRR = 1.
#
# Output: 06_Results/travel_duration_entryday_setting_{sae,death,daly}.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")

if (!file.exists("01_Data/psa_grid_bra_travel_finite_setting.RData")) {
  stop("01_Data/psa_grid_bra_travel_finite_setting.RData not found -- run 06_brazil_travel_final_finite.R first.")
}

library(dplyr)
library(ggplot2)
library(scales)

load("01_Data/psa_grid_bra_travel_finite_setting.RData")  # -> psa_grid_travel_setting

psa_grid_travel_setting <- psa_grid_travel_setting %>%
  mutate(entry_week = (entry_day - 1) %/% 7 + 1,
         age_group  = factor(age_group, levels = c("18-64", "65+")),
         setting    = factor(setting, levels = c("Low", "Moderate", "High")))

# duration_grid is unevenly spaced (7-day steps + 30/90/180 forced in), which
# breaks geom_raster's even-tile assumption -- interpolate onto a uniform
# 1-day duration axis (on the log10(BRR) scale) before rasterising, same fix
# as 12_travel_duration_entryday_heatmap.R.
interp_duration <- function(df, xout) {
  ok <- is.finite(df$duration) & is.finite(df$log_brr)
  if (sum(ok) < 2) return(tibble::tibble(duration = xout, log_brr = NA_real_))
  fit <- approx(x = df$duration[ok], y = df$log_brr[ok], xout = xout, rule = 2)
  tibble::tibble(duration = fit$x, log_brr = fit$y)
}

# ---- Nature Medicine-style chrome -------------------------------------------
# Diverging blue/red poles + neutral grey midpoint (BRR = 1), restrained ink
# tones, hairline panel borders, no chart title -- journal figures carry their
# explanation in an external caption, so the in-figure text stays a small,
# muted annotation rather than a headline.
ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
ink_gridline  <- "#e1e0d9"
ink_baseline  <- "#c3c2b7"
pole_blue     <- "#2a78d6"
pole_red      <- "#e34948"
mid_grey      <- "#f0efec"

theme_hm <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = ink_baseline, fill = NA, linewidth = 0.3),
      panel.spacing    = unit(6, "pt"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size + 0.5, colour = ink_primary),
      axis.title       = element_text(size = base_size, colour = ink_primary),
      axis.text        = element_text(size = base_size - 1, colour = ink_muted),
      axis.ticks       = element_line(colour = ink_baseline, linewidth = 0.3),
      legend.title      = element_text(size = base_size - 0.5, colour = ink_primary),
      legend.text       = element_text(size = base_size - 1.5, colour = ink_secondary),
      legend.key.width  = unit(9, "pt"),
      legend.key.height = unit(16, "pt"),
      legend.position   = "right",
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption  = element_text(size = base_size - 2.5, colour = ink_secondary,
                                    hjust = 0, lineheight = 1.15, margin = margin(t = 6)),
      plot.caption.position = "plot",
      plot.margin = margin(8, 10, 6, 6)
    )
}

brr_fill_scale <- function() {
  # Convention matched to the outbreak-response weeksweep heatmap
  # (finite_weeksweep_heatmap.R): BRR > 1 (benefit > risk) = blue,
  # BRR < 1 (risk > benefit) = red -- so low=red, high=blue here.
  scale_fill_gradient2(
    name = "BRR",
    trans = "log10",
    low = pole_red, mid = mid_grey, high = pole_blue,
    midpoint = 0,  # log10(1) = 0
    breaks = c(0.1, 1, 10),
    labels = label_number(accuracy = 0.1),
    na.value = ink_gridline
  )
}

build_setting_heatmap <- function(outcome_label, med_col) {
  duration_uniform <- seq(min(psa_grid_travel_setting$duration), max(psa_grid_travel_setting$duration), by = 1)

  d_raw <- psa_grid_travel_setting %>%
    mutate(brr_capped = pmin(pmax(.data[[med_col]], 0.01), 100),  # guard log-scale/contour against 0 or Inf
           log_brr = log10(brr_capped))

  d <- d_raw %>%
    group_by(setting, age_group, entry_week) %>%
    group_modify(~ interp_duration(.x, duration_uniform)) %>%
    ungroup() %>%
    mutate(brr_capped = 10^log_brr)

  p <- ggplot(d, aes(x = duration, y = entry_week, fill = brr_capped)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = brr_capped), breaks = 1, colour = ink_primary, linewidth = 0.35, linetype = "42") +
    facet_grid(age_group ~ setting) +
    brr_fill_scale() +
    scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
    scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 52, 13)) +
    labs(
      caption = paste0(
        outcome_label, " — median BRR, states pooled with draws within setting (Low/Moderate/High, case-rate classification).\n",
        "Dashed line = BRR = 1; blue side = benefit > risk. Finite-history baseline immunity, S0-denominator attack rate."
      )
    ) +
    theme_hm()

  out_file <- sprintf("06_Results/travel_duration_entryday_setting_%s.png",
                       tolower(gsub("[^A-Za-z]", "", outcome_label)))
  ggsave(out_file, p, width = 9.5, height = 5.9, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
  invisible(p)
}

outcomes <- list(SAE = "brr_sae_median", Death = "brr_death_median", DALY = "brr_daly_median")
for (oc_label in names(outcomes)) build_setting_heatmap(oc_label, outcomes[[oc_label]])
