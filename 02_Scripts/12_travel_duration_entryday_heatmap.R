# =============================================================================
# 12_travel_duration_entryday_heatmap.R
#
# x = travel duration (7-180 days), y = entry day (week of year, 1-52),
# fill/contour = median BRR, one panel per state, separate figures per
# outcome (SAE/Death/DALY) x age group (18-64/65+).
#
# Source: 01_Data/psa_grid_bra_travel_finite.RData (psa_grid_travel), built by
# the "Fine entry_day x duration grid" section of 06_brazil_travel_final_finite.R
# -- already a per-cell (state x age x duration x entry_day) summary across the
# 1000 PSA draws (median BRR + Pr(BRR>1)), not draw-level, so states are shown
# as separate panels rather than pooled by setting (pooling would need
# draw-level data, which this grid deliberately doesn't retain -- see 06's own
# comments on why).
#
# A dashed contour at BRR = 1 is overlaid on every panel: everything inside
# that contour (toward the fill scale's dark end) is where vaccinating before
# a trip of that duration, entered that week, has median benefit > risk.
#
# Run AFTER 06_brazil_travel_final_finite.R has finished with ar_source = "S0"
# (02b_setup_ar_by_state.R) -- i.e. after the 2026-08-13 rerun.
#
# Output: 06_Results/travel_duration_entryday_{sae,death,daly}_{18-64,65plus}.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")

if (!file.exists("01_Data/psa_grid_bra_travel_finite.RData")) {
  stop("01_Data/psa_grid_bra_travel_finite.RData not found -- run 06_brazil_travel_final_finite.R first.")
}

library(dplyr)
library(ggplot2)
library(scales)

load("01_Data/psa_grid_bra_travel_finite.RData")  # -> psa_grid_travel

# entry_day (1-364, weekly grid) -> calendar week, easier to read than day-of-year
psa_grid_travel <- psa_grid_travel %>%
  mutate(entry_week = (entry_day - 1) %/% 7 + 1)

# duration_grid (06's grid loop) is NOT evenly spaced -- mostly 7-day steps
# but with 30/90/180 forced in on top (so they line up with the discrete-
# duration table), which breaks geom_raster's even-tile assumption and shows
# up as thin comb-like stripes. Linearly interpolate each
# state x age x entry_week row onto a genuinely uniform 1-day duration axis
# (on the log10(BRR) scale, since BRR is plotted log-scaled) before rasterising.
interp_duration <- function(df, xout) {
  ok <- is.finite(df$duration) & is.finite(df$log_brr)
  if (sum(ok) < 2) return(tibble::tibble(duration = xout, log_brr = NA_real_))
  fit <- approx(x = df$duration[ok], y = df$log_brr[ok], xout = xout, rule = 2)
  tibble::tibble(duration = fit$x, log_brr = fit$y)
}

theme_hm <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
}

# Diverging, log-scaled fill centred on BRR = 1 (benefit = risk). Convention
# matched to the outbreak-response weeksweep heatmap (finite_weeksweep_heatmap.R):
# BRR > 1 (benefit > risk) = blue, BRR < 1 (risk > benefit) = red.
brr_fill_scale <- function(med_col) {
  scale_fill_gradient2(
    name = "Median BRR",
    trans = "log10",
    low = "#D73027", mid = "grey92", high = "#4575B4",
    midpoint = 0,  # log10(1) = 0
    labels = label_number(accuracy = 0.1),
    na.value = "grey80"
  )
}

build_heatmap <- function(outcome_label, med_col, age_sel, age_label) {
  duration_uniform <- seq(min(psa_grid_travel$duration), max(psa_grid_travel$duration), by = 1)

  d_raw <- psa_grid_travel %>%
    filter(age_group == age_sel) %>%
    mutate(brr_capped = pmin(pmax(.data[[med_col]], 0.01), 100),  # guard log-scale/contour against 0 or Inf
           log_brr = log10(brr_capped))

  d <- d_raw %>%
    group_by(state, entry_week) %>%
    group_modify(~ interp_duration(.x, duration_uniform)) %>%
    ungroup() %>%
    mutate(brr_capped = 10^log_brr)

  p <- ggplot(d, aes(x = duration, y = entry_week, fill = brr_capped)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = brr_capped), breaks = 1, colour = "black", linewidth = 0.4, linetype = "dashed") +
    facet_wrap(~state, ncol = 4) +
    brr_fill_scale(med_col) +
    scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
    scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 52, 13)) +
    labs(
      title = paste0("Travel vaccination — ", outcome_label, " BRR by duration and entry timing (", age_label, ")"),
      caption = paste0(
        "Fill = median BRR across 1000 PSA draws (log scale, blue = benefit > risk). ",
        "Dashed contour = BRR = 1. Finite-history baseline immunity, S0-denominator attack rate."
      )
    ) +
    theme_hm()

  out_file <- sprintf("06_Results/travel_duration_entryday_%s_%s.png",
                       tolower(gsub("[^A-Za-z]", "", outcome_label)),
                       gsub("[^A-Za-z0-9]", "", age_label))
  ggsave(out_file, p, width = 14, height = 9, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p)
}

outcomes <- list(
  SAE   = "brr_sae_median",
  Death = "brr_death_median",
  DALY  = "brr_daly_median"
)
ages <- list("18-64" = "18-64", "65plus" = "65+")

for (oc_label in names(outcomes)) {
  for (age_label in names(ages)) {
    # Death x 18-64 skipped: vaccine-attributable death risk for 18-64 is
    # fixed at exactly 0 (chat record 2026-08-27), so brr_death_median is NA
    # for every cell in this slice (not "very small" -- undefined, division
    # by an exact-0 denominator) -- geom_contour() errors on an all-NA z
    # field. Same exclusion already applied to Death/18-64 everywhere else
    # in this project (e.g. 10_streamline_mechanism_figure.R's outcome
    # tables); this is not a new gap, just the first heatmap script to hit it.
    if (oc_label == "Death" && age_label == "18-64") next
    build_heatmap(oc_label, outcomes[[oc_label]], ages[[age_label]], age_label)
  }
}
