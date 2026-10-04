# =============================================================================
# 15_travel_duration_entryday_probgt1_heatmap.R
#
# Uncertainty companion to 13_travel_duration_entryday_heatmap_setting.R:
# same x/y/facet layout (duration x entry-week, age x setting), but fill =
# Pr(BRR > 1) across posterior draws at each cell, instead of median BRR.
#
# The median-BRR heatmap shows WHERE the surface favours vaccination; this
# one shows HOW CONFIDENT that is -- a cell can have a favourable median BRR
# while still being close to a coin flip if its 95% UI is wide, and that
# distinction is invisible in the median-only figure. Read the two heatmaps
# side by side (see chat record 2026-08-18).
#
# Uses brr_*_prob_gt1 columns already saved in psa_grid_bra_travel_finite_setting.RData
# by 06_brazil_travel_final_finite.R's setting-level pooling block -- no
# rerun of the simulation needed.
#
# Output: 06_Results/travel_duration_entryday_setting_probgt1_{sae,death,daly}.png
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

# Same unevenly-spaced duration_grid issue as the median heatmap -- interpolate
# onto a uniform 1-day duration axis before rasterising. Linear on the raw
# [0,1] probability scale here (not log-space -- that was specifically a BRR
# ratio-scale fix, meaningless for an already-bounded probability).
interp_duration <- function(df, xout) {
  ok <- is.finite(df$duration) & is.finite(df$prob_gt1)
  if (sum(ok) < 2) return(tibble::tibble(duration = xout, prob_gt1 = NA_real_))
  fit <- approx(x = df$duration[ok], y = df$prob_gt1[ok], xout = xout, rule = 2)
  tibble::tibble(duration = fit$x, prob_gt1 = pmin(pmax(fit$y, 0), 1))
}

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

prob_fill_scale <- function() {
  # Same convention as the median heatmaps: red = more likely unfavourable,
  # blue = more likely favourable, neutral grey at 50% (as likely as not).
  scale_fill_gradient2(
    name = "Pr(BRR>1)",
    low = pole_red, mid = mid_grey, high = pole_blue,
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::percent_format(accuracy = 1),
    na.value = ink_gridline
  )
}

build_setting_prob_heatmap <- function(outcome_label, prob_col) {
  duration_uniform <- seq(min(psa_grid_travel_setting$duration), max(psa_grid_travel_setting$duration), by = 1)

  d_raw <- psa_grid_travel_setting %>%
    mutate(prob_gt1 = .data[[prob_col]])

  d <- d_raw %>%
    group_by(setting, age_group, entry_week) %>%
    group_modify(~ interp_duration(.x, duration_uniform)) %>%
    ungroup()

  p <- ggplot(d, aes(x = duration, y = entry_week, fill = prob_gt1)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = prob_gt1), breaks = 0.5, colour = ink_primary, linewidth = 0.35, linetype = "42") +
    facet_grid(age_group ~ setting) +
    prob_fill_scale() +
    scale_x_continuous(name = "Travel duration (days)", expand = c(0, 0)) +
    scale_y_continuous(name = "Entry week of year", expand = c(0, 0), breaks = seq(0, 52, 13)) +
    labs(
      caption = paste0(
        outcome_label, " — Pr(BRR>1) across posterior draws, states pooled with draws within setting.\n",
        "Dashed line = 50% (as likely as not); blue side = more likely benefit > risk. Companion to the median-BRR heatmap."
      )
    ) +
    theme_hm()

  out_file <- sprintf("06_Results/travel_duration_entryday_setting_probgt1_%s.png",
                       tolower(gsub("[^A-Za-z]", "", outcome_label)))
  ggsave(out_file, p, width = 9.5, height = 5.9, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
  invisible(p)
}

outcomes <- list(SAE = "brr_sae_prob_gt1", Death = "brr_death_prob_gt1", DALY = "brr_daly_prob_gt1")
for (oc_label in names(outcomes)) build_setting_prob_heatmap(oc_label, outcomes[[oc_label]])
