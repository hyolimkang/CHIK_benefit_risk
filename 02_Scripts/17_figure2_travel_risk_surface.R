# =============================================================================
# 17_figure2_travel_risk_surface.R
#
# Figure 2. Travel risk-benefit surfaces.
#   Panel A: Infection-risk surface. P(infection) = 1 - exp(-H_trip), where
#            H_trip is the cumulative force of infection accrued during the
#            trip (state-level daily FOI integrated over the entry_day ->
#            entry_day+duration window). x = travel duration (days), y =
#            entry week (week of year), one panel per transmission setting
#            (Low/Moderate/High), pooled across states within each setting.
#            NOT age-specific -- infection risk (unlike SAE/death/DALY risk
#            downstream) doesn't depend on age in this model: H_trip and
#            AR_travel are computed from state-level FOI alone, before the
#            age-specific SAE/death rates are ever applied
#            (06_brazil_travel_final_finite.R lines ~394-410: `cs`, `H_full`,
#            `AR_travel` are all computed OUTSIDE the age_idx loop, then just
#            replicated across age_levels when writing into the grid array).
#            Verified empirically: AR_median is bit-identical across
#            age_group for every (setting, duration, entry_day) cell.
#   Panel B: Pr(DALY BRR > 1) across posterior draws, same duration x
#            entry-week surface, faceted by age group x setting -- this is
#            15_travel_duration_entryday_probgt1_heatmap.R's DALY panel,
#            reproduced here (not sourced, to keep this script self-contained)
#            so it sits directly under panel A as one composite figure. Style
#            (ink-colour palette, geom_raster+interpolate, theme_hm) matches
#            that script exactly for visual consistency between A and B.
#
# Duration + entry week combinations that would run past the end of the
# observed 2022 FOI series (e.g. a 180-day trip starting in week 40) aren't
# computable and render as true grey/NA (interp_duration() uses
# approx(..., rule=1): NA outside the range of durations actually simulated
# for that entry_week -- NOT rule=2/flat-extrapolated, which would fabricate
# a value for a combination that was never simulated).
#
# Data source: 01_Data/psa_grid_bra_travel_finite_setting.RData
#   (psa_grid_travel_setting; fine entry_day x duration grid, built in
#   06_brazil_travel_final_finite.R, pooled by setting).
#
# Output: 06_Results/figure2_travel_risk_surfaces.png
# =============================================================================

setwd(getwd())
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

if (!file.exists("01_Data/psa_grid_bra_travel_finite_setting.RData")) {
  stop(
    "01_Data/psa_grid_bra_travel_finite_setting.RData not found. ",
    "Run the fine entry_day x duration grid block in 06_brazil_travel_final_finite.R first."
  )
}
load("01_Data/psa_grid_bra_travel_finite_setting.RData")  # -> psa_grid_travel_setting

psa_grid_travel_setting <- psa_grid_travel_setting %>%
  mutate(entry_week = (entry_day - 1) %/% 7 + 1,
         age_group  = factor(age_group, levels = c("18-64", "65+"),
                             labels = c("18-64", "\u226565")),
         setting    = factor(setting, levels = c("Low", "Moderate", "High")))

# ---- 0) Shared styling (matches 15_travel_duration_entryday_probgt1_heatmap.R) ----
ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
ink_muted     <- "#898781"
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

# Unevenly-spaced duration_grid (mostly 7-day steps, with 30/90/180 forced
# in -- see 06_brazil_travel_final_finite.R) -- interpolate onto a uniform
# 1-day duration axis before rasterising, per (setting, [age_group,] y-group).
# rule = 1 (NOT 2): xout values outside the range of durations actually
# computable for that entry_week (i.e. the trip would run past the observed
# 2022 FOI series) return NA and render as true grey -- rule = 2 would
# flat-extrapolate a value that was never simulated, which is fabricating
# data for a combination that's genuinely unknown, not just unsmoothed.
interp_duration <- function(df, yvar, xout) {
  ok <- is.finite(df$duration) & is.finite(df[[yvar]])
  if (sum(ok) < 2) return(tibble::tibble(duration = xout, value = NA_real_))
  fit <- approx(x = df$duration[ok], y = df[[yvar]][ok], xout = xout, rule = 1)
  tibble::tibble(duration = fit$x, value = fit$y)
}

# Validity is a strict diagonal (entry_day + duration <= 364, the same for
# every state -- phi_pred is always exactly 52 weeks/364 days long) -- no
# rectangular axis crop removes the NA triangle without cutting one of the
# two axes hard. Capping duration <=90d (364 - (40-1)*7 = 91 >= 90, so this
# is the largest duration that stays valid through entry week 40) and entry
# week <=40 keeps the full displayed rectangle grey-free with no fabricated
# values; 90d also matches the discrete duration category used elsewhere in
# this repo's travel figures (7d/14d/30d/90d), so nothing new is introduced.
psa_grid_travel_setting <- psa_grid_travel_setting %>%
  filter(duration <= 90, entry_week <= 40)

duration_uniform <- seq(min(psa_grid_travel_setting$duration), max(psa_grid_travel_setting$duration), by = 1)

# ---- A) Panel A -- infection-risk surface ----------------------------------
# age-independent (see header) -- one age_group avoids duplicated cells.
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

# ---- B) Panel B -- Pr(DALY BRR > 1) surface, age x setting -----------------
risk_B <- psa_grid_travel_setting %>%
  mutate(prob_gt1 = brr_daly_prob_gt1) %>%
  group_by(setting, age_group, entry_week) %>%
  group_modify(~ interp_duration(.x, "prob_gt1", duration_uniform)) %>%
  ungroup() %>%
  rename(prob_gt1 = value) %>%
  mutate(prob_gt1 = pmin(pmax(prob_gt1, 0), 1))

p_B <- ggplot(risk_B, aes(duration, entry_week, fill = prob_gt1)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = prob_gt1), breaks = 0.5, colour = ink_primary, linewidth = 0.35, linetype = "42") +
  facet_grid(age_group ~ setting) +
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
  labs(
    title = "B. DALY benefit-risk surface",
  ) +
  theme_hm()

# ---- C) Combine --------------------------------------------------------------
p_fig2 <- (p_A / p_B) +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, colour = ink_primary),
      plot.caption = element_text(size = 8, colour = ink_secondary, hjust = 0)
    )
  )

ggsave("06_Results/figure2_travel_risk_surfaces.png", p_fig2,
       width = 11, height = 11, dpi = 300, bg = "#fcfcfb")
message("Saved: 06_Results/figure2_travel_risk_surfaces.png")
