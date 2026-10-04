# =============================================================================
# 11_travel_headline_figures.R
#
# Builds the travel-vaccination-scenario headline benefit-risk plane figures,
# in the same visual format as the outbreak-response headline figures
# (10_streamline_mechanism_figure.R's headline_benefit_risk_plane_*.png):
# x = risk, y = benefit, both per 10,000 vaccinated (log-log), diagonal =
# BRR 1, colour = age group, shape = transmission setting (Low/Moderate/
# High), one panel per travel duration (7/14/30/90 days) instead of
# mechanism (the travel scenario applies a single vaccine efficacy, not a
# disease- vs infection-blocking split).
#
# Run AFTER 06_brazil_travel_final_finite.R has finished and saved
# 01_Data/psa_df_bra_travel_finite.RData.
#
# Setting classification reused from 01_Data/setting_key.RData (the
# case-rate-based Low/Moderate/High classification finalised in
# 08_attack_rate_classification_finite.R this session), so travel and
# outbreak-response results are stratified by the exact same state groups.
#
# Output: 06_Results/headline_benefit_risk_plane_travel_{DALY,Death,SAE}.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")

if (!file.exists("01_Data/psa_df_bra_travel_finite.RData")) {
  stop(
    "01_Data/psa_df_bra_travel_finite.RData not found. ",
    "Run 02_Scripts/06_brazil_travel_final_finite.R first (~1000-draw PSA, ",
    "takes roughly 1.5-2 hours)."
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

load("01_Data/psa_df_bra_travel_finite.RData")   # -> psa_df
load("01_Data/setting_key.RData")                # -> setting_key (case-rate-based)

# ---- 0) Shared styling (mirrors 10_streamline_mechanism_figure.R) ----------
theme_nm <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
}
log_num_labels <- scales::label_number(accuracy = NULL, big.mark = ",")
# Distinct hues (not a same-hue light->dark ramp -- those read as "all blue"
# at thin line widths) for the ordinal Low/Moderate/High setting tiers.
pal_setting <- c(Low = "#2A9D8F", Moderate = "#F4A261", High = "#C44E52")

# ---- 1) Attach setting, tidy to long format (outcome x risk/benefit) -------
psa_df <- psa_df %>%
  mutate(setting = unname(setting_key[state])) %>%
  filter(!is.na(setting)) %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    days = factor(days, levels = c("7d", "14d", "30d", "90d")),
    age_group = factor(age_group, levels = c("18-64", "65+"))
  )

# One row per draw x state x age x days x outcome, with a common
# benefit/risk/brr naming scheme.
travel_long <- bind_rows(
  psa_df %>% transmute(draw, state, setting, age_group, days, entry_day,
                        outcome = "DALY", benefit = daly_averted, risk = daly_sae, brr = brr_daly),
  psa_df %>% transmute(draw, state, setting, age_group, days, entry_day,
                        outcome = "SAE", benefit = averted_10k_sae, risk = excess_10k_sae, brr = brr_sae),
  psa_df %>% transmute(draw, state, setting, age_group, days, entry_day,
                        outcome = "Death", benefit = averted_10k_death, risk = excess_10k_death, brr = brr_death)
)

# ---- 2) Pool draws within each setting (same convention as the outbreak- ---
# response headline plane: brr_draw_summary_true / build_setting_synthesis_data)
build_travel_synthesis_data <- function(outcome_sel) {
  travel_long %>%
    filter(outcome == outcome_sel) %>%
    group_by(setting, age_group, days) %>%
    summarise(
      benefit_med = median(benefit, na.rm = TRUE),
      benefit_lo  = quantile(benefit, 0.025, na.rm = TRUE),
      benefit_hi  = quantile(benefit, 0.975, na.rm = TRUE),
      risk_med = median(risk, na.rm = TRUE),
      risk_lo  = quantile(risk, 0.025, na.rm = TRUE),
      risk_hi  = quantile(risk, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Floor at a small epsilon so log10(0) doesn't break the axis/CI bars
    # for outcomes (e.g. Death) whose 2.5th percentile can be exactly 0 --
    # same reasoning as the outbreak-response headline plane.
    mutate(across(c(benefit_med, benefit_lo, benefit_hi, risk_med, risk_lo, risk_hi),
                   ~ pmax(.x, 1e-4)))
}

# ---- 3) Benefit-risk plane, one outcome at a time --------------------------
build_travel_brplane <- function(outcome_sel) {
  d <- build_travel_synthesis_data(outcome_sel)

  rng <- pmax(range(c(d$risk_med, d$benefit_med, d$risk_lo, d$risk_hi,
                       d$benefit_lo, d$benefit_hi), na.rm = TRUE), 1e-3)

  # No position_dodge(): risk_med is driven by p_sae_vacc/p_death_vacc (age-
  # dependent only), so it's identical across settings within an age group --
  # dodging the x-position would misleadingly draw points away from their
  # true risk value (same bug/fix as headline_benefit_risk_plane_*.png in
  # 10_streamline_mechanism_figure.R).
  p <- ggplot(d, aes(risk_med, benefit_med, colour = age_group, shape = setting)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(ymin = benefit_lo, ymax = benefit_hi), width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = risk_lo, xmax = risk_hi), height = 0, alpha = 0.3) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~days, nrow = 1) +
    coord_equal(xlim = rng, ylim = rng) +
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"), name = "Age group") +
    scale_shape_manual(values = c(Low = 15, Moderate = 17, High = 16), name = "Setting") +
    scale_x_log10(labels = log_num_labels, name = paste(outcome_sel, "attributable risk per 10,000 vaccinated (log scale)")) +
    scale_y_log10(labels = log_num_labels, name = paste(outcome_sel, "averted per 10,000 vaccinated (log scale)")) +
    labs(
      title = paste0("Travel vaccination — ", outcome_sel, " benefit-risk plane, by travel duration"),
      caption = "Points below the dashed line: benefit > risk (BRR > 1). Points above: risk > benefit (BRR < 1).\nOne point per Setting x Age group x Duration (draws pooled across states within each setting, same setting classification as the outbreak-response results)."
    ) +
    theme_nm()

  out_file <- sprintf("06_Results/headline_benefit_risk_plane_travel_%s.png", outcome_sel)
  ggsave(out_file, p, width = 14, height = 5.5, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p)
}

for (oc in c("DALY", "Death", "SAE")) build_travel_brplane(oc)

# ---- 3b) Net benefit-or-risk bar chart -- same visual grammar as the ---
# outbreak-response net_benefit_risk_*.png (10_streamline_mechanism_figure.R):
# horizontal bars showing benefit MINUS risk (absolute difference per 10,000
# vaccinated, not the ratio), blue = net benefit, red = net risk, one row per
# Setting x Age group. Travel has no mechanism/serostatus-adjustment split
# (single VE assumption, no base-vs-adjusted risk computed here), so duration
# takes the place mechanism/serostatus occupy in the outbreak-response
# version -- it's the axis travel actually varies.
build_travel_net_benefit_data <- function(outcome_sel) {
  travel_long %>%
    filter(outcome == outcome_sel) %>%
    mutate(net_10k = benefit - risk) %>%
    group_by(setting, age_group, days) %>%
    summarise(net_med = median(net_10k, na.rm = TRUE),
              net_lo  = quantile(net_10k, 0.025, na.rm = TRUE),
              net_hi  = quantile(net_10k, 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      row_label = factor(paste0(setting, " | ", age_group),
                          levels = rev(as.vector(outer(c("Low", "Moderate", "High"),
                                                        c("18-64", "65+"), paste, sep = " | ")))),
      sign = ifelse(net_med >= 0, "Net benefit", "Net risk")
    )
}

build_and_save_travel_net_benefit_chart <- function() {
  d <- bind_rows(
    build_travel_net_benefit_data("DALY")  %>% mutate(outcome = "DALY"),
    build_travel_net_benefit_data("SAE")   %>% mutate(outcome = "SAE"),
    build_travel_net_benefit_data("Death") %>% mutate(outcome = "Death")
  ) %>%
    mutate(outcome = factor(outcome, levels = c("DALY", "SAE", "Death"))) %>%
    filter(!is.na(row_label))

  p <- ggplot(d, aes(row_label, net_med, fill = sign)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = net_lo, ymax = net_hi), width = 0.2, colour = "grey30", linewidth = 0.4) +
    coord_flip() +
    facet_grid(days ~ outcome, scales = "free_x") +
    scale_fill_manual(values = c("Net benefit" = "#274690", "Net risk" = "#C44E52"), name = NULL) +
    scale_y_continuous(labels = log_num_labels) +
    labs(x = NULL, y = "Difference in incidence per 10,000 vaccinated (benefit − risk)"
         ) +
    theme_nm() +
    theme(legend.position = "bottom")

  out_file <- "06_Results/net_benefit_risk_travel.png"
  ggsave(out_file, p, width = 11, height = 12, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p)
}

build_and_save_travel_net_benefit_chart()

# ---- 3c) Continuous net benefit-or-risk vs travel duration, with 95% UI ----
# The bar chart above only shows 4 discrete durations (7d/14d/30d/90d) as
# separate bars. This uses the SAME draw-level data (travel_long, which has
# real per-draw benefit/risk at those 4 durations, so real quantiles are
# available) but draws it as a smooth curve with an uncertainty ribbon:
# median/95%UI are computed at the 4 simulated durations, then natural-spline
# interpolated to a dense duration sequence purely for a smooth visual --- no
# new data are simulated. NOTE: this deliberately does NOT use the fine
# entry_day x duration grid (01_Data/psa_grid_bra_travel_finite_setting.RData,
# 7-180 days) built in 06_brazil_travel_final_finite.R, because that grid
# only stored per-cell MEDIANS (draws discarded to save memory -- see that
# script's comment) -- there's no quantile information to build a ribbon
# from at those finer duration points. So the curve here only extends to 90
# days (the largest duration actually simulated with retained draws);
# extrapolating the spline past 90d would not be supported by real data.
build_travel_net_benefit_smooth_data <- function(outcome_sel) {
  travel_long %>%
    filter(outcome == outcome_sel) %>%
    mutate(net_10k = benefit - risk,
           duration = as.integer(sub("d$", "", as.character(days)))) %>%
    group_by(setting, age_group, duration) %>%
    summarise(net_med = median(net_10k, na.rm = TRUE),
              net_lo  = quantile(net_10k, 0.025, na.rm = TRUE),
              net_hi  = quantile(net_10k, 0.975, na.rm = TRUE), .groups = "drop")
}

duration_interp_seq <- seq(7, 90, by = 1)

travel_net_smooth <- bind_rows(
  build_travel_net_benefit_smooth_data("DALY")  %>% mutate(outcome = "DALY"),
  build_travel_net_benefit_smooth_data("SAE")   %>% mutate(outcome = "SAE"),
  build_travel_net_benefit_smooth_data("Death") %>% mutate(outcome = "Death")
) %>%
  mutate(outcome = factor(outcome, levels = c("DALY", "SAE", "Death"))) %>%
  group_by(outcome, setting, age_group) %>%
  group_modify(~ tibble(
    duration = duration_interp_seq,
    net_med  = stats::spline(.x$duration, .x$net_med, xout = duration_interp_seq, method = "natural")$y,
    net_lo   = stats::spline(.x$duration, .x$net_lo,  xout = duration_interp_seq, method = "natural")$y,
    net_hi   = stats::spline(.x$duration, .x$net_hi,  xout = duration_interp_seq, method = "natural")$y
  )) %>%
  ungroup()

p_net_duration <- ggplot(travel_net_smooth, aes(duration, net_med, colour = setting, fill = setting)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = net_lo, ymax = net_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(age_group ~ outcome, scales = "free_y") +
  scale_colour_manual(values = pal_setting, name = "Setting") +
  scale_fill_manual(values = pal_setting, name = "Setting") +
  labs(x = "Travel duration (days)",
       y = "Net benefit − risk per 10,000 vaccinated"
       ) +
  theme_nm()

ggsave("06_Results/net_benefit_risk_travel_by_duration.png", p_net_duration,
       width = 11, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/net_benefit_risk_travel_by_duration.png")

# ---- 3d) Continuous net benefit-or-risk vs entry timing --------------------
# Same idea as 3c but for entry timing instead of duration: psa_df already
# carries entry_day per draw (median of the within-draw sampled entry days,
# 06_brazil_travel_final_finite.R line ~265) as a genuine per-draw quantity,
# so geom_smooth()'s loess + SE ribbon can be used directly on travel_long
# without any manual spline step. One file per duration (entry_day's valid
# range depends on duration -- a 90-day trip can't start after day ~275 --
# so pooling durations together would mix truncated ranges); within each
# file, facet by outcome x age group, colour by setting.
build_and_save_travel_net_benefit_by_entry <- function(days_label) {
  d <- travel_long %>%
    filter(days == days_label) %>%
    mutate(net_10k = benefit - risk,
           outcome = factor(outcome, levels = c("DALY", "SAE", "Death")))

  p <- ggplot(d, aes(entry_day, net_10k, colour = setting, fill = setting)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.9) +
    facet_grid(age_group ~ outcome, scales = "free_y") +
    scale_colour_manual(values = pal_setting, name = "Setting") +
    scale_fill_manual(values = pal_setting, name = "Setting") +
    labs(x = "Entry day (day of year)",
         y = "Net benefit − risk per 10,000 vaccinated",
         title = paste0("Travel vaccination — net benefit/risk vs entry timing (", days_label, " trips)"),
         caption = "Line = loess fit, ribbon = SE, across draw-level entry days (one median entry day per draw). Pooled across states within each setting.") +
    theme_nm()

  out_file <- sprintf("06_Results/net_benefit_risk_travel_by_entry_%s.png", days_label)
  ggsave(out_file, p, width = 11, height = 7, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p)
}

for (dl in c("7d", "14d", "30d", "90d")) build_and_save_travel_net_benefit_by_entry(dl)

# ---- 4) Compact summary table (median BRR + Pr(BRR>1), by setting x age x --
# duration x outcome) -- same spirit as brr_table_final_long on the
# outbreak-response side.
travel_summary_table <- travel_long %>%
  group_by(outcome, setting, age_group, days) %>%
  summarise(
    BRR_med = median(brr, na.rm = TRUE),
    BRR_lo  = quantile(brr, 0.025, na.rm = TRUE),
    BRR_hi  = quantile(brr, 0.975, na.rm = TRUE),
    prob_gt1 = mean(brr > 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(outcome, setting, age_group, days)

print(travel_summary_table, n = 50)

ft <- flextable::flextable(
  travel_summary_table %>%
    mutate(
      `Median BRR (95% UI)` = sprintf("%.2f (%.2f\u2013%.2f)", BRR_med, BRR_lo, BRR_hi),
      `Pr(BRR>1)` = sprintf("%.0f%%", 100 * prob_gt1)
    ) %>%
    select(Outcome = outcome, Setting = setting, `Age group` = age_group,
           Duration = days, `Median BRR (95% UI)`, `Pr(BRR>1)`)
) %>%
  flextable::set_caption(caption = paste(
    "Travel vaccination benefit-risk ratio (BRR) by outcome, setting, age group, and travel duration.",
    "Setting = case-rate-based Low/Moderate/High classification (same as outbreak-response results).",
    "Draws pooled across states within each setting; 1000 PSA draws."
  )) %>%
  flextable::merge_v(j = c("Outcome", "Setting", "Age group")) %>%
  flextable::valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
  flextable::autofit() %>%
  flextable::theme_vanilla()

doc <- officer::read_docx() %>%
  flextable::body_add_flextable(ft)
print(doc, target = "06_Results/travel_brr_summary_table.docx")
message("Saved: 06_Results/travel_brr_summary_table.docx")
