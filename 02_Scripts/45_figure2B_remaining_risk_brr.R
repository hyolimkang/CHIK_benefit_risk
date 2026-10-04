# =============================================================================
# 45_figure2B_remaining_risk_brr.R
#
# Figure 2B: Outbreak-response benefit-risk across REMAINING epidemic risk
# (chat record 2026-09-03) -- reframes campaign timing not as calendar week,
# but as "how much of the outbreak's total infection risk is still ahead"
# at the moment the campaign starts. This is the epidemiologically meaningful
# quantity a vaccine can actually act on (only future infections), unlike
# calendar week which is state-specific and not directly comparable across
# settings with different outbreak curves.
#
#   x = remaining symptomatic-infection risk after campaign initiation (%)
#       = 1 - exp(-sum(posterior lambda[week:end])), i.e. the cumulative
#         hazard NOT YET accrued at the campaign's start week, converted to
#         an attack rate. Same identity already used/validated for full-year
#         AR elsewhere in this project (16_figure1_epi_map_ar.R,
#         32_figure1_state_summary_table.R: rowSums(lambda) reproduces the
#         known full-year AR, e.g. Ceará ~6.5%) -- this is just a PARTIAL
#         (tail) sum starting at the campaign week instead of the full sum.
#   y = median DALY BRR (log scale), reference coverage 50% (interpolated --
#       the week-sweep grid only has cov 10/30/60/90, no exact 50% point;
#       linear interpolation between the two bracketing coverage levels,
#       no extrapolation).
#   colour = protection mechanism (DB / D+I)
#   linetype = base vs serostatus-adjusted risk
#   point shape = transmission setting (Low/Moderate/High, 2022 intensity)
#   columns = age group (18-64 / 65+)
#
# Inputs:
#   - CHIK_ORV_impact/00_Data/0_2_Processed/posterior_finite_all.RData
#     (posterior_ce$lambda etc. -- weekly hazard, 1000 draws x 54 weeks)
#   - 01_Data/setting_key.RData (state -> Low/Moderate/High)
#   - CHIK_ORV_impact/00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep{,_ve0}.RData
#     (week x coverage BRR grid, both mechanisms)
#
# Output: 06_Results/figure2B_remaining_risk_vs_brr.png
#         07_Final_Results/Main_Figures/figure2B_remaining_risk_vs_brr.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

ori_root <- "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"

# ---- 1) Remaining infection risk by state x campaign-start-week -----------
WEEKS_TESTED <- c(1, 8, 16, 24, 32, 42)  # matches the week-sweep grid

message("Loading posterior_finite_all.RData ...")
load(file.path(ori_root, "00_Data/0_2_Processed/posterior_finite_all.RData"))
load("01_Data/setting_key.RData")  # -> setting_key

posterior_list <- list(
  "Ceará" = posterior_ce, "Bahia" = posterior_bh, "Paraíba" = posterior_pa,
  "Pernambuco" = posterior_pn, "Rio Grande do Norte" = posterior_rg,
  "Piauí" = posterior_pi, "Alagoas" = posterior_ag, "Tocantins" = posterior_tc,
  "Minas Gerais" = posterior_mg, "Sergipe" = posterior_se, "Goiás" = posterior_go
)

remaining_ar_draws <- purrr::imap_dfr(posterior_list, function(post, region) {
  n_weeks <- ncol(post$lambda)
  purrr::map_dfr(WEEKS_TESTED, function(w) {
    remaining_hazard <- rowSums(post$lambda[, w:n_weeks, drop = FALSE])
    tibble::tibble(
      region = region,
      setting = setting_key[[region]],
      week = w,
      remaining_ar = 1 - exp(-remaining_hazard)
    )
  })
})

# Pool draws across states within a setting (same convention as
# 21_foi_risk_brr_peak_timing.R's setting-level pooling) -- one combined
# distribution per (setting, week), not an average of per-state medians.
remaining_ar_summary <- remaining_ar_draws %>%
  group_by(setting, week) %>%
  summarise(remaining_ar_med = median(remaining_ar, na.rm = TRUE), .groups = "drop") %>%
  mutate(remaining_pct = 100 * remaining_ar_med)

print(as.data.frame(remaining_ar_summary %>% arrange(setting, week)))

# ---- 2) BRR at each (setting, age, mechanism, week), interpolated to cov50 --
load(file.path(ori_root, "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep_ve0.RData"))
weeksweep_ve0 <- all_weeks_brr
load(file.path(ori_root, "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep.RData"))
weeksweep_ve989 <- all_weeks_brr
rm(all_weeks_brr)

weeksweep_summary <- bind_rows(weeksweep_ve0, weeksweep_ve989) %>%
  filter(outcome == "DALY", RR_seropos == 0, AgeCat %in% c("18-64", "65+"),
         week %in% WEEKS_TESTED) %>%
  group_by(setting, AgeCat, VE_label, week, Coverage) %>%
  summarise(
    brr_base_med = median(brr_base, na.rm = TRUE),
    brr_adj_med  = median(brr_adj,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(coverage_pct = as.numeric(sub("cov", "", Coverage)))

# Linear interpolation to coverage = 50% (within the observed grid only;
# no extrapolation -- rule = 1 gives NA outside the covered range, which
# can't happen here since 50 sits between the tested 30/60).
interpolate_to_cov50 <- function(data, value_col) {
  data %>%
    group_by(setting, AgeCat, VE_label, week) %>%
    group_modify(~ {
      x <- .x$coverage_pct
      y <- .x[[value_col]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 2) return(tibble(BRR = NA_real_))
      tibble(BRR = approx(x[ok], y[ok], xout = 50, rule = 1)$y)
    }) %>%
    ungroup()
}

brr_base_50 <- interpolate_to_cov50(weeksweep_summary, "brr_base_med") %>% mutate(risk_assumption = "Base risk")
brr_adj_50  <- interpolate_to_cov50(weeksweep_summary, "brr_adj_med")  %>% mutate(risk_assumption = "Serostatus-adjusted risk")

brr_50 <- bind_rows(brr_base_50, brr_adj_50) %>%
  mutate(
    risk_assumption = factor(risk_assumption, levels = c("Base risk", "Serostatus-adjusted risk")),
    mechanism = factor(VE_label, levels = c("Disease blocking only", "Disease and infection blocking")),
    age_group = factor(AgeCat, levels = c("18-64", "65+"), labels = c("18-64", "\u226565")),
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    BRR = ifelse(is.finite(BRR) & BRR > 0, BRR, NA_real_)
  )

# ---- 3) Join remaining infection risk onto the BRR grid --------------------
plot_df <- brr_50 %>%
  left_join(remaining_ar_summary %>% mutate(setting = factor(setting, levels = c("Low", "Moderate", "High"))),
            by = c("setting", "week"))

# ---- 4) Plot -----------------------------------------------------------------
# Redesign (chat record 2026-09-03): setting moves from point-shape to its own
# facet column (2x3 grid: rows = age group, columns = setting). With setting
# already split out as a facet, encoding it AGAIN via point shape was
# redundant/cluttered -- dropped entirely. Each panel now holds exactly 4
# lines (mechanism x risk assumption), each a clean monotonic curve in its
# own right (no more cross-setting zigzag to work around).
mechanism_cols <- c("Disease blocking only" = "#1F9E89", "Disease and infection blocking" = "#C44E52")

p <- ggplot(plot_df, aes(x = remaining_pct, y = BRR, colour = mechanism, linetype = risk_assumption,
                          group = interaction(mechanism, risk_assumption))) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  geom_path(linewidth = 0.9, lineend = "round", na.rm = TRUE) +
  geom_point(size = 1.5, alpha = 0.75, na.rm = TRUE) +
  scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
  scale_linetype_manual(values = c("Base risk" = "solid", "Serostatus-adjusted risk" = "dashed"), name = "Risk assumption") +
  scale_x_continuous(name = "Remaining infection risk after campaign initiation (%)",
                      labels = scales::label_number(suffix = "%", accuracy = 0.1),
                      breaks = scales::breaks_pretty(4),
                      expand = expansion(mult = 0.08)) +
  scale_y_log10(name = "Median DALY benefit-risk ratio (BRR)",
                breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")) +
  facet_grid(age_group ~ setting) +
  labs(title = "B. Outbreak-response benefit–risk across remaining epidemic risk") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey40", linewidth = 0.3),
    axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12.5),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom",
    panel.spacing = unit(12, "pt")
  )

# ---- Endpoint labels, ONE panel only (chat record 2026-09-03) --------------
# Earlier version repeated "Earlier/Later campaign start" italic text in all
# 6 (age x setting) panels; a figure-wide arrow subtitle was tried as a
# replacement but the user wants the original in-panel endpoint-label style
# back, just shown once (top-left panel: 18-64 x Low) instead of 6x.
endpoint_labels <- plot_df %>%
  filter(week %in% c(1, 42), age_group == "18-64", setting == "Low") %>%
  group_by(age_group, setting, week) %>%
  summarise(remaining_pct = first(remaining_pct), brr_at_week = median(BRR, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    label = ifelse(week == 1, "Earlier campaign start", "Later campaign start"),
    # Low panel's week-1 point sits close to the y-axis, so growing the text
    # leftward (hjust=1) clips it against the axis -- grow rightward instead.
    hjust = -0.05,
    y = ifelse(week == 1, 10^(log10(brr_at_week) + 0.32), 10^(log10(brr_at_week) - 0.32))
  )

p <- p +
  geom_text(data = endpoint_labels, aes(x = remaining_pct, y = y, label = label, hjust = hjust),
            inherit.aes = FALSE, size = 2.9, fontface = "italic", colour = "grey25")

dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)

for (out_file in c("06_Results/figure2B_remaining_risk_vs_brr.png",
                    "07_Final_Results/Main_Figures/figure2B_remaining_risk_vs_brr.png")) {
  ggsave(out_file, p, width = 11, height = 6.5, dpi = 300, bg = "white")
  message("Saved: ", out_file)
}
