# =============================================================================
# 34_seromix_seroneg_by_week_mechanism.R
#
# Supplementary line plot for the paper's timing-sensitivity result: for each
# INDEPENDENT campaign-start-week scenario (week 1/8/16/24/32/42; week 52
# excluded -- see below), the proportion of vaccine recipients who were
# seronegative AT THE MOMENT OF VACCINATION (i.e. campaign-wide, dose-weighted
# across the whole rollout), by mechanism (DB vs D+I). This is the same
# quantity behind the manuscript sentence "91.6% of vaccinees were
# seronegative under disease blocking only compared with 79.3% under disease
# and infection blocking" (reference scenario, week 2/cov50) -- this figure
# shows how that number moves across campaign start week instead of reporting
# only the single reference-week value.
#
# Design (chat record 2026-09-02, updated 2026-09-02): Mechanism (DB/D+I) is
# the primary visual comparison -> colour. Facets = Age group (rows) x
# Setting x Coverage (nested columns) -- all three transmission settings are
# shown (not just High), since a diagnostic check found coverage barely moves
# this metric under DB (<0.01pp across cov10-90) but moves it substantially
# under D+I, and MORE so in Low/Moderate than in High (up to ~19pp range
# across cov10-90 in Low/Moderate vs a smaller range in High) -- vaccinating
# deeper into the population measurably changes who's left to be naturally
# infected when infection-blocking immunity is also removing susceptibles,
# and that setting-dependence is exactly the kind of effect a single fixed
# setting would hide.
#
# q_seroneg_vacc isn't saved directly in the week-sweep per-draw export, but
# sae_10k_seroneg = 1e4 * q_seroneg_vacc * p_sae_vacc_base while
# sae_10k_base = 1e4 * p_sae_vacc_base unconditionally, so their ratio
# recovers q_seroneg_vacc exactly (same identity as script 28).
#
# Inputs (CHIK_ORV_impact project -- week-sweep per-draw RData live there):
#   - draw_level_xy_serostatus_finite_weeksweep.RData      (D+I / VE98.9)
#   - draw_level_xy_serostatus_finite_weeksweep_ve0.RData  (DB / VE0)
#
# Output: 06_Results/seromix_seroneg_by_week_mechanism.png
#         06_Results/seromix_seroneg_by_week_mechanism.xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(writexl)
  library(patchwork)
})

ori_root <- "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"

extract_seroneg <- function(rdata_path, mechanism_label) {
  e <- new.env()
  load(rdata_path, envir = e)
  e$all_weeks_brr %>%
    filter(outcome == "DALY", RR_seropos == 0,
           setting %in% c("Low", "Moderate", "High"),
           AgeCat %in% c("18-64", "65+"),
           week != 52) %>%   # truncated rollout -- not a comparable endpoint
    mutate(
      q_seroneg   = sae_10k_seroneg / sae_10k_base,
      seroneg_pct = 100 * q_seroneg
    ) %>%
    group_by(week, Coverage, AgeCat, setting) %>%
    summarise(
      seroneg_med = median(seroneg_pct, na.rm = TRUE),
      seroneg_lo  = quantile(seroneg_pct, 0.025, na.rm = TRUE),
      seroneg_hi  = quantile(seroneg_pct, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Mechanism = mechanism_label)
}

plot_df <- bind_rows(
  extract_seroneg(file.path(ori_root, "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep_ve0.RData"), "Disease blocking only"),
  extract_seroneg(file.path(ori_root, "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep.RData"), "Disease and infection blocking")
) %>%
  mutate(
    Mechanism = factor(Mechanism, levels = c("Disease blocking only", "Disease and infection blocking")),
    AgeCat    = factor(AgeCat, levels = c("18-64", "65+"), labels = c("18-64 years", "65+ years")),
    setting   = factor(setting, levels = c("Low", "Moderate", "High")),
    coverage_pct = as.integer(sub("cov", "", Coverage)),
    Coverage  = factor(paste0(coverage_pct, "%"), levels = paste0(sort(unique(coverage_pct)), "%"))
  )

pal_mechanism <- c("Disease blocking only" = "#4C72B0",
                    "Disease and infection blocking" = "#C44E52")

theme_seromix <- function(base_size = 20) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey85", colour = NA),
      strip.text = element_text(face = "bold", size = rel(1.0)),
      axis.title = element_text(face = "bold", size = rel(0.95)),
      axis.text = element_text(colour = "black", size = rel(0.85)),
      legend.position = "bottom",
      legend.text = element_text(size = rel(0.9)),
      legend.title = element_text(size = rel(0.95)),
      plot.title = element_text(face = "bold", size = rel(1.1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(0.7))
    )
}

# Stacked design (chat record 2026-09-02): a single 2 (age) x 12 (setting x
# coverage) grid was too wide to read at a legible font size. Stacking one
# block per setting -- each its own AgeCat x Coverage grid (matching the
# original 2x4 layout) -- keeps every panel a normal aspect ratio, so text
# can be sized generously while the overall figure grows taller, not wider.
build_setting_block <- function(setting_name) {
  d <- plot_df %>% filter(setting == setting_name)
  ggplot(d, aes(x = week, y = seroneg_med, colour = Mechanism, group = Mechanism)) +
    geom_ribbon(aes(ymin = seroneg_lo, ymax = seroneg_hi, fill = Mechanism), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.6) +
    facet_grid(AgeCat ~ Coverage) +
    scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
    scale_fill_manual(values = pal_mechanism, name = "Mechanism") +
    scale_x_continuous(name = "Campaign start week", breaks = c(1, 16, 32, 42)) +
    scale_y_continuous(name = "Susceptible at\nvaccination fraction (%)", labels = label_number(suffix = "%"), limits = c(0, 100)) +
    labs(title = setting_name) +
    theme_seromix() +
    theme(plot.title = element_text(face = "bold", size = rel(1.15)))
}

p_low <- build_setting_block("Low")
p_moderate <- build_setting_block("Moderate")
p_high <- build_setting_block("High")

p <- (p_low / p_moderate / p_high) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Modelled susceptible fraction among vaccine recipients by campaign start week",
    theme = theme(
      plot.title = element_text(face = "bold", size = 22)
    )
  ) &
  theme(legend.position = "bottom")

ggsave("06_Results/seromix_seroneg_by_week_mechanism.png", p, width = 13, height = 18, dpi = 300, bg = "white")
message("Saved: 06_Results/seromix_seroneg_by_week_mechanism.png")

write_xlsx(as.data.frame(plot_df %>% arrange(AgeCat, Coverage, Mechanism, week)),
           "06_Results/seromix_seroneg_by_week_mechanism.xlsx")
message("Saved: 06_Results/seromix_seroneg_by_week_mechanism.xlsx")

print(as.data.frame(plot_df %>% arrange(Mechanism, AgeCat, Coverage, week)))
