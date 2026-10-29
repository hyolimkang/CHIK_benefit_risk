# =============================================================================
# 50_headline_brplane_base_adj_overlay.R
#
# Single-row alternative to Figure_2A_new.png (script 36): base and
# serostatus-adjusted risk overlaid in the same DALY benefit-risk plane
# instead of two near-identical rows. Benefit (DALY averted) is the same
# under both assumptions, so each base/adjusted pair shares a y value; an
# arrow marks the x shift (adjusted = base x seronegative share among
# vaccinees, RR_seropos = 0).
#
# Also writes a linear-axis version (rows = age group, free scales) to check
# whether the shift reads more clearly without the log scale.
#
# Reuses build_setting_synthesis_data()/theme_nm() from 36 by evaluating that
# script's definitions without its final save call.
#
# Output: 06_Results/Figure_2A_base_adj_overlay.png
#         06_Results/Figure_2A_base_adj_overlay_linear.png
#         07_Final_Results/Supplementary_Figures/ (same two files)
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
src <- parse("02_Scripts/36_headline_brplane_ab_tags.R")
is_save_call <- vapply(src, function(e) {
  is.call(e) && identical(e[[1]], as.name("build_and_save_brplane"))
}, logical(1))
for (e in src[!is_save_call]) eval(e, envir = globalenv())

outcome_sel <- "DALY"

d <- bind_rows(
  build_setting_synthesis_data(outcome_sel, FALSE) %>% mutate(risk_type = "Base"),
  build_setting_synthesis_data(outcome_sel, TRUE)  %>% mutate(risk_type = "Serostatus-adjusted")
) %>%
  mutate(
    risk_type = factor(risk_type, levels = c("Base", "Serostatus-adjusted")),
    # Same per-setting log10 x offset as script 36, so settings don't overlap.
    x_offset = 10^case_when(setting == "Low" ~ -0.06, setting == "High" ~ 0.06, TRUE ~ 0),
    risk_med_plot = risk_med * x_offset,
    risk_lo_plot  = risk_lo  * x_offset,
    risk_hi_plot  = risk_hi  * x_offset
  )

d_shift <- d %>%
  select(setting, VE_label, AgeCat, risk_type, risk_med_plot, averted_med) %>%
  tidyr::pivot_wider(names_from = risk_type, values_from = risk_med_plot)

rng <- pmax(range(c(d$risk_lo, d$risk_hi, d$averted_lo, d$averted_hi), na.rm = TRUE), 1e-3)

age_cols <- c("18-64" = "#4C72B0", "65+" = "#DD8452")
# Fillable shapes: base = filled in age colour, adjusted = open (white fill).
setting_shapes <- c(Low = 22, Moderate = 24, High = 21)

plane_layers <- function(d, d_shift) {
  list(
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40"),
    geom_errorbar(data = d %>% filter(risk_type == "Base"),
                  aes(ymin = averted_lo, ymax = averted_hi),
                  width = 0, alpha = 0.3),
    geom_errorbarh(aes(xmin = risk_lo_plot, xmax = risk_hi_plot, linetype = risk_type),
                   height = 0, alpha = 0.3),
    geom_segment(data = d_shift,
                 aes(x = Base, xend = `Serostatus-adjusted`, y = averted_med, yend = averted_med,
                     colour = AgeCat),
                 inherit.aes = FALSE, linewidth = 0.4,
                 arrow = arrow(length = unit(4, "pt"), type = "closed")),
    geom_point(aes(fill = interaction(AgeCat, risk_type)), size = 3, stroke = 0.9),
    # Invisible layer: geom_point has no linetype aesthetic, so this one only
    # puts a filled/open point into the Risk assumption legend key.
    geom_pointrange(aes(ymin = averted_med, ymax = averted_med, linetype = risk_type),
                    alpha = 0, key_glyph = "point"),
    scale_colour_manual(values = age_cols, name = "Age group"),
    scale_fill_manual(values = c("18-64.Base" = age_cols[["18-64"]], "65+.Base" = age_cols[["65+"]],
                                 "18-64.Serostatus-adjusted" = "white",
                                 "65+.Serostatus-adjusted" = "white"),
                      guide = "none"),
    scale_shape_manual(values = setting_shapes, name = "Setting"),
    scale_linetype_manual(values = c("Base" = "solid", "Serostatus-adjusted" = "22"),
                          name = "Risk assumption"),
    guides(shape = guide_legend(order = 1, override.aes = list(fill = "grey40")),
           colour = guide_legend(order = 2),
           linetype = guide_legend(order = 3, override.aes = list(alpha = 0.8, colour = "grey30", shape = 21, fill = c("grey30", "white"),
                                                                 size = 3, stroke = 0.9))),
    theme_nm(),
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.spacing.y = unit(2, "pt"), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(4, "pt"))
  )
}

p_log <- ggplot(d, aes(risk_med_plot, averted_med, colour = AgeCat, shape = setting)) +
  plane_layers(d, d_shift) +
  facet_wrap(~VE_label) +
  coord_equal(xlim = rng, ylim = rng) +
  scale_x_log10(labels = log_num_labels,
                name = paste(outcome_sel, "caused by vaccination (per 10,000 vaccinated individuals, log scale)")) +
  scale_y_log10(labels = log_num_labels,
                name = paste(outcome_sel, "averted per 10,000 vaccinated (log scale)"))

# Linear version: no setting offset (it was a log-scale device) and free
# scales per panel, since 18-64 risk (~0.15) would otherwise sit at zero
# next to 65+ (~5).
d_lin <- d %>% mutate(risk_med_plot = risk_med, risk_lo_plot = risk_lo, risk_hi_plot = risk_hi)
d_shift_lin <- d_lin %>%
  select(setting, VE_label, AgeCat, risk_type, risk_med_plot, averted_med) %>%
  tidyr::pivot_wider(names_from = risk_type, values_from = risk_med_plot)

p_lin <- ggplot(d_lin, aes(risk_med_plot, averted_med, colour = AgeCat, shape = setting)) +
  plane_layers(d_lin, d_shift_lin) +
  facet_wrap(AgeCat ~ VE_label, ncol = 2, scales = "free") +
  scale_x_continuous(name = paste(outcome_sel, "caused by vaccination (per 10,000 vaccinated individuals, linear scale)"),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(name = paste(outcome_sel, "averted per 10,000 vaccinated (linear scale)"),
                     limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))

dir.create("07_Final_Results/Supplementary_Figures", showWarnings = FALSE, recursive = TRUE)
for (dir in c("06_Results", "07_Final_Results/Supplementary_Figures")) {
  out_log <- file.path(dir, "Figure_2A_base_adj_overlay.png")
  # Sized for a landscape Word page (A4, 2.54 cm margins: ~24.7 x 16 cm).
  ggsave(out_log, p_log, width = 9.5, height = 6.2, dpi = 300, bg = "white")
  message("Saved: ", out_log)
  out_lin <- file.path(dir, "Figure_2A_base_adj_overlay_linear.png")
  ggsave(out_lin, p_lin, width = 10, height = 9, dpi = 300, bg = "white")
  message("Saved: ", out_lin)
}
