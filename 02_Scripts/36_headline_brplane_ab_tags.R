# =============================================================================
# 36_headline_brplane_ab_tags.R
#
# Lightweight rebuild of headline_benefit_risk_plane_DALY.png with "A."/"B."
# panel tags (chat record 2026-09-02) instead of bare "Base risk"/
# "Serostatus-adjusted risk" row labels, matching the panel-letter convention
# used elsewhere (figure2_travel_risk_surface.R's "A. Infection-risk surface",
# 23_travel_brr_range_by_timing.R's "C. DALY BRR range by entry timing").
#
# Reuses build_setting_synthesis_data()/build_benefit_risk_plane()/
# build_and_save_brplane() verbatim from 10_streamline_mechanism_figure.R
# (source of truth for that logic, updated in place with the same tag change)
# -- but those functions only need draw_level_xy_serostatus_finite.RData, not
# the 1.3GB posterior_finite_all.RData that script also loads for its other
# (unrelated) figures, so this rebuilds just the one figure without that cost.
#
# Output (preserves the existing Figure_2A.png):
#         06_Results/Figure_2A_new.png
#         07_Final_Results/Main_Figures/Figure_2A_new.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

load("01_Data/draw_level_xy_serostatus_finite.RData")

log_num_labels <- scales::label_number(accuracy = NULL, big.mark = ",")

theme_nm <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = "") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.line = element_line(linewidth = 0.3, colour = "grey30"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
      strip.text = element_text(face = "bold", size = rel(0.95)),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.9)),
      legend.text = element_text(size = rel(0.85))
    )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

build_setting_synthesis_data <- function(outcome_sel = "DALY", use_adj = FALSE) {
  brr_col  <- if (use_adj) "brr_adj" else "brr_base"
  risk_col <- if (use_adj) "x_10k_adj" else "x_10k_base"

  draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel,
           AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64")) %>%
    group_by(setting, VE_label, AgeCat) %>%
    summarise(
      averted_med = median(averted_10k, na.rm = TRUE),
      averted_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
      averted_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),
      risk_med = median(.data[[risk_col]], na.rm = TRUE),
      risk_lo  = quantile(.data[[risk_col]], 0.025, na.rm = TRUE),
      risk_hi  = quantile(.data[[risk_col]], 0.975, na.rm = TRUE),
      brr_med = median(.data[[brr_col]], na.rm = TRUE),
      brr_lo  = quantile(.data[[brr_col]], 0.025, na.rm = TRUE),
      brr_hi  = quantile(.data[[brr_col]], 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(setting = factor(setting, levels = c("Low", "Moderate", "High"))) %>%
    mutate(across(c(averted_med, averted_lo, averted_hi,
                     risk_med, risk_lo, risk_hi,
                     brr_med, brr_lo, brr_hi),
                   ~ pmax(.x, 1e-4)))
}

build_benefit_risk_plane <- function(outcome_sel = "DALY", use_adj = FALSE,
                                      shared_rng = NULL) {
  d <- build_setting_synthesis_data(outcome_sel, use_adj) %>%
    # Risk is identical across transmission settings within each age group,
    # so the three point/error-bar glyphs otherwise sit on the same vertical
    # line. Apply a small, deterministic offset on the log10 x scale purely
    # for visibility; all x components are shifted together so each point
    # remains centred on its horizontal uncertainty interval.
    mutate(
      x_offset = 10^case_when(
        setting == "Low"      ~ -0.06,
        setting == "Moderate" ~  0,
        setting == "High"     ~  0.06,
        TRUE                  ~  0
      ),
      risk_med_plot = risk_med * x_offset,
      risk_lo_plot  = risk_lo  * x_offset,
      risk_hi_plot  = risk_hi  * x_offset
    )

  rng <- shared_rng %||% pmax(range(c(d$risk_med, d$averted_med, d$risk_lo, d$risk_hi,
                                       d$averted_lo, d$averted_hi), na.rm = TRUE), 1e-3)

  risk_label <- if (use_adj) paste0(outcome_sel, " caused by vaccination (per 10,000 vaccinated individuals)\n(serostatus-adjusted, log scale)") else
    paste0(outcome_sel, " caused by vaccination (per 10,000 vaccinated individuals)\n(base, log scale)")

  ggplot(d, aes(risk_med_plot, averted_med, colour = AgeCat, shape = setting)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(
      aes(ymin = averted_lo, ymax = averted_hi,
          linetype = "95% UI"),
      width = 0, alpha = 0.3,
      show.legend = c(colour = FALSE, linetype = TRUE)
    ) +
    geom_errorbarh(
      aes(xmin = risk_lo_plot, xmax = risk_hi_plot,
          linetype = "95% UI"),
      height = 0, alpha = 0.3,
      show.legend = c(colour = FALSE, linetype = TRUE)
    ) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~VE_label) +
    coord_equal(xlim = rng, ylim = rng) +
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"), name = "Age group") +
    scale_shape_manual(values = c(Low = 15, Moderate = 17, High = 16), name = "Setting") +
    scale_linetype_manual(
      values = c("95% UI" = "solid"),
      name = "Bars"
    ) +
    scale_x_log10(labels = log_num_labels, name = risk_label) +
    scale_y_log10(labels = log_num_labels, name = paste(outcome_sel, "averted per 10,000 vaccinated (log scale)")) +
    guides(
      shape = guide_legend(order = 1),
      colour = guide_legend(order = 2),
      linetype = guide_legend(
        order = 3,
        override.aes = list(colour = "grey40", alpha = 0.6)
      )
    ) +
    theme_nm()
}

build_and_save_brplane <- function(outcome_sel) {
  d_base_rng <- build_setting_synthesis_data(outcome_sel, FALSE)
  d_adj_rng  <- build_setting_synthesis_data(outcome_sel, TRUE)
  shared_rng <- pmax(range(c(d_base_rng$risk_med, d_base_rng$averted_med, d_base_rng$risk_lo, d_base_rng$risk_hi,
                              d_base_rng$averted_lo, d_base_rng$averted_hi,
                              d_adj_rng$risk_med, d_adj_rng$risk_lo, d_adj_rng$risk_hi), na.rm = TRUE), 1e-3)

  # Compact row titles avoid the wide tag column created by the longer
  # "Serostatus-adjusted risk" patchwork tag, leaving more room for panels.
  p_base <- build_benefit_risk_plane(outcome_sel, use_adj = FALSE, shared_rng = shared_rng) +
    labs(title = "A. Base risk") +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
          plot.title.position = "plot")
  p_adj  <- build_benefit_risk_plane(outcome_sel, use_adj = TRUE,  shared_rng = shared_rng) +
    labs(title = "Serostatus-adjusted risk") +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0),
          plot.title.position = "plot")

  p_combined <- (p_base / p_adj) +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(plot.margin = margin(t = 7, r = 5.5, b = 5.5, l = 5.5))
    ) &
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.spacing.x = unit(3, "pt"),
      legend.key.spacing.x = unit(2, "pt"),
      legend.key.width = unit(12, "pt")
    )

  # The figure letter is integrated into the first row title so no outer tag
  # column or wrapper consumes space that should be used by the data panels.
  p_tagged <- p_combined

  out_files <- c(
    "06_Results/Figure_2A_new.png",
    "07_Final_Results/Main_Figures/Figure_2A_new.png"
  )
  dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
  for (out_file in out_files) {
    ggsave(out_file, p_tagged, width = 10, height = 10, dpi = 300, bg = "white")
    message("Saved: ", out_file)
  }
  invisible(p_tagged)
}

build_and_save_brplane("DALY")
