# =============================================================================
# 38_headline_brplane_sae_death_supplementary.R
#
# Supplementary companions to the headline benefit-risk plane (Panel A of the
# main manuscript figure, DALY only -- 36_headline_brplane_ab_tags.R), same
# style/layout, for SAE and Death.
#
# Reuses build_setting_synthesis_data()/build_benefit_risk_plane() verbatim
# from 10_streamline_mechanism_figure.R (source of truth) / 36's copy -- those
# already exclude 18-64 for Death (`!(outcome == "Death" & AgeCat == "18-64")`),
# since vaccine-attributable death risk is fixed at exactly 0 for that age
# group (01_setup.R's p_death_vacc_u65 <- 0), so brr_death is undefined there.
#
# No outer "A"/"B" corner tag here (that pairing is specific to the DALY
# figure's role as Panel A of the main manuscript figure); otherwise same
# style: no overall title, "Base risk"/"Serostatus-adjusted risk" row tags,
# same enlarged font sizes as the DALY version.
#
# Output: 07_Final_Results/Supplementary_Figures/headline_benefit_risk_plane_SAE.png
#         07_Final_Results/Supplementary_Figures/headline_benefit_risk_plane_Death.png
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

build_setting_synthesis_data <- function(outcome_sel, use_adj = FALSE) {
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

build_benefit_risk_plane <- function(outcome_sel, use_adj = FALSE, shared_rng = NULL) {
  d <- build_setting_synthesis_data(outcome_sel, use_adj)

  rng <- shared_rng %||% pmax(range(c(d$risk_med, d$averted_med, d$risk_lo, d$risk_hi,
                                       d$averted_lo, d$averted_hi), na.rm = TRUE), 1e-3)

  risk_label <- if (use_adj) paste0(outcome_sel, " caused by vaccination (per 10,000 vaccinated individuals)\n(serostatus-adjusted, log scale)") else
    paste0(outcome_sel, " caused by vaccination (per 10,000 vaccinated individuals)\n(base, log scale)")

  age_levels_present <- levels(droplevels(d$AgeCat %>% factor(levels = c("18-64", "65+"))))
  age_cols <- c("18-64" = "#4C72B0", "65+" = "#DD8452")[age_levels_present]

  ggplot(d, aes(risk_med, averted_med, colour = AgeCat, shape = setting)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(ymin = averted_lo, ymax = averted_hi), width = 0, alpha = 0.3) +
    geom_errorbarh(aes(xmin = risk_lo, xmax = risk_hi), height = 0, alpha = 0.3) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~VE_label) +
    coord_equal(xlim = rng, ylim = rng) +
    scale_colour_manual(values = age_cols, name = "Age group", drop = TRUE) +
    scale_shape_manual(values = c(Low = 15, Moderate = 17, High = 16), name = "Setting") +
    scale_x_log10(labels = log_num_labels, name = risk_label) +
    scale_y_log10(labels = log_num_labels, name = paste(outcome_sel, "averted per 10,000 vaccinated (log scale)")) +
    theme_nm()
}

build_and_save_brplane_supp <- function(outcome_sel, panel_label) {
  d_base_rng <- build_setting_synthesis_data(outcome_sel, FALSE)
  d_adj_rng  <- build_setting_synthesis_data(outcome_sel, TRUE)
  shared_rng <- pmax(range(c(d_base_rng$risk_med, d_base_rng$averted_med, d_base_rng$risk_lo, d_base_rng$risk_hi,
                              d_base_rng$averted_lo, d_base_rng$averted_hi,
                              d_adj_rng$risk_med, d_adj_rng$risk_lo, d_adj_rng$risk_hi), na.rm = TRUE), 1e-3)

  p_base <- build_benefit_risk_plane(outcome_sel, use_adj = FALSE, shared_rng = shared_rng)
  p_adj  <- build_benefit_risk_plane(outcome_sel, use_adj = TRUE,  shared_rng = shared_rng)

  p_combined <- (p_base / p_adj) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = panel_label,
      tag_levels = list(c("Base risk", "Serostatus-adjusted risk")),
      # Extra top margin opens blank space above the title -- the "Base risk"
      # row tag (below) sits there instead of colliding with it (both anchor
      # near the same top-left corner by default).
      theme = theme(plot.title = element_text(face = "bold", size = 16),
                    plot.tag   = element_text(face = "bold", size = 15),
                    plot.margin = margin(t = 22, r = 5.5, b = 5.5, l = 5.5))
    ) &
    theme(legend.position = "bottom")

  dir.create("07_Final_Results/Supplementary_Figures", showWarnings = FALSE, recursive = TRUE)
  out_file <- sprintf("07_Final_Results/Supplementary_Figures/headline_benefit_risk_plane_%s.png", outcome_sel)
  ggsave(out_file, p_combined, width = 10, height = 10, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p_combined)
}

build_and_save_brplane_supp("SAE",   "A. SAE benefit-risk space")
build_and_save_brplane_supp("Death", "B. Death benefit-risk space")
