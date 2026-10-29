# =============================================================================
# 48_headline_brplane_three_outcomes.R
#
# DALY | SAE | Death benefit-risk planes side by side (one column per
# outcome, rows = protection mechanism), one figure per risk assumption.
#
# Reuses build_setting_synthesis_data()/build_benefit_risk_plane() from
# 38_headline_brplane_sae_death_supplementary.R (same data, same setting
# x-offset jitter) by evaluating that script's definitions without its
# final save calls. Each outcome keeps its own log axis range (shared between
# base and serostatus-adjusted) because the three outcomes differ by orders
# of magnitude; coord_equal keeps the 1:1 diagonal meaningful in each.
#
# Death excludes 18-64 (vaccine-attributable death risk fixed at 0, so BRR
# is undefined) -- same rule as scripts 36/38.
#
# Output: 07_Final_Results/Supplementary_Figures/
#           headline_benefit_risk_plane_DALY_SAE_Death_base.png
#           headline_benefit_risk_plane_DALY_SAE_Death_adj.png
#         07_Final_Results/Main_Figures/headline_benefit_risk_plane_DALY_SAE_Death_AB.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
src <- parse("02_Scripts/38_headline_brplane_sae_death_supplementary.R")
is_save_call <- vapply(src, function(e) {
  is.call(e) && identical(e[[1]], as.name("build_and_save_brplane_supp"))
}, logical(1))
for (e in src[!is_save_call]) eval(e, envir = globalenv())

OUTCOMES <- c("DALY", "SAE", "Death")

shared_range <- function(outcome_sel) {
  b <- build_setting_synthesis_data(outcome_sel, FALSE)
  a <- build_setting_synthesis_data(outcome_sel, TRUE)
  pmax(range(c(b$risk_lo, b$risk_hi, b$averted_lo, b$averted_hi,
               a$risk_lo, a$risk_hi, a$averted_lo, a$averted_hi), na.rm = TRUE), 1e-3)
}

outcome_panel <- function(outcome_sel, use_adj) {
  build_benefit_risk_plane(outcome_sel, use_adj = use_adj, shared_rng = shared_range(outcome_sel)) +
    facet_wrap(~VE_label, ncol = 1) +
    # Both colour levels always declared so the collected legend shows both
    # age groups even though the Death panel only contains 65+.
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"),
                        limits = c("18-64", "65+"), name = "Age group") +
    # Axis names live in the scales built by build_benefit_risk_plane(), so
    # replace the scales (labs() alone would not override them).
    scale_x_log10(labels = log_num_labels,
                  name = "Caused by vaccination\n(per 10,000 vaccinated, log scale)") +
    scale_y_log10(labels = log_num_labels,
                  name = "Averted by vaccination\n(per 10,000 vaccinated, log scale)") +
    labs(title = outcome_sel) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
}

for (use_adj in c(FALSE, TRUE)) {
  panels <- lapply(OUTCOMES, outcome_panel, use_adj = use_adj)
  # Y-axis title only on the first column to save horizontal space.
  panels[2:3] <- lapply(panels[2:3], function(p) p + theme(axis.title.y = element_blank()))
  # Death has no 18-64 data, so its colour key renders differently and is not
  # merged by guides = "collect"; the DALY/SAE legend already covers it.
  panels[[3]] <- panels[[3]] + guides(colour = "none")

  risk_title <- if (use_adj) "Serostatus-adjusted risk" else "Base risk"
  p <- wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = risk_title,
      caption = "Death: adults aged 18-64 years not shown (vaccine-attributable death risk assumed to be zero). Dashed line: benefit = risk.",
      theme = theme(plot.title = element_text(face = "bold", size = 18),
                    plot.caption = element_text(colour = "grey35", size = 10, hjust = 0))
    ) &
    theme(legend.position = "bottom")

  out_file <- sprintf("07_Final_Results/Supplementary_Figures/headline_benefit_risk_plane_DALY_SAE_Death_%s.png",
                      if (use_adj) "adj" else "base")
  ggsave(out_file, p, width = 14, height = 9.5, dpi = 300, bg = "white")
  message("Saved: ", out_file)
}

# ---- Combined: A = base risk, B = serostatus-adjusted risk (MAIN figure) ----
# Each risk-assumption row is rendered as its own block (wrap_elements) with
# the panel letter as title and the risk assumption as subtitle; the shared
# legend + caption are shown once, under row B. Text is enlarged for the
# main-figure version so it stays legible once scaled to a Word page.
main_text <- theme(
  plot.title  = element_text(face = "bold", size = 20, hjust = 0.5),
  strip.text  = element_text(face = "bold", size = 16),
  axis.title  = element_text(size = 16),
  axis.text   = element_text(size = 14),
  legend.title = element_text(size = 16),
  legend.text  = element_text(size = 15)
)

make_risk_row <- function(use_adj, letter, show_legend) {
  panels <- lapply(OUTCOMES, outcome_panel, use_adj = use_adj)
  # Larger strip text no longer fits "Disease and infection blocking" on one
  # line inside the square (coord_equal) panels -> break the label onto two
  # lines at the data level.
  panels <- lapply(panels, function(p) {
    p <- p + main_text
    lv <- levels(p$data$VE_label)
    p$data$VE_label <- factor(p$data$VE_label, levels = lv,
                              labels = sub("Disease and infection blocking",
                                           "Disease and infection\nblocking", lv, fixed = TRUE))
    p
  })
  panels[2:3] <- lapply(panels[2:3], function(p) p + theme(axis.title.y = element_blank()))
  panels[[3]] <- panels[[3]] + guides(colour = "none")
  row <- wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = letter,
      subtitle = if (use_adj) "Serostatus-adjusted risk" else "Base risk",
      # no in-figure caption: Death 18-64 note and diagonal meaning are in the figure legend text
      theme = theme(plot.title = element_text(face = "bold", size = 28),
                    plot.subtitle = element_text(face = "bold", size = 21),
                    plot.caption = element_text(colour = "grey35", size = 14, hjust = 0))
    ) &
    theme(legend.position = if (show_legend) "bottom" else "none")
  wrap_elements(full = row)
}

p_ab <- make_risk_row(FALSE, "A", show_legend = FALSE) /
        make_risk_row(TRUE,  "B", show_legend = TRUE) +
  plot_layout(heights = c(1, 1.1))

dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
out_file <- "07_Final_Results/Main_Figures/headline_benefit_risk_plane_DALY_SAE_Death_AB.png"
ggsave(out_file, p_ab, width = 15, height = 18, dpi = 300, bg = "white")
message("Saved: ", out_file)
