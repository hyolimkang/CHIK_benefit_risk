# =============================================================================
# 49_infection_surface_plus_threshold_curves.R
#
# Combined travel figure:
#   a = infection-risk surface (travel duration x entry week, by setting)
#       -- p_A from 17_figure2_travel_risk_surface.R
#   b = Pr(BRR > 1) vs cumulative infection risk, DALY | SAE | Death x age
#       -- p_all from 47_infection_risk_threshold_curve_by_outcome.R
# Both sources are run as-is (their own outputs are re-saved unchanged); only
# titles/tags and text sizes are adjusted here so the two blocks match.
#
# Output: 07_Final_Results/Main_Figures/travel_infection_surface_threshold_curves.png
#         07_Final_Results/Main_Figures/travel_infection_surface_threshold_curves_landscape.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(ggplot2)
  library(patchwork)
})

env17 <- new.env(parent = globalenv())
sys.source("02_Scripts/17_figure2_travel_risk_surface.R", envir = env17)
env47 <- new.env(parent = globalenv())
sys.source("02_Scripts/47_infection_risk_threshold_curve_by_outcome.R", envir = env47)

# Match panel a's text to panel b (base_size 17 in script 47).
text_17 <- theme(
  plot.title   = element_text(face = "bold", size = 20),
  strip.text   = element_text(face = "bold", size = 20),
  axis.title   = element_text(face = "bold", size = 17),
  axis.text    = element_text(size = 16),
  legend.title = element_text(face = "bold", size = 16),
  legend.text  = element_text(size = 15)
)

p_a <- env17$p_A + labs(title = NULL) + text_17 +
  theme(legend.key.height = unit(22, "pt"))
p_b <- env47$p_all + labs(title = NULL)

p <- (p_a / p_b) +
  plot_layout(heights = c(1, 2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 26))

out_file <- "07_Final_Results/Main_Figures/travel_infection_surface_threshold_curves.png"
ggsave(out_file, p, width = 16, height = 13.5, dpi = 300, bg = "#fcfcfb")
message("Saved: ", out_file)

# ---- Landscape version (Word landscape page) --------------------------------
# a: the three settings stacked vertically (strip labels on the right, like
#    b's age labels) with a horizontal colour bar underneath;
# b: unchanged 2 x 3 grid to the right.
p_a_land <- env17$p_A + labs(title = NULL) +
  facet_grid(setting ~ .) +
  text_17 +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.key.width = unit(62, "pt"),
        legend.text = element_text(size = 13),
        legend.key.height = unit(12, "pt"),
        axis.text.x = element_text(size = 15))

p_land <- (p_a_land | p_b) +
  plot_layout(widths = c(1, 2.7)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 26))

out_file <- "07_Final_Results/Main_Figures/travel_infection_surface_threshold_curves_landscape.png"
ggsave(out_file, p_land, width = 20, height = 10.5, dpi = 300, bg = "#fcfcfb")
message("Saved: ", out_file)
