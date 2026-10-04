# =============================================================================
# 42_figure2_grid_layout.R
#
# Alternative Figure 2 composite (chat record 2026-09-05): the previous
# top/bottom composite (40_figure2_travel_threshold_composite.R: A on top,
# threshold curve below, both full-width) squeezed every facet too small.
# This lays out the SAME three plots as a 2x2-ish grid instead:
#
#   A   B1
#   C   B2
#
# Left column: A (infection-risk surface, 1x3 facets) stacked above
# C (DALY benefit-risk surface, 2x3 facets, age x setting) -- these are
# 17_figure2_travel_risk_surface.R's own p_A/p_B, unchanged.
# Right column: the threshold curve (39_figure2_infection_risk_threshold_
# curve.R's p), with its age_group facet switched from side-by-side
# (nrow=1) to stacked (ncol=1) so its two panels (B1=18-64, B2=>=65) occupy
# the same column instead of being squeezed next to each other.
#
# Output: 06_Results/figure2_grid_layout.png
#         07_Final_Results/Main_Figures/figure2_grid_layout.png
# =============================================================================

project_root <- getwd()
if (!file.exists(file.path(project_root, "01_Data/psa_grid_bra_travel_finite_setting.RData"))) {
  stop("Run this script from the project root: ", project_root)
}

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

env17 <- new.env(parent = globalenv())
sys.source("02_Scripts/17_figure2_travel_risk_surface.R", envir = env17)
# -> env17$p_A ("A. Infection-risk surface"), env17$p_B ("B. DALY benefit-risk surface")

env39 <- new.env(parent = globalenv())
sys.source("02_Scripts/39_figure2_infection_risk_threshold_curve.R", envir = env39)
# -> env39$p ("B. Probability that DALY benefit exceeds risk...", faceted nrow=1)

p_A <- env17$p_A
p_C <- env17$p_B + labs(title = "C. DALY benefit-risk surface")

p_B_stacked <- env39$p +
  facet_wrap(~age_group, ncol = 1) +
  labs(title = "B. Probability that DALY benefit exceeds risk\nacross cumulative infection risk")

left_col  <- (p_A / p_C) + plot_layout(heights = c(1, 2))
# Give the probability curves a little more horizontal room than in the
# original 2:1 split, while retaining the surface panels as the larger column.
p_grid    <- (left_col | p_B_stacked) + plot_layout(widths = c(1.8, 1.2))

# Wider canvas + larger text throughout (chat record 2026-09-05) -- the
# previous version's font sizes (inherited from theme_hm()/theme_threshold(),
# tuned for each panel standalone) read small once assembled into this wider
# grid. `&` applies the override to every sub-plot in the patchwork at once.
p_grid <- p_grid & theme(
  plot.title    = element_text(size = 19, face = "bold"),
  plot.subtitle = element_text(size = 13),
  strip.text    = element_text(size = 16, face = "bold"),
  axis.title    = element_text(size = 16, face = "bold"),
  axis.text     = element_text(size = 14),
  legend.title  = element_text(size = 15, face = "bold"),
  legend.text   = element_text(size = 13)
)

out_files <- c(
  "06_Results/figure2_grid_layout.png",
  "07_Final_Results/Main_Figures/figure2_grid_layout.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p_grid, width = 22, height = 14, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}
