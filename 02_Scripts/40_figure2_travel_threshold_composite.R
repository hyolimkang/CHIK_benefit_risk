# =============================================================================
# 40_figure2_travel_threshold_composite.R
#
# Separate composite figure:
#   A. Existing infection-risk surface
#   B. Probability that vaccination benefit exceeds risk across cumulative
#      infection risk, with the 90% pooled threshold
# The DALY surface is intentionally excluded and will be sent separately.
#
# Panel B is the threshold curve from 39_figure2...; the threshold is the
# interpolated crossing of the pooled binned curve, with no extrapolation.
# =============================================================================

project_root <- getwd()
if (!file.exists(file.path(project_root, "01_Data/psa_grid_bra_travel_finite_setting.RData"))) {
  stop("Run this script from the project root: ", project_root)
}

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

surface_env <- new.env(parent = globalenv())
sys.source("02_Scripts/17_figure2_travel_risk_surface.R", envir = surface_env)

threshold_env <- new.env(parent = globalenv())
sys.source("02_Scripts/39_figure2_infection_risk_threshold_curve.R", envir = threshold_env)

p_A <- surface_env$p_A +
  labs(title = "A. Infection-risk surface")

p_B <- threshold_env$p +
  labs(title = "B. Probability that vaccination benefit exceeds risk across cumulative infection risk")

p_fig <- (p_A / p_B) +
  plot_layout(heights = c(1, 1.25))

out_files <- c(
  "06_Results/figure2_travel_risk_surfaces_with_infection_threshold.png",
  "07_Final_Results/Main_Figures/figure2_travel_risk_surfaces_with_infection_threshold.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p_fig, width = 11, height = 11, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}