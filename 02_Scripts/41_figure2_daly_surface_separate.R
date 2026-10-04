# =============================================================================
# 41_figure2_daly_surface_separate.R
#
# Standalone C panel from Figure 2: DALY benefit-risk surface.
# =============================================================================

project_root <- getwd()
if (!file.exists(file.path(project_root, "01_Data/psa_grid_bra_travel_finite_setting.RData"))) {
  stop("Run this script from the project root: ", project_root)
}

suppressMessages({
  library(ggplot2)
})

surface_env <- new.env(parent = globalenv())
sys.source("02_Scripts/17_figure2_travel_risk_surface.R", envir = surface_env)

p_C <- surface_env$p_B +
  labs(title = "C. DALY benefit-risk surface")

out_files <- c(
  "06_Results/figure2_daly_benefit_risk_surface.png",
  "07_Final_Results/Main_Figures/figure2_daly_benefit_risk_surface.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p_C, width = 11, height = 8.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}
