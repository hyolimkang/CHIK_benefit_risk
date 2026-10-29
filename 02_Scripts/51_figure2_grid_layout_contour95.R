# =============================================================================
# 51_figure2_grid_layout_contour95.R
#
# Figure 2 (42_figure2_grid_layout.R) with the panel B dotted contour drawn at
# Pr(DALY BRR > 1) = 0.95 instead of 0.5, i.e. the entry timing x travel
# duration combinations above which benefit exceeds risk in >=95% of
# posterior draws (rather than the median BRR = 1 line).
#
# Output: 06_Results/figure2_grid_layout_contour95.png
#         07_Final_Results/Main_Figures/figure2_grid_layout_contour95.png
# =============================================================================

contour_prob <- 0.95
source("02_Scripts/42_figure2_grid_layout.R")
