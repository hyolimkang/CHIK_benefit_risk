# =============================================================================
# 46_daly_weeksweep_probgt1_panelE.R
#
# Supplementary Pr(BRR>1) week x coverage heatmaps (outbreak response), one
# figure per outcome, relabelled a/b/c (previously C/D/E):
#   a = SAE, b = Death, c = DALY
# Each has an outcome title and the two mechanism rows, same style as the
# SAE/Death branch of make_prob_heatmap_combined() in
# finite_weeksweep_heatmap_probgt1.R. The title-less "B"-tagged DALY version
# used as a main figure (07_Final_Results/Main_Figures) is left untouched.
#
# Reuses finite_weeksweep_heatmap_probgt1.R (CHIK_ORV_impact copy, the one
# actually executed) by evaluating every top-level expression EXCEPT the final
# pwalk() that saves the combined B/C/D figures, so no existing output there
# is overwritten. `plots` (per-mechanism heatmaps) comes from that script.
#
# Output: 07_Final_Results/Supplementary_Figures/
#           brr_weeksweep_heatmap_probgt1_{sae,death,daly}_{base,adj}_mechanism_combined.png
# =============================================================================

br_root  <- "c:/Users/user/OneDrive/CHIK_benefit_risk"
ori_root <- "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"

setwd(ori_root)
src_exprs <- parse("01_Script/1_4_Analysis/finite_weeksweep_heatmap_probgt1.R")
for (e in src_exprs[-length(src_exprs)]) eval(e, envir = globalenv())

PANEL_LETTERS <- c(SAE = "a", Death = "b", DALY = "c")

make_supp_panel <- function(outcome, version) {
  p_ve0   <- plots[[paste(outcome, version, "ve0",   sep = "_")]] +
    labs(title = "Disease blocking only", caption = NULL)
  p_ve989 <- plots[[paste(outcome, version, "ve989", sep = "_")]] +
    labs(title = "Disease and infection blocking", caption = NULL)

  # Same title/margin settings as the SAE/Death branch of
  # make_prob_heatmap_combined().
  p_combined <- (p_ve0 / p_ve989) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = outcome,
      theme = theme(
        plot.title  = element_text(face = "bold", size = 17, colour = ink_primary,
                                    margin = margin(t = 12)),
        plot.margin = margin(t = 34, r = 5.5, b = 5.5, l = 5.5)
      )
    ) &
    theme(legend.position = "right")

  p_tagged <- patchwork::wrap_elements(full = p_combined) +
    plot_annotation(tag_levels = list(PANEL_LETTERS[[outcome]])) &
    theme(plot.tag = element_text(face = "bold", size = 21, colour = ink_primary))

  out_file <- file.path(br_root, sprintf(
    "07_Final_Results/Supplementary_Figures/brr_weeksweep_heatmap_probgt1_%s_%s_mechanism_combined.png",
    tolower(outcome), version))
  ggsave(out_file, p_tagged, width = 9.6, height = 9.6, dpi = 600, bg = "white")
  message("Saved: ", out_file)
}

for (oc in names(PANEL_LETTERS)) for (v in names(PROB_VERSIONS)) make_supp_panel(oc, v)
