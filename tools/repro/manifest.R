# =============================================================================
# tools/repro/manifest.R
#
# Publication manifest: every output in the Nature Medicine submission that
# can be regenerated from archived model outputs (Mode B), with
#   script : the unmodified scientific script that produces it
#   inputs : files it reads (area "data" = CHIK_DATA_DIR, "orv" = CHIK_ORV_ROOT)
#   checks : reproduced file (relative to CHIK_OUTPUT_DIR) vs reference file
#            (relative to the project root), and how to compare them
#
# References are the files in the submission folder
# (03_Manuscript/Submission/submission_naturemed/), which the corresponding
# author identified as the final versions. They are used ONLY as comparison
# targets, never as inputs. See docs/FIGURE_TABLE_REPRODUCTION_MAP.md.
# =============================================================================

SUB <- "03_Manuscript/Submission/submission_naturemed"

inp <- function(area, file, role) list(area = area, file = file, role = role)
chk <- function(reproduced, reference, type) list(reproduced = reproduced, reference = reference, type = type)

PUBLICATION_MANIFEST <- list(
  list(
    id = "fig1", output = "Figure 1 (travel infection-risk surface, DALY Pr(BRR>1) surface with 95% contour, threshold curves)",
    script = "02_Scripts/51_figure2_grid_layout_contour95.R",
    inputs = list(inp("data", "psa_grid_bra_travel_finite_setting.RData", "setting-level travel grid (entry day x duration)")),
    checks = list(chk("07_Final_Results/Main_Figures/figure2_grid_layout_contour95.png",
                      file.path(SUB, "Main_Figures/Figure_1_Final.png"), "png"))
  ),
  list(
    id = "fig2", output = "Figure 2 (ORI DALY benefit-risk plane, base vs serostatus-adjusted risk)",
    script = "02_Scripts/50_headline_brplane_base_adj_overlay.R",
    inputs = list(inp("data", "draw_level_xy_serostatus_finite.RData", "ORI draw-level benefit/risk by serostatus scenario")),
    checks = list(chk("07_Final_Results/Supplementary_Figures/Figure_2A_base_adj_overlay.png",
                      file.path(SUB, "Main_Figures/Figure_2_Final.png"), "png"))
  ),
  list(
    id = "tab1", output = "Table 1 (travel DALY benefit, risk, BRR by duration)",
    script = "02_Scripts/20_main_table2_travel_daly.R",
    inputs = list(inp("data", "psa_df_bra_travel_finite.RData", "travel draw-level PSA output"),
                  inp("data", "setting_key.RData", "state -> transmission setting")),
    checks = list(chk("06_Results/Main_Table_2_Travel_DALY.xlsx",
                      file.path(SUB, "Main_Tables/Main_Table_2_Travel_DALY.xlsx"), "table"))
  ),
  list(
    id = "tab2", output = "Table 2 (ORI DALY benefit, risk, BRR, baseline scenario)",
    script = "02_Scripts/19_main_table1_daly.R",
    inputs = list(inp("data", "draw_level_xy_serostatus_finite.RData", "ORI draw-level benefit/risk")),
    checks = list(chk("06_Results/Main_Table_1_DALY.xlsx",
                      file.path(SUB, "Main_Tables/Main_Table_1_DALY.xlsx"), "table"))
  ),
  list(
    id = "stab1", output = "Supplementary Table 1 (reporting rate and attack rate by state)",
    script = "02_Scripts/33_state_reporting_rate_table.R",
    inputs = list(inp("data", "setting_key.RData", "state -> transmission setting"),
                  inp("data", "posterior_finite_all.RData", "transmission-model posterior draws (all states)"),
                  inp("data", "state_ar_hazard_summary_finite.RData", "state attack-rate summary (cache; rebuilt from the posterior with --recompute-caches)")),
    checks = list(chk("06_Results/Table_S_state_reporting_rate.xlsx",
                      file.path(SUB, "Main_Tables/Table_S_state_reporting_rate.xlsx"), "table"))
  ),
  list(
    id = "stab2", output = "Supplementary Table 2 (travel SAE)",
    script = "02_Scripts/26_main_table2_travel_sae.R",
    inputs = list(inp("data", "psa_df_bra_travel_finite.RData", "travel draw-level PSA output"),
                  inp("data", "setting_key.RData", "state -> transmission setting")),
    checks = list(chk("06_Results/Main_Table_2_Travel_SAE.xlsx",
                      file.path(SUB, "Main_Tables/Main_Table_2_Travel_SAE.xlsx"), "table"))
  ),
  list(
    id = "stab3", output = "Supplementary Table 3 (travel deaths)",
    script = "02_Scripts/27_main_table2_travel_death.R",
    inputs = list(inp("data", "psa_df_bra_travel_finite.RData", "travel draw-level PSA output"),
                  inp("data", "setting_key.RData", "state -> transmission setting")),
    checks = list(chk("06_Results/Main_Table_2_Travel_Death.xlsx",
                      file.path(SUB, "Main_Tables/Main_Table_2_Travel_Death.xlsx"), "table"))
  ),
  list(
    id = "stab4", output = "Supplementary Table 4 (ORI SAE)",
    script = "02_Scripts/24_main_table1_sae.R",
    inputs = list(inp("data", "draw_level_xy_serostatus_finite.RData", "ORI draw-level benefit/risk")),
    checks = list(chk("06_Results/Main_Table_1_SAE.xlsx",
                      file.path(SUB, "Main_Tables/Main_Table_1_SAE.xlsx"), "table"))
  ),
  list(
    id = "stab5", output = "Supplementary Table 5 (ORI deaths)",
    script = "02_Scripts/25_main_table1_death.R",
    inputs = list(inp("data", "draw_level_xy_serostatus_finite.RData", "ORI draw-level benefit/risk")),
    checks = list(chk("06_Results/Main_Table_1_Death.xlsx",
                      file.path(SUB, "Main_Tables/Main_Table_1_Death.xlsx"), "table"))
  ),
  list(
    id = "sdata2", output = "Supplementary Data 2 (full travel grid: infection risk, SAE/death/DALY benefit, risk, BRR, Pr(BRR>1))",
    script = "02_Scripts/22_export_infection_risk_full_grid.R",
    inputs = list(inp("data", "psa_grid_bra_travel_finite_setting.RData", "setting-level travel grid")),
    checks = list(chk("06_Results/infection_risk_full_grid.xlsx",
                      file.path(SUB, "Manuscript/Supplementary Data table 2.xlsx"), "grid"))
  ),
  list(
    id = "sfig1", output = "Supplementary Figure 1 (posterior predictive fit by state)",
    script = "02_Scripts/35_model_fit_figure_polished.R",
    inputs = list(inp("orv", "00_Data/0_2_Processed/observed_2022.RData", "observed 2022 reported cases by state (SINAN-derived)"),
                  inp("orv", "00_Data/0_2_Processed/posterior_finite_all.RData", "transmission-model posterior draws")),
    checks = list(chk("07_Final_Results/Main_Figures/figure_model_fit_by_state.png",
                      file.path(SUB, "Supplementary_Figures/Extended_fig_1.png"), "png"))
  ),
  list(
    id = "sfig2", output = "Supplementary Figure 2 (state map, attack rate, peak case rate)",
    script = "02_Scripts/16_figure1_epi_map_ar.R",
    inputs = list(inp("data", "setting_key.RData", "state -> transmission setting"),
                  inp("orv", "00_Data/0_2_Processed/observed_2022.RData", "observed 2022 reported cases"),
                  inp("data", "pop_by_state.RData", "state population"),
                  inp("data", "posterior_finite_all.RData", "transmission-model posterior draws")),
    checks = list(chk("06_Results/figure1_epi_map_ar.png",
                      file.path(SUB, "Supplementary_Figures/Extended_fig_2.png"), "png"))
  ),
  list(
    id = "sfig3", output = "Supplementary Figure 3 (travel SAE and death Pr(BRR>1) surfaces)",
    script = "02_Scripts/37_figure2_travel_risk_surface_sae_death.R",
    inputs = list(inp("data", "psa_grid_bra_travel_finite_setting.RData", "setting-level travel grid")),
    checks = list(chk("07_Final_Results/Supplementary_Figures/figure2_travel_risk_surfaces_SAE.png",
                      file.path(SUB, "Supplementary_Figures/figure2_travel_risk_surfaces_SAE.png"), "png"))
  ),
  list(
    id = "sfig4ab", output = "Supplementary Figure 4A-B (ORI SAE and death benefit-risk planes)",
    script = "02_Scripts/38_headline_brplane_sae_death_supplementary.R",
    inputs = list(inp("data", "draw_level_xy_serostatus_finite.RData", "ORI draw-level benefit/risk")),
    checks = list(chk("07_Final_Results/Supplementary_Figures/headline_benefit_risk_plane_SAE.png",
                      file.path(SUB, "Supplementary_Figures/headline_benefit_risk_plane_SAE.png"), "png"),
                  chk("07_Final_Results/Supplementary_Figures/headline_benefit_risk_plane_Death.png",
                      file.path(SUB, "Supplementary_Figures/headline_benefit_risk_plane_Death.png"), "png"))
  ),
  list(
    id = "sfig4cde", output = "Supplementary Figure 4C-E (ORI Pr(BRR>1) by campaign start week and coverage)",
    script = "02_Scripts/46_daly_weeksweep_probgt1_panelE.R",
    inputs = list(inp("orv", "01_Script/1_4_Analysis/finite_weeksweep_heatmap_probgt1.R", "plotting code sourced from the companion project"),
                  inp("orv", "00_Data/0_2_Processed/brr_weeksweep_summary_finite.RData", "week-sweep BRR summary (disease and infection blocking)"),
                  inp("orv", "00_Data/0_2_Processed/brr_weeksweep_summary_finite_ve0.RData", "week-sweep BRR summary (disease blocking only)")),
    checks = list(
      chk("07_Final_Results/Supplementary_Figures/brr_weeksweep_heatmap_probgt1_sae_base_mechanism_combined.png",
          file.path(SUB, "Supplementary_Figures/brr_weeksweep_heatmap_probgt1_sae_base_mechanism_combined.png"), "png"),
      chk("07_Final_Results/Supplementary_Figures/brr_weeksweep_heatmap_probgt1_death_base_mechanism_combined.png",
          file.path(SUB, "Supplementary_Figures/brr_weeksweep_heatmap_probgt1_death_base_mechanism_combined.png"), "png"),
      chk("07_Final_Results/Supplementary_Figures/brr_weeksweep_heatmap_probgt1_daly_base_mechanism_combined.png",
          file.path(SUB, "Supplementary_Figures/brr_weeksweep_heatmap_probgt1_daly_base_mechanism_combined.png"), "png"))
  )
)

# Upstream stages (Mode C). Listed for preflight and documentation; the
# runner executes them only with --mode full and explicit --allow-full.
UPSTREAM_STAGES <- list(
  list(id = "orv_fit", stage = "Transmission-model fitting (Stan) and outbreak-response simulations",
       where = "CHIK_ORV_impact (companion project, not in this repository)",
       produces = c("posterior_finite_all.RData", "sim_results_vc_ixchiq_model_finite.RData",
                    "postsim_vc_ixchiq_model_finite.RData", "combined_nnv_df_region_coverage_model_finite.RData",
                    "brr_weeksweep_summary_finite*.RData")),
  list(id = "setting_key", stage = "Transmission-setting classification",
       where = "02_Scripts/08_attack_rate_classification_finite.R", produces = "01_Data/setting_key.RData"),
  list(id = "ori", stage = "Outbreak-response benefit-risk (draw level)",
       where = "02_Scripts/03_brazil_all_draws_ori_v3.R (sources 01_setup.R, 02_setup_age_props.R)",
       produces = "01_Data/draw_level_xy_serostatus_finite.RData (+ all_draws_*_finite.RData)"),
  list(id = "travel", stage = "Travel benefit-risk (draw level and entry-day x duration grid)",
       where = "02_Scripts/06_brazil_travel_final_finite.R (sources 01_setup.R, 02_setup_age_props.R)",
       produces = c("01_Data/psa_df_bra_travel_finite.RData", "01_Data/psa_grid_bra_travel_finite.RData",
                    "01_Data/psa_grid_bra_travel_finite_setting.RData"))
)
