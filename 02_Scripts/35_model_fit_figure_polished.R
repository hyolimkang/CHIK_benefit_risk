# =============================================================================
# 35_model_fit_figure_polished.R
#
# Publication-ready (Nature Medicine style) version of the model-fit figure
# from CHIK_ORV_impact/01_Script/1_4_Analysis/age_struc_fitting_region_finite_2022.R
# (observed vs predicted weekly reported symptomatic cases, faceted by state).
#
# Rebuilt from ALREADY-SAVED fit output -- no Stan refit. The original script's
# `pred_all` (Median/Lower/Upper per week per state) was never itself saved to
# disk (only assembled transiently from df_xx_summ objects derived from the
# raw stanfit list, fits_prevacc_finite.RData, 2.5GB). posterior_finite_all.RData
# already contains, per state, `pred_cases` (1000 draws x 52 weeks -- the
# age-aggregated posterior-predictive REPORTED case count, i.e. exactly the
# quantity fit against observed_all$Observed) -- summarising that directly
# reproduces pred_all's Median/Lower/Upper without touching the 2.5GB stanfit
# object or resampling anything.
#
# Inputs (CHIK_ORV_impact project):
#   - 00_Data/0_2_Processed/observed_2022.RData      -> observed_all
#   - 00_Data/0_2_Processed/posterior_finite_all.RData -> posterior_<abbrev>$pred_cases
#
# Output: 07_Final_Results/Main_Figures/figure_model_fit_by_state.png
#         07_Final_Results/Main_Figures/figure_model_fit_by_state.pdf
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

ori_root <- "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"

load(file.path(ori_root, "00_Data/0_2_Processed/observed_2022.RData"))  # -> observed_all

message("Loading posterior_finite_all.RData (~1.3GB, this is slow)...")
load(file.path(ori_root, "00_Data/0_2_Processed/posterior_finite_all.RData"))

state_abbrev_map <- c(
  ce = "Ceará", bh = "Bahia", pa = "Paraíba", pn = "Pernambuco",
  rg = "Rio Grande do Norte", pi = "Piauí", ag = "Alagoas",
  tc = "Tocantins", mg = "Minas Gerais", se = "Sergipe", go = "Goiás"
)

pred_all <- dplyr::bind_rows(lapply(names(state_abbrev_map), function(ab) {
  post <- get(paste0("posterior_", ab))
  pc <- post$pred_cases  # 1000 draws x 52 weeks
  tibble::tibble(
    Week   = seq_len(ncol(pc)),
    Median = apply(pc, 2, median, na.rm = TRUE),
    Lower  = apply(pc, 2, quantile, probs = 0.025, na.rm = TRUE),
    Upper  = apply(pc, 2, quantile, probs = 0.975, na.rm = TRUE),
    Type   = "Predicted",
    region = state_abbrev_map[[ab]]
  )
}))

# ---- Nature Medicine-style polish -------------------------------------------
ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
grid_col      <- "#e7e6e0"
axis_col      <- "#c3c2b7"
col_predicted <- "#C0392B"

model_fit <- ggplot() +
  geom_ribbon(
    data = pred_all,
    aes(x = Week, ymin = Lower, ymax = Upper, fill = Type),
    alpha = 0.18, colour = NA
  ) +
  geom_line(
    data = pred_all,
    aes(x = Week, y = Median, colour = Type),
    linewidth = 0.8
  ) +
  geom_point(
    data = observed_all,
    aes(x = Week, y = Observed, shape = Type),
    size = 1.5, colour = ink_primary, alpha = 0.85
  ) +
  facet_wrap(~region, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = c(Predicted = col_predicted), name = NULL) +
  scale_fill_manual(values = c(Predicted = col_predicted), name = NULL) +
  scale_shape_manual(values = c(Observed = 16), name = NULL) +
  scale_x_continuous(name = "Epidemiological week (2022)", breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(name = "Reported symptomatic cases", labels = scales::comma,
                      expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Model fit: predicted versus observed weekly reported cases, by state",
    subtitle = "Points = observed reported cases; line/ribbon = posterior median and 95% credible interval (age-structured Bayesian SEIR fit, finite-history baseline immunity)"
  ) +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 1.8)),
    shape  = guide_legend(override.aes = list(size = 3.2))
  ) +
  theme_minimal(base_size = 17) +
  theme(
    text                 = element_text(colour = ink_primary),
    plot.title           = element_text(face = "bold", size = 20, colour = ink_primary),
    plot.subtitle        = element_text(size = 13, colour = ink_secondary, margin = margin(b = 12), lineheight = 1.15),
    plot.title.position  = "plot",
    axis.title           = element_text(size = 15, colour = ink_secondary),
    axis.text            = element_text(colour = ink_primary, size = 12.5),
    axis.line.x          = element_line(colour = axis_col, linewidth = 0.3),
    axis.ticks           = element_blank(),
    panel.grid.major     = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor     = element_blank(),
    panel.spacing        = unit(14, "pt"),
    strip.text           = element_text(face = "bold", size = 14.5, colour = ink_primary),
    strip.background     = element_blank(),
    legend.position      = "bottom",
    legend.title         = element_blank(),
    legend.text          = element_text(size = 14, colour = ink_secondary),
    legend.key           = element_blank(),
    plot.margin          = margin(16, 20, 14, 16)
  )

dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)

ggsave("07_Final_Results/Main_Figures/figure_model_fit_by_state.png", model_fit,
       width = 15, height = 10.5, dpi = 400, bg = "white")
ggsave("07_Final_Results/Main_Figures/figure_model_fit_by_state.pdf", model_fit,
       width = 15, height = 10.5, device = cairo_pdf)
ggsave("07_Final_Results/Main_Figures/figure_model_fit_by_state.jpg", model_fit,
       width = 15, height = 10.5, dpi = 400, bg = "white")

message("Saved: 07_Final_Results/Main_Figures/figure_model_fit_by_state.png / .pdf / .jpg")
