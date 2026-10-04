# =============================================================================
# 10_streamline_mechanism_figure.R
#
# "Why does the BRR result look the way it does?" -- a single composite,
# age-resolved figure connecting: fitted transmission dynamics (beta, R_eff)
# -> age-specific attack rate during the 2022 fitting window (from the
# posterior S compartment) -> age-specific averted burden -> age-specific
# vaccine-attributable risk -> age-specific BRR.
#
# Required inputs (must exist on disk before running):
#   01_Data/draw_level_xy_serostatus_finite.RData  (saved from a DATA_MODE <-
#     "flat" run of 03_...v3.R: save(draw_level_xy_serostatus, file =
#     "01_Data/draw_level_xy_serostatus_finite.RData"))
#   01_Data/posterior_finite_all.RData      (the finite-history Stan fits,
#     loaded as 11 individual objects posterior_ce/posterior_bh/.../posterior_go
#     -- NOT a single posterior_list. Reassembled into posterior_list below
#     using the same region-name mapping as 03_brazil_all_draws_ori_v3.R.
#     Confirmed fields (posterior_ce, 1000 draws): beta (1000x54, unused --
#     wrong week count), beta_observed (1000x52, the one to use), R_eff
#     (1000x52), S_pred (1000x20x55, age-resolved), susceptible_fraction
#     (1000x52, NOT age-resolved -- aggregate only, not used here).)
#
# Output: 06_Results/mechanism_streamline_<region>.png (one per region) and
# a combined small-multiple version across all regions.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(purrr)
library(tibble)
library(scales)
library(ggrepel)

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")

# ---- 0) Load ---------------------------------------------------------------
load("01_Data/draw_level_xy_serostatus_finite.RData")
load("01_Data/posterior_finite_all.RData")

posterior_list <- list(
  "Ceará"               = posterior_ce,
  "Bahia"               = posterior_bh,
  "Paraíba"             = posterior_pa,
  "Pernambuco"          = posterior_pn,
  "Rio Grande do Norte" = posterior_rg,
  "Piauí"               = posterior_pi,
  "Alagoas"             = posterior_ag,
  "Tocantins"           = posterior_tc,
  "Minas Gerais"        = posterior_mg,
  "Sergipe"             = posterior_se,
  "Goiás"               = posterior_go
)

age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
                 mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
                 mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
                 mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))

age_bin_agecat <- c(rep("1-11", 4), "12-17", rep("18-64", 10), rep("65+", 5))
stopifnot(length(age_bin_agecat) == 20)

qsum <- function(x, probs = c(0.5, 0.025, 0.975)) {
  q <- stats::quantile(x, probs, na.rm = TRUE)
  setNames(as.list(q), c("med", "lo", "hi"))
}

# ---- Nature-style theme -----------------------------------------------------
theme_nm <- function(base_size = 10) {
  theme_minimal(base_size = base_size, base_family = "") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      axis.line = element_line(linewidth = 0.3, colour = "grey30"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey30"),
      strip.text = element_text(face = "bold", size = rel(0.95)),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      legend.position = "bottom",
      legend.title = element_text(size = rel(0.85)),
      legend.text = element_text(size = rel(0.8))
    )
}

pal_mechanism <- c("Disease blocking only" = "#4C72B0",
                    "Disease and infection blocking" = "#C44E52")
pal_outcome <- c("DALY" = "#2A9D8F", "SAE" = "#E9762B", "Death" = "#6C5B7B")

# =============================================================================
# Panel A -- fitted weekly transmission rate (beta) + R_eff, 2022 season
# =============================================================================
build_panel_beta <- function(region) {
  p <- posterior_list[[region]]
  beta_df <- as.data.frame(p$beta_observed)
  names(beta_df) <- paste0("w", seq_len(ncol(beta_df)))
  beta_long <- beta_df %>%
    mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = "week", values_to = "beta") %>%
    mutate(week = as.integer(sub("w", "", week)))

  reff_df <- as.data.frame(p$R_eff)
  names(reff_df) <- paste0("w", seq_len(ncol(reff_df)))
  reff_long <- reff_df %>%
    mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = "week", values_to = "R_eff") %>%
    mutate(week = as.integer(sub("w", "", week)))

  beta_summ <- beta_long %>% group_by(week) %>%
    summarise(med = median(beta, na.rm = TRUE),
              lo = quantile(beta, 0.025, na.rm = TRUE),
              hi = quantile(beta, 0.975, na.rm = TRUE), .groups = "drop")
  reff_summ <- reff_long %>% group_by(week) %>%
    summarise(med = median(R_eff, na.rm = TRUE),
              lo = quantile(R_eff, 0.025, na.rm = TRUE),
              hi = quantile(R_eff, 0.975, na.rm = TRUE), .groups = "drop")

  scale_r <- max(beta_summ$hi, na.rm = TRUE) / max(reff_summ$hi, na.rm = TRUE)

  ggplot() +
    geom_ribbon(data = beta_summ, aes(week, ymin = lo, ymax = hi),
                fill = "grey70", alpha = 0.35) +
    geom_line(data = beta_summ, aes(week, med), colour = "grey20", linewidth = 0.6) +
    geom_line(data = reff_summ, aes(week, med * scale_r), colour = "#C44E52",
              linewidth = 0.6, linetype = "22") +
    geom_hline(yintercept = 1 * scale_r, linetype = "dotted", colour = "#C44E52",
               linewidth = 0.4) +
    scale_y_continuous(
      name = "Weekly transmission rate (\u03b2)",
      sec.axis = sec_axis(~ . / scale_r, name = "Effective reproduction number (R\u2091\u2091, dashed)")
    ) +
    labs(x = "Fitting week (2022 season)",
         title = paste0(region, " \u2014 fitted transmission dynamics")) +
    theme_nm()
}

# =============================================================================
# Panel B -- age-specific BASELINE immunity at the start of the fitting
# window (week 1), i.e. how the sero_finite prior actually landed in the
# fitted model. Self-normalised: N_age is reconstructed per draw as
# S+E+I+R at week 1 (no external age-population vector needed), so
# 1 - S[week1]/N_age is the immune fraction at t0.
#
# NOTE: using proportional depletion (1 - S[last week]/S[week1]) instead
# comes out flat across ages in this model, because per-susceptible
# infection hazard isn't age-differentiated -- a uniform hazard depletes
# every age group's *pool* by the same proportion regardless of its
# starting size. That's a real feature of the model, not a bug, but it
# doesn't show baseline immunity -- hence using week-1 levels instead.
# =============================================================================
build_panel_age_attack <- function(region) {
  p <- posterior_list[[region]]
  S <- p$S_pred; E <- p$E_pred; I <- p$I_pred; R <- p$R_pred  # draw x age x week

  N_age <- S[, , 1] + E[, , 1] + I[, , 1] + R[, , 1]  # draw x age
  immune_frac <- 1 - S[, , 1] / N_age                 # draw x age

  im_df <- as.data.frame(immune_frac)
  names(im_df) <- paste0("a", seq_len(ncol(im_df)))
  im_long <- im_df %>% mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = "age_bin", values_to = "immune_frac") %>%
    mutate(age_bin = as.integer(sub("a", "", age_bin)),
           age = age_groups[age_bin]) %>%
    filter(age_bin >= 6)  # drop bins covering ages 1-17 (bins 1-5); manuscript scope is adults only

  im_summ <- im_long %>% group_by(age_bin, age) %>%
    summarise(med = median(immune_frac, na.rm = TRUE),
              lo = quantile(immune_frac, 0.025, na.rm = TRUE),
              hi = quantile(immune_frac, 0.975, na.rm = TRUE), .groups = "drop")

  ggplot(im_summ, aes(age, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4C72B0", alpha = 0.25) +
    geom_line(colour = "#4C72B0", linewidth = 0.7) +
    geom_point(colour = "#4C72B0", size = 1.2) +
    scale_y_continuous(labels = scales::percent, name = "Baseline immune fraction\n(start of 2022 fitting window)",
                        limits = c(0, NA)) +  # anchor at 0 so the near-flat
    # plateau (true variation across adult ages is <0.2 percentage points --
    # MCMC noise, not a real age effect) doesn't get visually exaggerated by
    # auto-scaling into a tiny y-range
    labs(x = "Age (years)", title = paste0(region, " \u2014 age-specific baseline immunity")) +
    theme_nm()
}

# =============================================================================
# Panels C/D/E -- age-specific averted burden, risk, BRR (from
# draw_level_xy_serostatus; cov50 pipeline default, RR_seropos = 0 = base risk
# unless *_adj requested)
# =============================================================================
prep_age_summary <- function(region, use_adj = FALSE) {
  df <- draw_level_xy_serostatus %>%
    filter(Region == region, RR_seropos == 0,
           AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64"))  # adults only; manuscript excludes ages 1-17

  averted <- df %>%
    group_by(AgeCat, outcome, VE_label, setting) %>%
    summarise(med = median(averted_10k, na.rm = TRUE),
              lo = quantile(averted_10k, 0.025, na.rm = TRUE),
              hi = quantile(averted_10k, 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(panel = "Averted burden\n(per 10,000 vaccinated)")

  risk_col <- if (use_adj) "x_10k_adj" else "x_10k_base"
  risk <- df %>%
    group_by(AgeCat, outcome, VE_label, setting) %>%
    summarise(med = median(.data[[risk_col]], na.rm = TRUE),
              lo = quantile(.data[[risk_col]], 0.025, na.rm = TRUE),
              hi = quantile(.data[[risk_col]], 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(panel = "Vaccine-attributable risk\n(per 10,000 vaccinated)")

  brr_col <- if (use_adj) "brr_adj" else "brr_base"
  brr <- df %>%
    group_by(AgeCat, outcome, VE_label, setting) %>%
    summarise(med = median(.data[[brr_col]], na.rm = TRUE),
              lo = quantile(.data[[brr_col]], 0.025, na.rm = TRUE),
              hi = quantile(.data[[brr_col]], 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(panel = "Benefit-risk ratio")

  list(averted = averted, risk = risk, brr = brr)
}

# setting is a region-level constant (each region is classified Low/
# Moderate/High), so faceting by it within a single-region plot is a no-op
# facet -- omitted here. It's shown once in the figure title instead.
build_panel_averted <- function(region_summary) {
  ggplot(region_summary$averted,
         aes(AgeCat, med, colour = outcome, group = interaction(outcome, VE_label))) +
    geom_pointrange(aes(ymin = lo, ymax = hi, shape = VE_label),
                     position = position_dodge(width = 0.5), fatten = 2) +
    scale_colour_manual(values = pal_outcome, name = "Outcome") +
    scale_shape_discrete(name = "Mechanism") +
    scale_y_continuous(name = "Averted per 10,000 vaccinated") +
    labs(x = "Age group", title = "Age-specific benefit") +
    theme_nm()
}

build_panel_risk <- function(region_summary) {
  ggplot(region_summary$risk,
         aes(AgeCat, med, colour = outcome, group = interaction(outcome, VE_label))) +
    geom_pointrange(aes(ymin = lo, ymax = hi, shape = VE_label),
                     position = position_dodge(width = 0.5), fatten = 2) +
    scale_colour_manual(values = pal_outcome, name = "Outcome") +
    scale_shape_discrete(name = "Mechanism") +
    scale_y_continuous(name = "Attributable per 10,000 vaccinated") +
    labs(x = "Age group", title = "Age-specific vaccine-attributable risk") +
    theme_nm()
}

build_panel_brr <- function(region_summary) {
  ggplot(region_summary$brr,
         aes(AgeCat, med, colour = VE_label, group = interaction(outcome, VE_label))) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(ymin = lo, ymax = hi, shape = outcome),
                     position = position_dodge(width = 0.5), fatten = 2) +
    scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
    scale_shape_discrete(name = "Outcome") +
    scale_y_log10(name = "BRR (log scale)") +  # BRR is a ratio spanning orders of
    # magnitude (e.g. Ceará DALY/infection-blocking CI reaches ~3000) -- a linear
    # axis squashes every other point near 0, so this panel keeps log scale
    # while C/D (benefit, risk) stay linear per your request.
    labs(x = "Age group", title = "Age-specific benefit-risk ratio") +
    theme_nm()
}

# =============================================================================
# Assemble one composite per region
# =============================================================================
build_mechanism_figure <- function(region, use_adj = FALSE) {
  p_beta   <- build_panel_beta(region)
  p_attack <- build_panel_age_attack(region)
  rs       <- prep_age_summary(region, use_adj = use_adj)
  p_avert  <- build_panel_averted(rs)
  p_risk   <- build_panel_risk(rs)
  p_brr    <- build_panel_brr(rs)

  region_setting <- unique(draw_level_xy_serostatus$setting[draw_level_xy_serostatus$Region == region])[1]

  top    <- p_beta + p_attack + plot_layout(ncol = 2)
  bottom <- p_avert / p_risk / p_brr

  (top / bottom) +
    plot_layout(heights = c(1, 2.4), guides = "collect") +
    plot_annotation(
      title = paste0(region, " (", region_setting, " transmission setting): from transmission dynamics to benefit-risk"),
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    ) &
    theme(legend.position = "bottom")
}

# ---- Save each panel individually (for one-at-a-time review) ---------------
save_individual_panels <- function(region, use_adj = FALSE, out_dir = "06_Results") {
  p_beta   <- build_panel_beta(region)
  p_immun  <- build_panel_age_attack(region)
  rs       <- prep_age_summary(region, use_adj = use_adj)
  p_avert  <- build_panel_averted(rs)
  p_risk   <- build_panel_risk(rs)
  p_brr    <- build_panel_brr(rs)

  panels <- list(A_beta = p_beta, B_immunity = p_immun, C_averted = p_avert,
                 D_risk = p_risk, E_brr = p_brr)
  for (nm in names(panels)) {
    ggsave(sprintf("%s/panel_%s_%s.png", out_dir, nm, region),
           panels[[nm]], width = 7, height = 5, dpi = 300, bg = "white")
  }
  message("Saved 5 individual panels for ", region, " to ", out_dir, "/panel_*_", region, ".png")
  invisible(panels)
}

# ---- Run for one region first as a sanity check -----------------------------
region_check <- "Ceará"
save_individual_panels(region_check)

fig_check <- build_mechanism_figure(region_check)
ggsave(sprintf("06_Results/mechanism_streamline_%s.png", region_check),
       fig_check, width = 11, height = 13, dpi = 300, bg = "white")

message("Saved: 06_Results/mechanism_streamline_", region_check, ".png")
message("Once this looks right, loop build_mechanism_figure() (or save_individual_panels()) over all 11 regions.")

# =============================================================================
# HEADLINE FIGURE -- 11-state synthesis: does baseline immunity level predict
# the benefit-risk result? One point per state x mechanism x age group.
# x = state's baseline immune fraction (adult ages, flat plateau value,
#     reconstructed the same way as build_panel_age_attack but summarised to
#     a single number per state -- one draw-level mean across adult age bins,
#     then median + 95% CI across draws).
# y = averted burden (fig 1) / BRR (fig 2), per state, DALY outcome, cov50.
# =============================================================================
compute_state_immunity <- function(region) {
  p <- posterior_list[[region]]
  S <- p$S_pred; E <- p$E_pred; I <- p$I_pred; R <- p$R_pred
  N_age <- S[, , 1] + E[, , 1] + I[, , 1] + R[, , 1]
  immune_frac <- 1 - S[, , 1] / N_age            # draw x age
  adult_bins <- 6:ncol(immune_frac)              # ages >= 18 (bin 6 = 18-19)
  im_by_draw <- rowMeans(immune_frac[, adult_bins], na.rm = TRUE)
  tibble::tibble(
    Region = region,
    immunity_med = median(im_by_draw, na.rm = TRUE),
    immunity_lo  = quantile(im_by_draw, 0.025, na.rm = TRUE),
    immunity_hi  = quantile(im_by_draw, 0.975, na.rm = TRUE)
  )
}

state_immunity <- purrr::map_dfr(names(posterior_list), compute_state_immunity)

build_state_synthesis_data <- function(outcome_sel = "DALY", use_adj = FALSE) {
  brr_col     <- if (use_adj) "brr_adj" else "brr_base"
  risk_col    <- if (use_adj) "x_10k_adj" else "x_10k_base"

  state_outcome <- draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel,
           AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64")) %>%
    group_by(Region, VE_label, AgeCat) %>%
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
    )

  state_outcome %>% left_join(state_immunity, by = "Region")
}

log_num_labels <- scales::label_number(accuracy = NULL, big.mark = ",")

build_headline_averted <- function(outcome_sel = "DALY", use_adj = FALSE) {
  d <- build_state_synthesis_data(outcome_sel, use_adj)
  ggplot(d, aes(immunity_med, averted_med, colour = VE_label, shape = AgeCat)) +
    geom_errorbar(aes(ymin = averted_lo, ymax = averted_hi), width = 0, alpha = 0.4) +
    geom_errorbarh(aes(xmin = immunity_lo, xmax = immunity_hi), height = 0, alpha = 0.4) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
    scale_shape_manual(values = c("18-64" = 16, "65+" = 17), name = "Age group") +
    scale_x_continuous(labels = scales::percent, name = "State baseline immune fraction (adults)") +
    scale_y_continuous(labels = log_num_labels, name = paste(outcome_sel, "averted per 10,000 vaccinated")) +
    labs(title = paste0("Does baseline immunity predict benefit? (", outcome_sel, ", 11 states)")) +
    theme_nm()
}

build_headline_brr <- function(outcome_sel = "DALY", use_adj = FALSE) {
  d <- build_state_synthesis_data(outcome_sel, use_adj)
  ggplot(d, aes(immunity_med, brr_med, colour = VE_label, shape = AgeCat)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_errorbar(aes(ymin = brr_lo, ymax = brr_hi), width = 0, alpha = 0.4) +
    geom_errorbarh(aes(xmin = immunity_lo, xmax = immunity_hi), height = 0, alpha = 0.4) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
    scale_shape_manual(values = c("18-64" = 16, "65+" = 17), name = "Age group") +
    scale_x_continuous(labels = scales::percent, name = "State baseline immune fraction (adults)") +
    scale_y_log10(labels = log_num_labels, breaks = c(0.1, 1, 10, 100, 1000),
                  name = paste(outcome_sel, "benefit-risk ratio (log scale)")) +
    labs(title = paste0("Does baseline immunity predict BRR? (", outcome_sel, ", 11 states)")) +
    theme_nm()
}

p_headline_averted <- build_headline_averted("DALY")
p_headline_brr     <- build_headline_brr("DALY")

ggsave("06_Results/headline_immunity_vs_averted_DALY.png", p_headline_averted,
       width = 8, height = 5.5, dpi = 300, bg = "white")
ggsave("06_Results/headline_immunity_vs_brr_DALY.png", p_headline_brr,
       width = 8, height = 5.5, dpi = 300, bg = "white")

message("Saved headline 11-state synthesis figures: ",
        "06_Results/headline_immunity_vs_{averted,brr}_DALY.png")

# =============================================================================
# BENEFIT-RISK PLANE -- x = risk, y = benefit, both per 10,000 vaccinated,
# log-log, diagonal = BRR 1. One point per state, coloured by age group,
# faceted by mechanism (so the two mechanisms no longer sit on the same
# vertical line). This directly shows the 65+ vs 18-64 story: does 65+ sit
# further right (higher risk) relative to how far up it sits (higher
# benefit) than 18-64 -- i.e. is the BRR gap coming from risk outpacing
# benefit, not from benefit being lower.
# =============================================================================
build_setting_synthesis_data <- function(outcome_sel = "DALY", use_adj = FALSE) {
  brr_col  <- if (use_adj) "brr_adj" else "brr_base"
  risk_col <- if (use_adj) "x_10k_adj" else "x_10k_base"

  # Pool draws across every state within a setting (Low/Moderate/High) --
  # this is the paper's own primary reporting unit (matches
  # brr_table_final_long), not individual states. 11 states -> 3 groups.
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
    # Death risk/benefit can hit exactly 0 at the 2.5th percentile (rare
    # events, e.g. Reviewer #4's "zero death risk" point) -- log10(0) = -Inf,
    # which breaks position_dodge()'s interval-overlap check. Floor at a
    # small epsilon (1 per million vaccinated) so the log axis stays finite;
    # the point/CI still reads as "very close to zero".
    mutate(across(c(averted_med, averted_lo, averted_hi,
                     risk_med, risk_lo, risk_hi,
                     brr_med, brr_lo, brr_hi),
                   ~ pmax(.x, 1e-4)))
}

build_benefit_risk_plane <- function(outcome_sel = "DALY", use_adj = FALSE,
                                      shared_rng = NULL) {
  d <- build_setting_synthesis_data(outcome_sel, use_adj) %>%
    # Small deterministic horizontal offset by setting on the log10 scale.
    # This separates otherwise coincident point/error-bar glyphs; it is a
    # display adjustment only, not a change to the underlying estimates.
    mutate(
      x_offset = 10^case_when(
        setting == "Low"      ~ -0.06,
        setting == "Moderate" ~  0,
        setting == "High"     ~  0.06,
        TRUE                  ~  0
      ),
      risk_med_plot = risk_med * x_offset,
      risk_lo_plot  = risk_lo  * x_offset,
      risk_hi_plot  = risk_hi  * x_offset
    )

  rng <- shared_rng %||% pmax(range(c(d$risk_med, d$averted_med, d$risk_lo, d$risk_hi,
                                       d$averted_lo, d$averted_hi), na.rm = TRUE), 1e-3)

  risk_label <- if (use_adj) paste0(outcome_sel, " caused by vaccination (per 10,000 vaccinated individuals)\n(serostatus-adjusted, log scale)") else
    paste0(outcome_sel, " caused by vaccination (per 10,000 vaccinated individuals)\n(base, log scale)")

  ggplot(d, aes(risk_med_plot, averted_med, colour = AgeCat, shape = setting)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(
      aes(ymin = averted_lo, ymax = averted_hi,
          linetype = "95% UI"),
      width = 0, alpha = 0.3,
      show.legend = c(colour = FALSE, linetype = TRUE)
    ) +
    geom_errorbarh(
      aes(xmin = risk_lo_plot, xmax = risk_hi_plot,
          linetype = "95% UI"),
      height = 0, alpha = 0.3,
      show.legend = c(colour = FALSE, linetype = TRUE)
    ) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~VE_label) +
    coord_equal(xlim = rng, ylim = rng) +
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"), name = "Age group") +
    scale_shape_manual(values = c(Low = 15, Moderate = 17, High = 16), name = "Setting") +
    scale_linetype_manual(
      values = c("95% UI" = "solid"),
      name = "Bars"
    ) +
    scale_x_log10(labels = log_num_labels, name = risk_label) +
    scale_y_log10(labels = log_num_labels, name = paste(outcome_sel, "averted per 10,000 vaccinated (log scale)")) +
    guides(
      shape = guide_legend(order = 1),
      colour = guide_legend(order = 2),
      linetype = guide_legend(
        order = 3,
        override.aes = list(colour = "grey40", alpha = 0.6)
      )
    ) +
    theme_nm()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Base and serostatus-adjusted risk share one y-range (benefit doesn't change
# between the two) so the two rows are visually comparable.
build_and_save_brplane <- function(outcome_sel) {
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
      title = paste0("Benefit-risk plane by age group, mechanism, setting, and serostatus adjustment (", outcome_sel, ")"),
      tag_levels = list(c("Base risk", "Serostatus-adjusted risk")),
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    ) &
    theme(legend.position = "bottom")

  out_file <- sprintf("06_Results/headline_benefit_risk_plane_%s.png", outcome_sel)
  ggsave(out_file, p_combined, width = 10, height = 10, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p_combined)
}

for (oc in c("DALY", "Death", "SAE")) build_and_save_brplane(oc)

# =============================================================================
# HEADLINE DECISION MAP -- probability(BRR > 1) by setting x age group,
# faceted by mechanism. Answers Reviewer #5's explicit request directly: "a
# summary figure ... clearly highlighting the combinations of transmission
# setting, vaccination strategy, and age group for which vaccination provides
# a net benefit." Complements (not replaces) the benefit-risk plane above,
# which shows absolute magnitudes; this shows the decision + confidence in
# one glance, at the paper's actual reporting resolution (setting, not state).
# =============================================================================
build_decision_map_data <- function(outcome_sel = "DALY", use_adj = FALSE) {
  brr_col <- if (use_adj) "brr_adj" else "brr_base"
  draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel, AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64")) %>%
    group_by(setting, AgeCat, VE_label) %>%
    summarise(
      prob_gt1 = mean(.data[[brr_col]] > 1, na.rm = TRUE),
      brr_med  = median(.data[[brr_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(setting = factor(setting, levels = c("Low", "Moderate", "High")))
}

fmt_brr <- function(x) ifelse(x >= 100, sprintf("%.0f", x),
                        ifelse(x >= 10,  sprintf("%.1f", x), sprintf("%.2f", x)))

build_decision_map <- function(outcome_sel = "DALY", use_adj = FALSE) {
  d <- build_decision_map_data(outcome_sel, use_adj)
  ggplot(d, aes(setting, AgeCat, fill = prob_gt1)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(sprintf("%.0f%%", prob_gt1 * 100), "\nBRR=", fmt_brr(brr_med))),
              size = 3, lineheight = 0.9, colour = "grey10") +
    facet_wrap(~VE_label) +
    scale_fill_gradient2(low = "#C44E52", mid = "white", high = "#55A868",
                          midpoint = 0.5, limits = c(0, 1),
                          labels = scales::percent, name = "Pr(BRR > 1)") +
    labs(x = "Setting", y = "Age group") +
    theme_nm() +
    theme(panel.grid = element_blank())
}

# Base risk and serostatus-adjusted risk shown as a 2-row composite, same
# convention as build_and_save_brplane() -- risk denominator changes between
# rows, benefit (numerator) does not.
build_and_save_decision_map <- function(outcome_sel) {
  p_base <- build_decision_map(outcome_sel, use_adj = FALSE)
  p_adj  <- build_decision_map(outcome_sel, use_adj = TRUE)

  p_combined <- (p_base / p_adj) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(outcome_sel, " — probability of net benefit"),
      caption = paste0("Tile = probability that ", outcome_sel,
                        " BRR exceeds 1 across 1000 posterior draws; label shows that probability and the median BRR."),
      tag_levels = list(c("Base risk", "Serostatus-adjusted risk")),
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    ) &
    theme(legend.position = "bottom")

  out_file <- sprintf("06_Results/headline_decision_map_%s.png", outcome_sel)
  ggsave(out_file, p_combined, width = 7, height = 9.5, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p_combined)
}

for (oc in c("DALY", "Death", "SAE")) build_and_save_decision_map(oc)

# =============================================================================
# MECHANISM EXPLAINER -- one representative state per setting (Low/Moderate/
# High), showing age-specific baseline immunity (top row) directly above the
# age-specific DALY BRR it produces (bottom row). This is the "why" companion
# to the decision map above (the "what"): it makes the causal chain from
# baseline-immunity assumption -> age-specific benefit -> BRR visible, instead
# of only reporting the resulting numbers.
# =============================================================================
build_panel_brr_daly <- function(region_summary) {
  d <- region_summary$brr %>% filter(outcome == "DALY")
  ggplot(d, aes(AgeCat, med, colour = VE_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.3), fatten = 2) +
    scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
    scale_y_log10(name = "DALY BRR (log scale)") +
    labs(x = "Age group") +
    theme_nm()
}

build_mechanism_explainer <- function(regions = c("Bahia", "Pernambuco", "Ceará"),
                                       setting_labels = c("Low", "Moderate", "High")) {
  panels <- purrr::map2(regions, setting_labels, function(region, setting_lab) {
    p_immun <- build_panel_age_attack(region) +
      ggtitle(paste0(region, " (", setting_lab, ")"))
    rs <- prep_age_summary(region, use_adj = FALSE)
    p_brr <- build_panel_brr_daly(rs)
    list(immun = p_immun, brr = p_brr)
  })

  top_row <- wrap_plots(purrr::map(panels, "immun"), nrow = 1)
  bot_row <- wrap_plots(purrr::map(panels, "brr"), nrow = 1, guides = "collect")

  p_combined <- (top_row / bot_row) +
    plot_annotation(
      title = "From baseline immunity to benefit-risk ratio: one representative state per transmission setting",
      caption = paste(
        "Top row: age-specific baseline immunity (start of 2022 fitting window). Bottom row: resulting age-specific DALY BRR by mechanism (dashed line = BRR 1, log scale).",
        "Same underlying susceptible-pool level, read through each state's own fitted transmission intensity, produces a different age gradient in averted burden --",
        "translated through a fixed vaccine-attributable risk into the BRR pattern shown in the decision map.",
        sep = "\n"
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    ) &
    theme(legend.position = "bottom")

  ggsave("06_Results/mechanism_explainer_by_setting.png", p_combined,
         width = 13, height = 8, dpi = 300, bg = "white")
  message("Saved: 06_Results/mechanism_explainer_by_setting.png")
  invisible(p_combined)
}

build_mechanism_explainer()

# =============================================================================
# STATE-LEVEL BASELINE IMMUNITY / SUSCEPTIBLE POOL SUMMARY -- one point per
# state x age band, to see (a) whether the 11 states form natural tiers on
# baseline immunity itself (not just downstream attack rate), and (b) whether
# the current AR-based Low/Moderate/High setting classification (colour)
# tracks that. Reuses the same S/E/I/R week-1 computation as
# build_panel_age_attack(), summarised to one number per state x age band
# instead of a full age curve.
# =============================================================================
compute_state_immunity_summary <- function(region) {
  p <- posterior_list[[region]]
  S <- p$S_pred; E <- p$E_pred; I <- p$I_pred; R <- p$R_pred  # draw x age x week
  N_age <- S[, , 1] + E[, , 1] + I[, , 1] + R[, , 1]          # draw x age
  immune_frac <- 1 - S[, , 1] / N_age                         # draw x age

  age_band <- list("18-64" = 6:15, "65+" = 16:20)
  purrr::imap_dfr(age_band, function(bins, band_name) {
    draw_mean <- rowMeans(immune_frac[, bins, drop = FALSE])  # one number per draw
    tibble(
      Region = region,
      AgeCat = band_name,
      immune_med = median(draw_mean, na.rm = TRUE),
      immune_lo  = quantile(draw_mean, 0.025, na.rm = TRUE),
      immune_hi  = quantile(draw_mean, 0.975, na.rm = TRUE)
    )
  })
}

state_immunity_summary <- purrr::map_dfr(names(posterior_list), compute_state_immunity_summary) %>%
  mutate(susceptible_med = 1 - immune_med,
         susceptible_lo  = 1 - immune_hi,
         susceptible_hi  = 1 - immune_lo)

# Attach each state's current AR-based transmission setting for comparison.
region_setting_map <- draw_level_xy_serostatus %>% distinct(Region, setting)
state_immunity_summary <- state_immunity_summary %>%
  left_join(region_setting_map, by = "Region") %>%
  mutate(setting = factor(setting, levels = c("Low", "Moderate", "High")))

# Order states by 18-64 susceptible pool (ascending) so natural tiers, if any,
# are visible without relying on the AR-based colour grouping.
state_order <- state_immunity_summary %>%
  filter(AgeCat == "18-64") %>%
  arrange(susceptible_med) %>%
  pull(Region)
state_immunity_summary <- state_immunity_summary %>%
  mutate(Region = factor(Region, levels = state_order))

p_state_immunity <- ggplot(state_immunity_summary,
                            aes(susceptible_med, Region, colour = setting, shape = AgeCat)) +
  geom_errorbarh(aes(xmin = susceptible_lo, xmax = susceptible_hi),
                  height = 0, position = position_dodge(width = 0.5), alpha = 0.6) +
  geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(Low = "#55A868", Moderate = "#DD8452", High = "#C44E52"),
                       name = "AR-based setting\n(current classification)") +
  scale_shape_manual(values = c("18-64" = 16, "65+" = 17), name = "Age group") +
  scale_x_continuous(labels = scales::percent,
                      name = "Susceptible pool remaining\n(start of 2022 fitting window)") +
  labs(y = NULL,
       title = "Baseline immunity / susceptible pool by state (finite-history model)",
       caption = paste("One point per state x age band; error bars = 95% UI across 1000 posterior draws.",
                        "Colour = current AR-based Low/Moderate/High transmission-setting classification, for comparison.",
                        "States ordered by 18-64 susceptible pool (ascending).", sep = "\n")) +
  theme_nm()

ggsave("06_Results/state_baseline_immunity_summary.png", p_state_immunity,
       width = 8, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/state_baseline_immunity_summary.png")

# =============================================================================
# BASELINE IMMUNITY TIER (Low/Mid/High) x BENEFIT/BRR -- discretises the 11
# states into 3 baseline-immunity tiers (terciles of 18-64 immune fraction),
# analogous to the discrete seroprevalence-tier design (e.g. SP9 20/40/60/80%)
# used in dengue vaccine-impact benefit figures, instead of the continuous
# state-level scatter above. Two-row composite: benefit (top) and BRR
# (bottom), by tier x mechanism x age group x serostatus adjustment.
# =============================================================================
immunity_tier_map <- state_immunity_summary %>%
  filter(AgeCat == "18-64") %>%
  mutate(Region = as.character(Region),
         immunity_tier = cut(immune_med,
                              breaks = quantile(immune_med, probs = c(0, 1/3, 2/3, 1)),
                              labels = c("Low", "Mid", "High"),
                              include.lowest = TRUE)) %>%
  select(Region, immunity_tier)

build_tier_data <- function(outcome_sel = "DALY", use_adj = FALSE) {
  brr_col <- if (use_adj) "brr_adj" else "brr_base"
  draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel, AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64")) %>%
    inner_join(immunity_tier_map, by = "Region") %>%
    group_by(immunity_tier, AgeCat, VE_label) %>%
    summarise(
      averted_med = median(averted_10k, na.rm = TRUE),
      averted_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
      averted_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),
      brr_med = median(.data[[brr_col]], na.rm = TRUE),
      brr_lo  = quantile(.data[[brr_col]], 0.025, na.rm = TRUE),
      brr_hi  = quantile(.data[[brr_col]], 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(adj = if (use_adj) "Serostatus-adjusted" else "Base",
           immunity_tier = factor(immunity_tier, levels = c("Low", "Mid", "High")))
}

build_tier_figure <- function(outcome_sel = "DALY") {
  d <- bind_rows(build_tier_data(outcome_sel, FALSE), build_tier_data(outcome_sel, TRUE)) %>%
    mutate(adj = factor(adj, levels = c("Base", "Serostatus-adjusted")))

  p_benefit <- ggplot(d, aes(immunity_tier, averted_med, colour = AgeCat,
                              group = interaction(AgeCat, adj))) +
    geom_line(aes(linetype = adj), position = position_dodge(width = 0.3)) +
    geom_pointrange(aes(ymin = averted_lo, ymax = averted_hi, shape = adj),
                     position = position_dodge(width = 0.3), fatten = 2) +
    facet_wrap(~VE_label) +
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"), name = "Age group") +
    scale_linetype_manual(values = c("Base" = "solid", "Serostatus-adjusted" = "22"), name = "Risk basis") +
    scale_shape_manual(values = c("Base" = 16, "Serostatus-adjusted" = 17), name = "Risk basis") +
    labs(x = NULL, y = paste(outcome_sel, "averted per 10,000 vaccinated")) +
    theme_nm()

  p_brr <- ggplot(d, aes(immunity_tier, brr_med, colour = AgeCat,
                          group = interaction(AgeCat, adj))) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_line(aes(linetype = adj), position = position_dodge(width = 0.3)) +
    geom_pointrange(aes(ymin = brr_lo, ymax = brr_hi, shape = adj),
                     position = position_dodge(width = 0.3), fatten = 2) +
    facet_wrap(~VE_label) +
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"), name = "Age group") +
    scale_linetype_manual(values = c("Base" = "solid", "Serostatus-adjusted" = "22"), name = "Risk basis") +
    scale_shape_manual(values = c("Base" = 16, "Serostatus-adjusted" = 17), name = "Risk basis") +
    scale_y_log10(labels = log_num_labels, name = paste(outcome_sel, "BRR (log scale)")) +
    labs(x = "Baseline immunity tier (Low = most susceptible pool remaining, High = least)") +
    theme_nm()

  p_combined <- (p_benefit / p_brr) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(outcome_sel, " — benefit and BRR by baseline-immunity tier"),
      caption = "Tiers = terciles of 18-64 baseline immune fraction across the 11 states. Points/lines = median with 95% UI across 1000 posterior draws.",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    ) &
    theme(legend.position = "bottom")

  ggsave(sprintf("06_Results/immunity_tier_benefit_brr_%s.png", outcome_sel), p_combined,
         width = 9, height = 9, dpi = 300, bg = "white")
  message("Saved: 06_Results/immunity_tier_benefit_brr_", outcome_sel, ".png")
  invisible(p_combined)
}

build_tier_figure("DALY")

# =============================================================================
# SUPPLEMENTARY FIGURE -- baseline immunity vs fitted transmission intensity
# (beta), one point per state. A face-validity check on the fitted model: if
# beta and the remaining susceptible pool are correlated in the expected
# direction (more historical transmission -> more of the pool already
# consumed), that is evidence the finite-history fit recovers a mechanistically
# sensible relationship, not an artefact. This is a plausibility check, not a
# claim that immunity alone predicts BRR (state-to-state BRR variation
# reflects the joint, partly-offsetting effect of beta and susceptible pool
# together, which an 11-state sample cannot cleanly decompose into separate
# marginal effects -- see the setting-tier and dose-response figures above).
# =============================================================================
compute_state_beta_summary <- function(region) {
  p <- posterior_list[[region]]
  draw_mean_beta <- rowMeans(p$beta_observed, na.rm = TRUE)  # one number per draw
  tibble(
    Region = region,
    beta_med = median(draw_mean_beta, na.rm = TRUE),
    beta_lo  = quantile(draw_mean_beta, 0.025, na.rm = TRUE),
    beta_hi  = quantile(draw_mean_beta, 0.975, na.rm = TRUE)
  )
}

state_beta_summary <- purrr::map_dfr(names(posterior_list), compute_state_beta_summary)

immunity_beta_df <- state_immunity_summary %>%
  filter(AgeCat == "18-64") %>%
  mutate(Region = as.character(Region)) %>%
  inner_join(state_beta_summary, by = "Region")

cor_test <- cor.test(immunity_beta_df$susceptible_med, immunity_beta_df$beta_med, method = "spearman")
cor_lab <- sprintf("Spearman ρ = %.2f (p = %.3f, n = %d states)",
                    cor_test$estimate, cor_test$p.value, nrow(immunity_beta_df))

p_immunity_beta <- ggplot(immunity_beta_df, aes(susceptible_med, beta_med)) +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey80", alpha = 0.3, linewidth = 0.6) +
  geom_errorbar(aes(ymin = beta_lo, ymax = beta_hi), width = 0, colour = "grey50", alpha = 0.5) +
  geom_errorbarh(aes(xmin = susceptible_lo, xmax = susceptible_hi), height = 0, colour = "grey50", alpha = 0.5) +
  geom_point(aes(colour = setting), size = 2.6) +
  ggrepel::geom_text_repel(aes(label = Region), size = 2.8, colour = "grey20",
                            segment.size = 0.25, segment.colour = "grey60",
                            min.segment.length = 0, max.overlaps = Inf,
                            box.padding = 0.4, seed = 1) +
  scale_colour_manual(values = c(Low = "#55A868", Moderate = "#DD8452", High = "#C44E52"),
                       name = "Setting") +
  scale_x_continuous(labels = scales::percent,
                      name = "Susceptible pool remaining, 18–64y (start of 2022 fitting window)") +
  scale_y_continuous(name = "Fitted transmission rate, β (posterior median, mean over 52-week fitting window)") +
  labs(title = "Baseline immunity correlates with fitted transmission intensity across states",
       subtitle = "States with more historical transmission (higher β) show a correspondingly smaller remaining susceptible pool",
       caption = cor_lab) +
  theme_nm() +
  theme(plot.subtitle = element_text(colour = "grey30", size = 9))

ggsave("06_Results/supp_baseline_immunity_vs_beta.png", p_immunity_beta,
       width = 8.5, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/supp_baseline_immunity_vs_beta.png")

# =============================================================================
# SUSCEPTIBLE POOL -> BRR DOSE-RESPONSE (state level) -- does BRR functionally
# track baseline immunity/susceptible pool, independent of the discrete
# AR-based setting bucket? One point per state x age group x mechanism (DALY,
# base risk); x = susceptible pool from state_immunity_summary, y = the
# state's median BRR. A fitted trend line shows the dose-response shape
# directly, instead of only plotting the already-computed BRR on a plane.
# =============================================================================
state_brr_summary <- draw_level_xy_serostatus %>%
  filter(RR_seropos == 0, outcome == "DALY", AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64")) %>%
  group_by(Region, AgeCat, VE_label) %>%
  summarise(brr_med = median(brr_base, na.rm = TRUE),
            brr_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
            brr_hi  = quantile(brr_base, 0.975, na.rm = TRUE),
            .groups = "drop")

dose_response_df <- state_brr_summary %>%
  inner_join(state_immunity_summary %>%
               select(Region, AgeCat, susceptible_med, susceptible_lo, susceptible_hi, setting),
             by = c("Region", "AgeCat"))

p_dose_response <- ggplot(dose_response_df, aes(susceptible_med, brr_med)) +
  geom_smooth(method = "lm", colour = "grey40", fill = "grey70", alpha = 0.25, linewidth = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = brr_lo, ymax = brr_hi), width = 0, alpha = 0.3) +
  geom_errorbarh(aes(xmin = susceptible_lo, xmax = susceptible_hi), height = 0, alpha = 0.3) +
  geom_point(aes(colour = setting), size = 2.8) +
  geom_text(aes(label = Region), size = 2.4, colour = "grey30", vjust = -0.9, check_overlap = TRUE) +
  facet_grid(AgeCat ~ VE_label) +
  scale_colour_manual(values = c(Low = "#55A868", Moderate = "#DD8452", High = "#C44E52"),
                       name = "AR-based setting\n(current classification)") +
  scale_x_continuous(labels = scales::percent,
                      name = "Susceptible pool remaining (start of 2022 fitting window)") +
  scale_y_log10(labels = log_num_labels, name = "DALY BRR, base risk (log scale)") +
  labs(title = "Does baseline immunity functionally predict the BRR? (DALY, state level)",
       caption = paste("One point per state; error bars = 95% UI across 1000 posterior draws.",
                        "Line = linear trend (log BRR ~ susceptible pool), shaded band = 95% CI of the trend.",
                        "Point colour (AR-based setting) does not track this axis cleanly -- see previous figure.",
                        sep = "\n")) +
  theme_nm()

ggsave("06_Results/susceptible_pool_vs_brr_DALY.png", p_dose_response,
       width = 9, height = 8, dpi = 300, bg = "white")
message("Saved: 06_Results/susceptible_pool_vs_brr_DALY.png")

# =============================================================================
# BRR SENSITIVITY TORNADO -- does the baseline-immunity assumption (old
# long-term-average FOI vs new finite-history/short-term FOI) actually move
# the headline BRR? This is the figure Reviewer #2 asked for directly:
# "refitting the model under alternative baseline immunity assumptions and
# showing downstream vaccination results are not sensitive to that choice."
#
# Requires 01_Data/draw_level_xy_serostatus_original.RData -- the SAME
# object, saved from a DATA_MODE <- "original" run of
# 03_brazil_all_draws_ori_v3.R (long-term-average-FOI S0), analogous to how
# draw_level_xy_serostatus.RData was saved for the finite/short-term run.
# Not run automatically below (the file doesn't exist yet) -- call
# build_brr_tornado() once it's been saved.
# =============================================================================
build_brr_tornado <- function(outcome_sel = "DALY", use_adj = FALSE,
                               old_file = "01_Data/draw_level_xy_serostatus_longterm.RData") {
  if (!file.exists(old_file)) {
    stop("Missing ", old_file, " -- save draw_level_xy_serostatus from a ",
         'DATA_MODE <- "original" run first: save(draw_level_xy_serostatus, file = "',
         old_file, '")')
  }
  brr_col <- if (use_adj) "brr_adj" else "brr_base"

  summarise_brr <- function(df) {
    df %>%
      filter(RR_seropos == 0, outcome == outcome_sel, AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64")) %>%
      group_by(setting, VE_label, AgeCat) %>%
      summarise(brr_med = median(.data[[brr_col]], na.rm = TRUE),
                brr_lo  = quantile(.data[[brr_col]], 0.025, na.rm = TRUE),
                brr_hi  = quantile(.data[[brr_col]], 0.975, na.rm = TRUE),
                .groups = "drop")
  }

  new_summ <- summarise_brr(draw_level_xy_serostatus) %>%
    rename(brr_med_new = brr_med, brr_lo_new = brr_lo, brr_hi_new = brr_hi)

  e <- new.env()
  load(old_file, envir = e)
  old_summ <- summarise_brr(e$draw_level_xy_serostatus) %>%
    rename(brr_med_old = brr_med, brr_lo_old = brr_lo, brr_hi_old = brr_hi)

  d <- new_summ %>%
    inner_join(old_summ, by = c("setting", "VE_label", "AgeCat")) %>%
    mutate(
      label = paste(setting, AgeCat, VE_label, sep = " | "),
      fold_change = brr_med_new / brr_med_old,
      crosses_1 = (brr_med_old > 1) != (brr_med_new > 1)  # did the conclusion flip?
    ) %>%
    arrange(fold_change) %>%
    mutate(label = factor(label, levels = label))

  ggplot(d) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_segment(aes(x = brr_med_old, xend = brr_med_new, y = label, yend = label,
                      colour = crosses_1), linewidth = 1) +
    geom_point(aes(brr_med_old, label), shape = 21, fill = "white", colour = "grey30", size = 2.2) +
    geom_point(aes(brr_med_new, label, colour = crosses_1), size = 2.6) +
    scale_colour_manual(values = c("FALSE" = "#4C72B0", "TRUE" = "#C44E52"),
                         labels = c("FALSE" = "Conclusion unchanged", "TRUE" = "Conclusion flips"),
                         name = NULL) +
    scale_x_log10(labels = log_num_labels,
                   name = paste(outcome_sel, "BRR (log scale) — hollow = old long-term FOI, filled = new finite-history FOI")) +
    labs(y = NULL,
         title = paste0("Does the baseline-immunity assumption change the ", outcome_sel, " BRR conclusion?"),
         caption = "Each row: median BRR under the old (long-term average FOI) vs new (finite-history, short-term FOI) baseline immunity assumption.") +
    theme_nm()
}

for (oc in c("DALY", "Death", "SAE")) {
  p_tornado <- build_brr_tornado(oc)
  out_file <- sprintf("06_Results/brr_sensitivity_tornado_%s.png", oc)
  ggsave(out_file, p_tornado, width = 9, height = 8, dpi = 300, bg = "white")
  message("Saved: ", out_file)
}

# Long-term-vs-finite-history comparison figures (state-level age-immunity
# overlay, decision-map comparison) removed per project decision to report
# finite-history results only, now that the finite-history fitting is final.
# Benefit-risk-plane comparison and serostatus-mix explainer (both compared
# finite-history against long-term) removed for the same reason as above.
# =============================================================================
# NET BENEFIT-OR-RISK BAR CHART -- same visual grammar as Cracknell Daniels
# et al. (Nature Health 2025) Fig 3 panels d-f: horizontal bars showing
# benefit MINUS risk (absolute difference per 10,000 vaccinated, not the
# ratio), blue = net benefit, red = net risk, one row per Setting x Age
# group. Uses finite-history results only (this project's primary results,
# not the assumption-comparison). Death is scoped to 65+ only.
# =============================================================================
build_net_benefit_data <- function(outcome_sel, use_adj = FALSE) {
  risk_col <- if (use_adj) "x_10k_adj" else "x_10k_base"
  draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel,
           AgeCat %in% c("18-64", "65+"),
           !(outcome_sel == "Death" & AgeCat == "18-64")) %>%
    mutate(net_10k = averted_10k - .data[[risk_col]]) %>%
    group_by(setting, AgeCat, VE_label) %>%
    summarise(net_med = median(net_10k, na.rm = TRUE),
              net_lo  = quantile(net_10k, 0.025, na.rm = TRUE),
              net_hi  = quantile(net_10k, 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      setting = factor(setting, levels = c("Low", "Moderate", "High")),
      row_label = factor(paste0(setting, " | ", AgeCat),
                          levels = rev(as.vector(outer(c("Low", "Moderate", "High"),
                                                        c("18-64", "65+"), paste, sep = " | ")))),
      sign = ifelse(net_med >= 0, "Net benefit", "Net risk")
    )
}

build_net_benefit_panel <- function(mechanism_label, use_adj) {
  d <- bind_rows(
    build_net_benefit_data("DALY", use_adj) %>% mutate(outcome = "DALY"),
    build_net_benefit_data("SAE",  use_adj) %>% mutate(outcome = "SAE"),
    build_net_benefit_data("Death", use_adj) %>% mutate(outcome = "Death")
  ) %>%
    filter(VE_label == mechanism_label) %>%
    mutate(outcome = factor(outcome, levels = c("DALY", "SAE", "Death"))) %>%
    filter(!is.na(row_label))

  ggplot(d, aes(row_label, net_med, fill = sign)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = net_lo, ymax = net_hi), width = 0.2, colour = "grey30", linewidth = 0.4) +
    coord_flip() +
    facet_wrap(~outcome, scales = "free_x", nrow = 1) +
    scale_fill_manual(values = c("Net benefit" = "#274690", "Net risk" = "#C44E52"), name = NULL) +
    scale_y_continuous(labels = log_num_labels) +
    labs(x = NULL, y = "Difference in incidence per 10,000 vaccinated (benefit − risk)") +
    theme_nm() +
    theme(legend.position = "bottom")
}

build_and_save_net_benefit_chart <- function(use_adj) {
  serostatus_label <- if (use_adj) "Serostatus-adjusted risk" else "Base risk"

  p_disease_only <- build_net_benefit_panel("Disease blocking only", use_adj)
  p_both         <- build_net_benefit_panel("Disease and infection blocking", use_adj)

  p_combined <- (p_disease_only / p_both) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Net benefit or risk: ", serostatus_label),
      tag_levels = list(c("Disease blocking only", "Disease and infection blocking")),
      caption = "Bar = median, error bar = 95% UI across posterior draws, pooled across states within each setting. Finite-history results. Death restricted to 65+.",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    ) &
    theme(legend.position = "bottom")

  out_file <- sprintf("06_Results/net_benefit_risk_%s.png",
                       gsub("[^A-Za-z]+", "_", serostatus_label))
  ggsave(out_file, p_combined, width = 11, height = 9, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p_combined)
}

for (adj in c(FALSE, TRUE)) build_and_save_net_benefit_chart(adj)

# =============================================================================
# SEROSTATUS-ADJUSTMENT IMPACT SUMMARY -- collapses the setting x age
# granularity of the base-vs-adjusted comparison into one compact figure.
# averted_10k (benefit) is identical between base/adj (verified elsewhere in
# this script); serostatus adjustment only moves the risk side, so this is
# just % change in vaccine-attributable risk, adjusted vs base, pooled
# across states/settings/age groups, by outcome and mechanism.
# =============================================================================
build_serostatus_diff_data <- function(outcome_sel) {
  draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel,
           AgeCat %in% c("18-64", "65+"),
           !(outcome_sel == "Death" & AgeCat == "18-64")) %>%
    mutate(pct_diff = (x_10k_adj - x_10k_base) / x_10k_base * 100) %>%
    filter(is.finite(pct_diff)) %>%
    group_by(VE_label) %>%
    summarise(med = median(pct_diff, na.rm = TRUE),
              lo  = quantile(pct_diff, 0.025, na.rm = TRUE),
              hi  = quantile(pct_diff, 0.975, na.rm = TRUE), .groups = "drop") %>%
    mutate(outcome = outcome_sel)
}

build_and_save_serostatus_diff_chart <- function() {
  d <- bind_rows(
    build_serostatus_diff_data("DALY"),
    build_serostatus_diff_data("SAE"),
    build_serostatus_diff_data("Death")
  ) %>%
    mutate(outcome = factor(outcome, levels = c("DALY", "SAE", "Death")))

  p <- ggplot(d, aes(outcome, med, colour = VE_label)) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_pointrange(aes(ymin = lo, ymax = hi),
                     position = position_dodge(width = 0.4), size = 0.6, linewidth = 0.8) +
    scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
    labs(x = NULL, y = "Change in vaccine-attributable risk,\nserostatus-adjusted vs base (%)",
         title = "Impact of serostatus adjustment on attributable risk",
         caption = "Point = median, error bar = 95% UI across posterior draws, pooled across states, settings, and age groups. Finite-history results.") +
    theme_nm()

  out_file <- "06_Results/serostatus_adjustment_impact_summary.png"
  ggsave(out_file, p, width = 7, height = 5, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p)
}

build_and_save_serostatus_diff_chart()

# =============================================================================
# BRR ACCEPTABILITY CURVE (CEAC-style) -- Pr(BRR>1) as a single dot per
# setting/age (the old prob_dotplot) put almost every point at 99-100%,
# since BRR is comfortably >1 nearly everywhere here -- not informative,
# nothing to compare. Sweeping the threshold instead (Pr(BRR > t) as a
# function of t, not just t=1) shows how fast that confidence erodes as the
# bar is raised, which is where settings/ages actually separate. Reuses the
# make_brr_ceac_outbreak()/plot_brr_ceac_outbreak_ve() logic already
# validated in 03_brazil_all_draws_ori_v3.R, adapted to this file's
# finite-only data and Death-excludes-18-64 convention.
# =============================================================================
make_brr_ceac <- function(df, thresholds, group_vars) {
  df %>%
    tidyr::crossing(threshold = thresholds) %>%
    group_by(across(all_of(group_vars)), threshold) %>%
    summarise(p_accept = mean(brr > threshold), .groups = "drop")
}

build_brr_ceac_curve <- function(outcome_sel) {
  brr_long <- draw_level_xy_serostatus %>%
    filter(RR_seropos == 0, outcome == outcome_sel,
           AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64"),
           is.finite(brr_base), brr_base > 0, is.finite(brr_adj), brr_adj > 0) %>%
    mutate(setting = factor(setting, levels = c("Low", "Moderate", "High")),
           AgeCat  = factor(AgeCat,  levels = c("18-64", "65+"))) %>%
    tidyr::pivot_longer(cols = c(brr_base, brr_adj), names_to = "brr_type", values_to = "brr") %>%
    mutate(brr_type = factor(brr_type, levels = c("brr_base", "brr_adj"),
                              labels = c("Base", "Serostatus-adjusted")))

  thresholds <- 10^seq(floor(log10(min(brr_long$brr))), ceiling(log10(max(brr_long$brr))), by = 0.02)

  d <- make_brr_ceac(brr_long, thresholds,
                      group_vars = c("setting", "VE_label", "AgeCat", "brr_type"))

  p <- ggplot(d, aes(threshold, p_accept, colour = AgeCat, linetype = brr_type)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", alpha = 0.6) +
    geom_line(linewidth = 0.9) +
    facet_grid(setting ~ VE_label) +
    scale_x_log10(labels = log_num_labels, name = paste(outcome_sel, "BRR threshold (log scale)")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1), name = "Probability BRR exceeds threshold") +
    scale_colour_manual(values = c("18-64" = "#4C72B0", "65+" = "#DD8452"), name = "Age group") +
    scale_linetype_manual(values = c("Base" = "solid", "Serostatus-adjusted" = "dashed"), name = "Risk") +
    labs(title = paste0(outcome_sel, " — BRR acceptability curve"),
         caption = paste0("Dashed vertical line = BRR 1. Pr(BRR>threshold) across 1000 posterior draws, pooled across states within each setting.",
                           if (outcome_sel == "Death") " Restricted to 65+." else "")) +
    theme_nm()

  out_file <- sprintf("06_Results/brr_ceac_%s.png", outcome_sel)
  ggsave(out_file, p, width = 10, height = if (outcome_sel == "Death") 4.5 else 8, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p)
}

for (oc in c("DALY", "Death", "SAE")) build_brr_ceac_curve(oc)

# =============================================================================
# MECHANISM/SEROSTATUS BRR GRID -- same 3-column, one-representative-
# state-per-setting layout as build_mechanism_explainer(), but replacing the
# dense scatter-plane comparison with plain point-range line graphs (age
# group on x, BRR on y, colour = mechanism) -- one row per {base,
# serostatus-adjusted}, finite-history results only.
# =============================================================================
build_mechanism_brr_grid <- function(regions = c("Bahia", "Pernambuco", "Ceará"),
                                      setting_labels = c("Low", "Moderate", "High"),
                                      outcome_sel = "DALY") {

  row_specs <- list(
    list(label = "Base risk",               df = draw_level_xy_serostatus, use_adj = FALSE),
    list(label = "Serostatus-adjusted",      df = draw_level_xy_serostatus, use_adj = TRUE)
  )

  # Row identity is embedded directly in the leftmost panel's y-axis title
  # (bold, own line) rather than via patchwork's tag_levels -- tag_levels
  # tags every LEAF panel once nested wrap_plots() objects are stacked
  # (12 panels here, not the 4 row-groups intended), which produced
  # overlapping/misplaced labels.
  build_one_panel <- function(df, region, use_adj, title = NULL, row_label = NULL) {
    brr_col <- if (use_adj) "brr_adj" else "brr_base"
    d <- df %>%
      filter(Region == region, RR_seropos == 0, AgeCat %in% c("18-64", "65+"), !(outcome == "Death" & AgeCat == "18-64"), outcome == outcome_sel) %>%
      group_by(AgeCat, VE_label) %>%
      summarise(med = median(.data[[brr_col]], na.rm = TRUE),
                lo  = quantile(.data[[brr_col]], 0.025, na.rm = TRUE),
                hi  = quantile(.data[[brr_col]], 0.975, na.rm = TRUE), .groups = "drop")

    y_name <- paste(outcome_sel, "BRR (log)")
    if (!is.null(row_label)) y_name <- paste0(row_label, "\n", y_name)

    p <- ggplot(d, aes(AgeCat, med, colour = VE_label)) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.3),
                       fatten = 2.4, linewidth = 0.9) +
      scale_colour_manual(values = pal_mechanism, name = "Mechanism") +
      scale_y_log10(name = y_name, labels = log_num_labels) +
      labs(x = "Age group") +
      theme_nm(base_size = 13) +
      theme(axis.title.y = element_text(face = "bold", size = 11))
    if (!is.null(title)) p <- p + ggtitle(title)
    p
  }

  row_plots <- purrr::imap(row_specs, function(spec, row_i) {
    plots <- purrr::map2(regions, setting_labels, function(region, setting_lab) {
      ttl <- if (row_i == 1) paste0(region, " (", setting_lab, ")") else NULL
      is_leftmost <- region == regions[1]
      build_one_panel(spec$df, region, spec$use_adj, title = ttl,
                       row_label = if (is_leftmost) spec$label else NULL)
    })
    wrap_plots(plots, nrow = 1)
  })

  p_combined <- wrap_plots(row_plots, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(outcome_sel, " BRR by age group: base vs serostatus-adjusted risk, one state per setting"),
      caption = "Point = median, line = 95% UI across posterior draws. Dashed line = BRR 1 (log scale). Row identity labelled on each row's left-hand y-axis title. Finite-history results.",
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    ) &
    theme(legend.position = "bottom", legend.text = element_text(size = 11))

  out_file <- sprintf("06_Results/mechanism_brr_grid_%s.png", outcome_sel)
  ggsave(out_file, p_combined, width = 15, height = 9, dpi = 300, bg = "white")
  message("Saved: ", out_file)
  invisible(p_combined)
}

build_mechanism_brr_grid(outcome_sel = "DALY")
