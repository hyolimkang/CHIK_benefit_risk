### -----------------------------------------------------------------------
### Uncertainty companion to finite_weeksweep_heatmap.R: same x/y/facet
### layout (week x coverage, AgeCat x setting), but fill = Pr(BRR > 1)
### across posterior draws at each cell, instead of median BRR.
###
### The median-BRR heatmap shows WHERE the surface favours vaccination; this
### one shows HOW CONFIDENT that is. A cell can have a favourable median BRR
### while still being a coin flip (Pr(BRR>1) near 50%) if its 95% UI is wide
### -- that distinction is invisible in the median-only figure. Read the two
### heatmaps side by side, not this one alone (see chat record 2026-08-18).
###
### Same GAM-smoothing-of-a-sparse-grid caveat as finite_weeksweep_heatmap.R:
### visual aid for locating the 50% boundary, not new evidence between
### the 7 week x 4 coverage grid points.
###
### NOTE: this repo's copy is a mirror kept in sync with the CHIK_ORV_impact
### copy for reference/version control -- actual execution happens from
### CHIK_ORV_impact (00_Data/0_2_Processed relative paths), same convention
### as finite_weeksweep_heatmap.R.
### -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(scales)
library(purrr)

# ---- Config -----------------------------------------------------------
OUTCOMES        <- c("DALY", "SAE", "Death")
PROB_VERSIONS   <- c(base = "brr_base_prob_gt1", adj = "brr_adj_prob_gt1")
PROB_LABELS     <- c(base = "Pr(BRR>1)\n(unadjusted)", adj = "Pr(BRR>1)\n(seroadjusted)")
SCENARIOS_KEEP  <- c(3, 4)   # 18-64 and 65+
GRID_RES        <- 120

load("00_Data/0_2_Processed/brr_weeksweep_summary_finite.RData")      # brr_weeksweep_summary (VE98.9)
brr_weeksweep_summary_ve989 <- brr_weeksweep_summary
load("00_Data/0_2_Processed/brr_weeksweep_summary_finite_ve0.RData")  # brr_weeksweep_summary (VE0)
brr_weeksweep_summary <- dplyr::bind_rows(brr_weeksweep_summary_ve989, brr_weeksweep_summary)

VE_LABELS_KEEP <- c("Disease blocking only", "Disease and infection blocking")
VE_FILE_TAG    <- c("Disease blocking only" = "ve0", "Disease and infection blocking" = "ve989")

scenario_labels <- c(`3` = "18-64 years", `4` = "65+ years")

ink_primary   <- "#1a1a1a"
ink_secondary <- "#5a5a5a"
panel_border  <- "#c8c8c8"

smooth_one_facet <- function(df, week_range, cov_range) {
  grid <- expand.grid(
    week         = seq(week_range[1], week_range[2], length.out = GRID_RES),
    coverage_pct = seq(cov_range[1],  cov_range[2],  length.out = GRID_RES)
  )

  n_week <- length(unique(df$week))
  n_cov  <- length(unique(df$coverage_pct))
  k_week <- max(3, min(4, n_week - 1))
  k_cov  <- max(3, min(3, n_cov  - 1))

  fit <- tryCatch(
    mgcv::gam(prob_gt1 ~ te(week, coverage_pct, k = c(k_week, k_cov)), data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      mgcv::gam(prob_gt1 ~ s(week, k = k_week) + s(coverage_pct, k = k_cov), data = df),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) {
    warning("Facet with ", nrow(df), " points has too few points to smooth -- returning NA surface.")
    grid$prob_gt1 <- NA_real_
    return(grid)
  }

  grid$prob_gt1 <- pmin(pmax(predict(fit, newdata = grid), 0), 1)
  grid
}

make_prob_heatmap <- function(outcome, version, ve_label) {

  prob_column <- PROB_VERSIONS[[version]]

  # Death: vaccine-attributable death risk is fixed at exactly 0 for 18-64
  # (01_setup.R's p_death_vacc_u65 <- 0), so brr_death is structurally
  # undefined there -- keep 65+ only, for BOTH mechanisms consistently (chat
  # record 2026-09-02; the raw brr_weeksweep_summary data is inconsistent
  # about whether 18-64/Death survives the is.finite() filter below across
  # DB vs D+I, which would otherwise silently render a 1-row DB panel next
  # to a 2-row D+I panel -- forcing 65+-only here avoids that asymmetry).
  scenarios_keep_outcome <- if (outcome == "Death") 4 else SCENARIOS_KEEP

  plot_data <- brr_weeksweep_summary %>%
    filter(
      outcome == !!outcome,
      RR_seropos == 0,
      VE_label == !!ve_label,
      Scenario %in% scenarios_keep_outcome
    ) %>%
    mutate(
      coverage_pct = as.numeric(gsub("cov", "", Coverage)),
      setting      = factor(setting, levels = c("Low", "Moderate", "High")),
      age_label    = factor(scenario_labels[as.character(Scenario)],
                            levels = scenario_labels[as.character(scenarios_keep_outcome)]),
      prob_gt1     = .data[[prob_column]]
    ) %>%
    filter(is.finite(prob_gt1)) %>%
    select(setting, age_label, week, coverage_pct, prob_gt1)

  stopifnot(
    "No rows survived filtering -- check outcome/SCENARIOS_KEEP/prob_column against brr_weeksweep_summary's actual values" =
      nrow(plot_data) > 0
  )

  week_range <- range(plot_data$week)
  cov_range  <- range(plot_data$coverage_pct)

  smoothed <- plot_data %>%
    group_by(setting, age_label) %>%
    group_modify(~ smooth_one_facet(.x, week_range, cov_range)) %>%
    ungroup()

  # Diverging, continuous, centred at 50% (as-likely-as-not) -- same RdBu
  # convention as the median-BRR heatmap: red = more likely risk > benefit,
  # blue = more likely benefit > risk.
  ramp_colours <- RColorBrewer::brewer.pal(11, "RdBu")

  base_case_point <- data.frame(week = 2, coverage_pct = 50)

  p_prob_heatmap <- ggplot(smoothed, aes(x = week, y = coverage_pct, fill = prob_gt1)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = prob_gt1), breaks = 0.5, colour = "#1a1a1a", linewidth = 0.5, linetype = "22") +
    geom_point(
      data = base_case_point, aes(x = week, y = coverage_pct),
      inherit.aes = FALSE, shape = 23, size = 1.5, stroke = 0.35,
      fill = "#f5c518", colour = "#1a1a1a"
    ) +
    facet_grid(age_label ~ setting) +
    scale_fill_gradientn(
      colours = ramp_colours,
      breaks  = c(0, 0.25, 0.5, 0.75, 1),
      limits  = c(0, 1),
      labels  = scales::percent_format(accuracy = 1),
      name    = PROB_LABELS[[version]],
      guide   = guide_colorbar(barheight = unit(110, "pt"), barwidth = unit(7, "pt"))
    ) +
    scale_x_continuous(breaks = c(1, 13, 26, 39, 52), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = c(10, 30, 60, 90),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = week_range, ylim = cov_range) +
    labs(
      title = sprintf("%s  |  %s", ve_label, if (version == "adj") "Seroadjusted risk" else "Base risk (unadjusted for serostatus)"),
      subtitle = "◆ Reference scenario: week 2, 50% coverage",
      x = "Vaccination campaign start week",
      y = "Coverage"
    ) +
    theme_minimal(base_size = 15, base_family = "sans") +
    theme(
      text                  = element_text(colour = ink_primary),
      plot.title            = element_text(face = "bold", size = 16, colour = ink_primary, margin = margin(b = 2)),
      plot.subtitle         = element_text(size = 12.5, face = "bold", colour = "#8a6d00", margin = margin(b = 6)),
      plot.caption          = element_text(size = 10, colour = ink_secondary, hjust = 0,
                                            lineheight = 1.15, margin = margin(t = 8)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title            = element_text(size = 13, colour = ink_secondary),
      axis.text             = element_text(size = 11.5, colour = ink_primary),
      axis.ticks            = element_line(colour = panel_border, linewidth = 0.3),
      axis.ticks.length     = unit(3, "pt"),
      panel.grid            = element_blank(),
      panel.spacing         = unit(18, "pt"),
      panel.border          = element_rect(colour = panel_border, fill = NA, linewidth = 0.4),
      strip.text.x          = element_text(face = "bold", size = 12.5, colour = ink_primary, margin = margin(b = 5)),
      strip.text.y          = element_text(face = "bold", size = 12.5, colour = ink_primary, margin = margin(l = 5)),
      strip.background      = element_blank(),
      legend.title          = element_text(size = 12, colour = ink_secondary, lineheight = 1.05),
      legend.text           = element_text(size = 11, colour = ink_primary),
      legend.position        = "right",
      plot.margin            = margin(12, 14, 10, 12)
    )

  dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)

  file_stub <- sprintf("02_Outputs/2_1_Figures/brr_weeksweep_heatmap_probgt1_%s_%s_%s",
                        tolower(outcome), version, VE_FILE_TAG[[ve_label]])

  ggsave(paste0(file_stub, ".png"), plot = p_prob_heatmap, width = 7.2, height = 5.2, dpi = 600, bg = "white")

  message("Saved: ", file_stub, ".{pdf,png}")
  p_prob_heatmap
}

plot_grid_spec <- expand.grid(
  outcome = OUTCOMES, version = names(PROB_VERSIONS), ve_label = VE_LABELS_KEEP,
  stringsAsFactors = FALSE
)
plots <- purrr::pmap(plot_grid_spec, make_prob_heatmap)
names(plots) <- paste(plot_grid_spec$outcome, plot_grid_spec$version, VE_FILE_TAG[plot_grid_spec$ve_label], sep = "_")

print(plots[["DALY_adj_ve989"]])

# ---- Combined-mechanism figure: Disease blocking only stacked with Disease
# and infection blocking in ONE file, so the two can be read side by side
# without flipping between separate PNGs. Same stacking convention as
# net_benefit_risk_{Base,Serostatus_adjusted}_risk.png in
# 10_streamline_mechanism_figure.R (mechanism panels stacked via patchwork,
# shared legend collected). Both panels share the identical Pr(BRR>1) colour
# scale (same breaks/limits/name), so plot_layout(guides = "collect") merges
# them into a single legend rather than showing it twice.
library(patchwork)

make_prob_heatmap_combined <- function(outcome, version) {
  # Strip each panel's own title/caption (they duplicated the outer title
  # 3x and the caption was getting clipped at the figure edge) -- keep just
  # a short mechanism-name title per panel; no caption on the combined figure.
  p_ve0   <- plots[[paste(outcome, version, "ve0",   sep = "_")]] +
    labs(title = "Disease blocking only", caption = NULL)
  p_ve989 <- plots[[paste(outcome, version, "ve989", sep = "_")]] +
    labs(title = "Disease and infection blocking", caption = NULL)

  # DALY (Panel B, main manuscript) stays title-less, per chat record
  # 2026-09-02. SAE/Death (Panels C/D, supplementary) get their outcome name
  # as a short title -- "Disease blocking only"/"Disease and infection
  # blocking" sub-panel titles (set above) still identify each row within it.
  overall_title <- if (outcome %in% c("SAE", "Death")) outcome else NULL

  p_combined <- (p_ve0 / p_ve989) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = overall_title,
      # Extra top margin opens blank space above the title/top row title --
      # the outer corner tag (added below) sits there instead of colliding
      # with it (both anchor to the same top-left corner by default). Wider
      # than the title-less DALY case since a title is present here too --
      # tag and title otherwise sit almost flush against each other.
      theme = theme(
        plot.title  = element_text(face = "bold", size = 17, colour = ink_primary,
                                    margin = margin(t = if (!is.null(overall_title)) 12 else 0)),
        plot.margin = margin(t = if (!is.null(overall_title)) 34 else 24, r = 5.5, b = 5.5, l = 5.5)
      )
    ) &
    theme(legend.position = "right")

  # Whole-figure corner tag: "B" for DALY (Panel B of the main manuscript
  # figure; the headline benefit-risk plane is Panel A -- 10_streamline_
  # mechanism_figure.R). SAE/Death are supplementary companions, tagged "C"/
  # "D" to match their own headline_benefit_risk_plane_{SAE,Death} pairing
  # ("A. SAE benefit-risk space" / "B. Death benefit-risk space" -- chat
  # record 2026-09-02; note that lettering is local to those two files, so
  # "C"/"D" here don't collide with it). wrap_elements() turns the
  # already-composed DB/D+I 2-row plot into a single patch so
  # plot_annotation(tag_levels) labels the WHOLE figure once, not each
  # mechanism row again.
  tag_letter <- switch(outcome, DALY = "B", SAE = "C", Death = "D")
  p_tagged <- patchwork::wrap_elements(full = p_combined) +
    plot_annotation(tag_levels = list(tag_letter)) &
    theme(plot.tag = element_text(face = "bold", size = 21, colour = ink_primary))

  dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)
  out_file <- sprintf("02_Outputs/2_1_Figures/brr_weeksweep_heatmap_probgt1_%s_%s_mechanism_combined.png",
                       tolower(outcome), version)
  ggsave(out_file, p_tagged, width = 9.6, height = 9.6, dpi = 600, bg = "white")
  message("Saved: ", out_file)
  invisible(p_tagged)
}

combined_spec <- expand.grid(outcome = OUTCOMES, version = names(PROB_VERSIONS), stringsAsFactors = FALSE)
purrr::pwalk(combined_spec, make_prob_heatmap_combined)
