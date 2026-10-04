### -----------------------------------------------------------------------
### Publication-quality BRR heatmap: x = vaccination start week,
### y = coverage, fill = BRR (smoothed), faceted by AgeCat x setting --
### mirrors the style of the reference travel-BRR contour figure (blue/white/
### red diverging, dashed BRR=1 contour), rendered from
### brr_weeksweep_summary_finite.RData (finite_weeksweep_brr.R output).
###
### The raw grid is sparse (7 weeks x 4 coverage = 28 points per facet), so a
### 2D tensor-product GAM smooth is fit per facet and predicted onto a fine
### grid for a continuous-looking surface -- the same reason the reference
### figure's underlying grid is almost certainly coarser than its rendered
### resolution. This does NOT re-run any simulation; it only interpolates
### between already-computed BRR values, so read the smoothed surface as a
### visual aid for locating the BRR=1 boundary, not as new evidence between
### grid points.
###
### Colour is STEPPED (binned), not a continuous gradient: a continuous
### diverging fill fades to near-white across a wide band around BRR = 1,
### which washes out exactly the region readers most need to compare against
### the BRR = 1 contour. Discrete bins keep every band visibly distinct.
### -----------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)
library(scales)
library(purrr)

# ---- Config -----------------------------------------------------------
OUTCOMES       <- c("DALY", "SAE", "Death")
# Both the unadjusted ("base") and serostatus-adjusted ("adj", RR_seropos = 0
# -- seropositive vaccinees assumed to carry zero vaccine-attributable risk)
# BRR are rendered, one heatmap file per outcome x version.
BRR_VERSIONS   <- c(base = "brr_base_med", adj = "brr_adj_med")
BRR_LABELS     <- c(base = "BRR\n(unadjusted)", adj = "BRR\n(seroadjusted)")
SCENARIOS_KEEP <- c(3, 4)                         # 18-64 and 65+ -- the two adult/policy scenarios
GRID_RES       <- 120                              # smoothed prediction grid resolution per axis
# BRR is a ratio, so it's colour-mapped on log10(BRR) with the diverging
# midpoint at BRR = 1. But BRR can also go NEGATIVE (averted burden < 0 --
# vaccination increases burden on the median draw, not just "risk exceeds
# benefit" but "there is no benefit at all"), which log10() can't represent.
# Those cells -- and any BRR beyond this cap on the positive side -- are
# clipped to [BRR_FLOOR, BRR_CEIL] for display purposes only (all clipping
# happens after the med/lo/hi draw-level summary already computed upstream,
# so this never touches the underlying estimate, only how far the colour
# scale's tail extends). A BRR at or below BRR_FLOOR reads as "no net
# benefit," which is the correct qualitative message for both a tiny
# positive BRR and a negative one.
BRR_FLOOR      <- 0.01
BRR_CEIL       <- 100

# Both VE mechanisms: VE98.9 = "Disease and infection blocking" (the
# original sweep run), VE0 = "Disease blocking only" (added later so the
# headline heatmap isn't a mechanism-selective half-picture -- see chat
# record 2026-08-17). Load both and stack; VE_label already distinguishes
# them downstream (set in finite_weeksweep_brr.R from the VE column).
load("00_Data/0_2_Processed/brr_weeksweep_summary_finite.RData")      # brr_weeksweep_summary (VE98.9)
brr_weeksweep_summary_ve989 <- brr_weeksweep_summary
load("00_Data/0_2_Processed/brr_weeksweep_summary_finite_ve0.RData")  # brr_weeksweep_summary (VE0) -- overwrites the name above
brr_weeksweep_summary <- dplyr::bind_rows(brr_weeksweep_summary_ve989, brr_weeksweep_summary)

VE_LABELS_KEEP <- c("Disease blocking only", "Disease and infection blocking")
VE_FILE_TAG    <- c("Disease blocking only" = "ve0", "Disease and infection blocking" = "ve989")

scenario_labels <- c(`3` = "18-64 years", `4` = "65+ years")

week_range_global <- range(brr_weeksweep_summary$week)

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
  # te() marginal bases need k strictly less than the number of unique
  # values on that margin; cap conservatively and fall back to a coarser
  # smooth (then to no smoothing at all) if a facet has too few finite BRR
  # points to support even a minimal tensor smooth.
  k_week <- max(3, min(4, n_week - 1))
  k_cov  <- max(3, min(3, n_cov  - 1))

  fit <- tryCatch(
    mgcv::gam(log10(brr) ~ te(week, coverage_pct, k = c(k_week, k_cov)), data = df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      mgcv::gam(log10(brr) ~ s(week, k = k_week) + s(coverage_pct, k = k_cov), data = df),
      error = function(e) NULL
    )
  }
  if (is.null(fit)) {
    warning("Facet with ", nrow(df), " points (", n_week, " weeks x ", n_cov,
            " coverage levels) has too few points to smooth -- returning NA surface.")
    grid$log10_brr <- NA_real_
    return(grid)
  }

  grid$log10_brr <- predict(fit, newdata = grid)
  grid
}

make_brr_heatmap <- function(outcome, version, ve_label) {

  brr_column <- BRR_VERSIONS[[version]]

  plot_data <- brr_weeksweep_summary %>%
    filter(
      outcome == !!outcome,
      RR_seropos == 0,
      VE_label == !!ve_label,
      Scenario %in% SCENARIOS_KEEP
    ) %>%
    mutate(
      coverage_pct = as.numeric(gsub("cov", "", Coverage)),
      setting      = factor(setting, levels = c("Low", "Moderate", "High")),
      age_label    = factor(scenario_labels[as.character(Scenario)],
                            levels = scenario_labels[as.character(SCENARIOS_KEEP)]),
      brr_raw      = .data[[brr_column]]
    ) %>%
    # Only drop genuinely undefined cells (NA/Inf/NaN, e.g. zero vaccine-
    # attributable risk in the denominator). Negative/near-zero BRR is kept
    # and clipped below, not dropped.
    filter(is.finite(brr_raw)) %>%
    mutate(brr = pmax(pmin(brr_raw, BRR_CEIL), BRR_FLOOR)) %>%
    select(setting, age_label, week, coverage_pct, brr, brr_raw)

  stopifnot(
    "No rows survived filtering -- check outcome/SCENARIOS_KEEP/brr_column against brr_weeksweep_summary's actual values" =
      nrow(plot_data) > 0
  )

  week_range <- range(plot_data$week)
  cov_range  <- range(plot_data$coverage_pct)

  smoothed <- plot_data %>%
    group_by(setting, age_label) %>%
    group_modify(~ smooth_one_facet(.x, week_range, cov_range)) %>%
    ungroup()

  # ---- Colour scale: diverging, CONTINUOUS, centred at BRR = 1 ------------
  # Continuous, not binned -- and only clean power-of-ten breaks
  # (0.01/0.1/1/10/100), no in-between values like 0.3/3/30.
  #
  # The fix for "1-100 all looks like the same blue" is NOT a transform
  # that compresses the scale (that's what caused it: the earlier
  # sign-preserving-sqrt-of-log10 remap pulls big |log10(BRR)| values
  # closer together, so it was squeezing exactly the far end of the range
  # where DALY's real values sit). Instead, use the full ~11-stop RdBu
  # ramp as a continuous gradientn (not just a 3-colour low/mid/high
  # gradient2) on a plain log10 scale -- every decade gets an equal,
  # undistorted share of the ramp, so distinct multiples of BRR actually
  # look like distinct colours across the whole 0.01-100 range.
  brr_break_vals <- c(0.01, 0.1, 1, 10, 100)
  # brewer.pal("RdBu") returns red (low index) -> blue (high index) already,
  # which is exactly what we want ascending against fill_val: low BRR
  # (risk > benefit) = red, high BRR (benefit > risk) = blue -- no rev() here.
  ramp_colours <- RColorBrewer::brewer.pal(11, "RdBu")
  smoothed <- smoothed %>% mutate(fill_val = pmax(pmin(10^log10_brr, BRR_CEIL), BRR_FLOOR))

  # Base-case marker: week 2 / 50% coverage, the fixed point the previous
  # headline_benefit_risk_plane_*.png figures were built at. Same data frame
  # (no facetting columns) drawn on every panel via inherit.aes = FALSE, so it
  # shows up in all 6 setting x age facets at the identical (week, coverage).
  base_case_point <- data.frame(week = 2, coverage_pct = 50)

  p_brr_heatmap <- ggplot(smoothed, aes(x = week, y = coverage_pct, fill = fill_val)) +
    geom_raster(interpolate = TRUE) +
    geom_contour(aes(z = log10_brr), breaks = 0, colour = "#1a1a1a", linewidth = 0.5, linetype = "22") +
    geom_point(
      data = base_case_point, aes(x = week, y = coverage_pct),
      inherit.aes = FALSE, shape = 23, size = 1.5, stroke = 0.35,
      fill = "#f5c518", colour = "#1a1a1a"
    ) +
    facet_grid(age_label ~ setting) +
    scale_fill_gradientn(
      colours   = ramp_colours,
      breaks    = brr_break_vals,
      trans     = "log10",
      limits    = c(BRR_FLOOR, BRR_CEIL),
      labels    = brr_break_vals,
      oob       = scales::squish,
      name      = BRR_LABELS[[version]],
      guide     = guide_colorbar(barheight = unit(110, "pt"), barwidth = unit(7, "pt"))
    ) +
    scale_x_continuous(breaks = c(1, 13, 26, 39, 52), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = c(10, 30, 60, 90),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = week_range, ylim = cov_range) +
    labs(
      # Mechanism (VE_label) and scenario (base/adj) as their own explicit,
      # bold first line -- not buried in the methodological footnote --
      # since this heatmap now covers both VE mechanisms and both risk
      # scenarios across separate files, and the mechanism/scenario is the
      # single most important thing distinguishing one file from another.
      title = sprintf("%s  |  %s", ve_label, if (version == "adj") "Seroadjusted risk (RR_seropos = 0)" else "Base risk (unadjusted for serostatus)"),
      caption = sprintf(
        "%s — dashed line = BRR 1; diamond = base-case scenario (week 2, 50%% coverage).\n%s; colour capped at BRR in [0.01, 100].",
        outcome,
        if (version == "adj") {
          "Seropositive vaccinees assumed zero vaccine-attributable risk beyond this scenario's baseline"
        } else {
          "Unadjusted for recipient serostatus"
        }
      ),
      x        = "Vaccination campaign start week",
      y        = "Coverage"
    ) +
    theme_minimal(base_size = 10.5, base_family = "sans") +
    theme(
      text                 = element_text(colour = ink_primary),
      plot.title           = element_text(face = "bold", size = 11.5, colour = ink_primary, margin = margin(b = 4)),
      plot.subtitle        = element_blank(),
      plot.caption         = element_text(size = 7.3, colour = ink_secondary, hjust = 0,
                                           lineheight = 1.15, margin = margin(t = 8)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      axis.title           = element_text(size = 9.5, colour = ink_secondary),
      axis.text            = element_text(size = 8, colour = ink_primary),
      axis.ticks           = element_line(colour = panel_border, linewidth = 0.3),
      axis.ticks.length    = unit(3, "pt"),
      panel.grid           = element_blank(),
      panel.spacing        = unit(9, "pt"),
      panel.border         = element_rect(colour = panel_border, fill = NA, linewidth = 0.4),
      strip.text.x         = element_text(face = "bold", size = 9, colour = ink_primary, margin = margin(b = 5)),
      strip.text.y         = element_text(face = "bold", size = 9, colour = ink_primary, margin = margin(l = 5)),
      strip.background     = element_blank(),
      legend.title         = element_text(size = 8.5, colour = ink_secondary, lineheight = 1.05),
      legend.text          = element_text(size = 7.8, colour = ink_primary),
      legend.position       = "right",
      plot.margin           = margin(12, 14, 10, 12)
    )

  dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)

  file_stub <- sprintf("02_Outputs/2_1_Figures/brr_weeksweep_heatmap_%s_%s", tolower(outcome), version)

  ggsave(paste0(file_stub, ".png"), plot = p_brr_heatmap, width = 7.2, height = 5.2, dpi = 600, bg = "white")

  message("Saved: ", file_stub, ".{pdf,png}")
  p_brr_heatmap
}

plot_grid_spec <- expand.grid(outcome = OUTCOMES, version = names(BRR_VERSIONS), stringsAsFactors = FALSE)
plots <- Map(make_brr_heatmap, plot_grid_spec$outcome, plot_grid_spec$version)
names(plots) <- paste(plot_grid_spec$outcome, plot_grid_spec$version, sep = "_")

print(plots[["DALY_adj"]])
