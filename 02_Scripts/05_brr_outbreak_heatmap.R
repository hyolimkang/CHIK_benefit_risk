# =============================================================================
# 05_brr_outbreak_heatmap.R
#
# Pr(BRR > 1) + median BRR tile figure (outbreak response immunisation style)
# using objects created by 03_brazil_all_draws_ori_v3.R
#
# Prerequisites (same R session after v3 finishes, or load saved objects):
#   draw_level_xy_serostatus   — draw-level benefit/risk (required)
#
# Usage:
#   source("02_Scripts/03_brazil_all_draws_ori_v3.R")   # long run
#   source("02_Scripts/05_brr_outbreak_heatmap.R")
#
# Or, if you saved draw-level data earlier:
#   load("01_Data/draw_level_xy_serostatus.RData")
# =============================================================================

if (!exists("draw_level_xy_serostatus", inherits = TRUE)) {
  stop(
    "Object 'draw_level_xy_serostatus' not found. ",
    "Run source('02_Scripts/03_brazil_all_draws_ori_v3.R') first ",
    "or load('01_Data/draw_level_xy_serostatus.RData').",
    call. = FALSE
  )
}

# ---- User options ------------------------------------------------------------
OUTCOME   <- "SAE"          # "SAE", "Death", or "DALY"
COVERAGE  <- "cov50"
AGE_LEVELS <- c("18-64", "65+")

# Brazil v3 uses epidemic settings Low / Moderate / High (from setting_key).
# Relabel rows here to match your figure legend if needed.
SETTING_LEVELS <- c("Low", "Moderate", "High")
SETTING_LABELS <- c(
  Low      = "Low transmission",
  Moderate = "Moderate post-outbreak",
  High     = "High post-outbreak"
)

SAVE_PATH <- file.path("06_Results", paste0("brr_outbreak_heatmap_", OUTCOME, ".png"))
PLOT_WIDTH  <- 11
PLOT_HEIGHT <- 7

# ---- Build long-format heatmap data ------------------------------------------
prep_brr_heatmap_data <- function(
    draw_df,
    outcome = "SAE",
    coverage = "cov50",
    age_levels = c("18-64", "65+"),
    setting_levels = c("Low", "Moderate", "High"),
    setting_labels = SETTING_LABELS
) {
  scen_for_age <- c("18-64" = 3L, "65+" = 4L)

  base_df <- draw_df %>%
    dplyr::filter(
      outcome == !!outcome,
      Coverage == !!coverage,
      RR_seropos == 0,
      AgeCat %in% age_levels,
      setting %in% setting_levels
    ) %>%
    dplyr::mutate(
      Scenario = as.integer(Scenario),
      setting  = factor(setting, levels = setting_levels),
      AgeCat   = factor(AgeCat, levels = age_levels),
      VE_label = factor(
        VE,
        levels = c("VE0", "VE98.9"),
        labels = c("Disease blocking only", "Disease and infection blocking")
      )
    ) %>%
    dplyr::filter(Scenario == scen_for_age[as.character(AgeCat)])

  prob_df <- base_df %>%
    dplyr::group_by(setting, AgeCat, VE_label) %>%
    dplyr::summarise(
      prob_base = mean(is.finite(brr_base) & brr_base > 1, na.rm = TRUE),
      prob_adj  = mean(is.finite(brr_adj)  & brr_adj  > 1, na.rm = TRUE),
      .groups = "drop"
    )

  med_df <- base_df %>%
    dplyr::group_by(setting, AgeCat, VE_label) %>%
    dplyr::summarise(
      brr_base_med = stats::median(brr_base[is.finite(brr_base) & brr_base > 0],
                                   na.rm = TRUE),
      brr_adj_med  = stats::median(brr_adj[is.finite(brr_adj) & brr_adj > 0],
                                   na.rm = TRUE),
      .groups = "drop"
    )

  row_levels <- unlist(
    lapply(setting_levels, function(s) {
      paste0(setting_labels[[s]], "\n", age_levels)
    }),
    use.names = FALSE
  )

  prob_df %>%
    dplyr::left_join(med_df, by = c("setting", "AgeCat", "VE_label")) %>%
    tidyr::pivot_longer(
      cols = c(prob_base, prob_adj),
      names_to = "risk_type",
      values_to = "prob"
    ) %>%
    dplyr::mutate(
      risk_type = dplyr::recode(
        risk_type,
        prob_base = "Base risk",
        prob_adj  = "Adjusted risk"
      ),
      brr_med = dplyr::if_else(
        risk_type == "Base risk", brr_base_med, brr_adj_med
      ),
      setting_label = unname(setting_labels[as.character(setting)]),
      row_key = factor(
        paste0(setting_label, "\n", AgeCat),
        levels = row_levels
      ),
      col_key = factor(
        paste0(VE_label, "\n", risk_type),
        levels = c(
          "Disease blocking only\nBase risk",
          "Disease blocking only\nAdjusted risk",
          "Disease and infection blocking\nBase risk",
          "Disease and infection blocking\nAdjusted risk"
        )
      ),
      prob_pct = 100 * prob,
      label = dplyr::if_else(
        is.finite(brr_med),
        sprintf("%0.0f%%\n(%s)", 100 * prob, formatC(brr_med, digits = 2, format = "f")),
        sprintf("%0.0f%%\n(—)", 100 * prob)
      )
    )
}

# ---- Plot --------------------------------------------------------------------
plot_brr_outbreak_heatmap <- function(
    heat_df,
    title = "Outbreak response immunisation",
    subtitle = NULL,
    fill_mid = 50
) {
  heat_df <- heat_df %>%
    dplyr::filter(!is.na(row_key), !is.na(col_key))

  if (!nrow(heat_df)) {
    stop("No rows to plot after filtering. Check OUTCOME / settings / coverage.",
         call. = FALSE)
  }

  ggplot2::ggplot(
    heat_df,
    ggplot2::aes(x = col_key, y = row_key, fill = prob_pct)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 1.1) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 3.2,
      lineheight = 0.95,
      color = "grey15"
    ) +
    ggplot2::scale_fill_gradient2(
      name = "Pr(BRR > 1)",
      low = "#f4a582",
      mid = "#ffffbf",
      high = "#008080",
      midpoint = fill_mid,
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      guide = ggplot2::guide_colorbar(
        barwidth = ggplot2::unit(0.35, "npc"),
        barheight = ggplot2::unit(0.55, "npc")
      )
    ) +
    ggplot2::scale_x_discrete(
      position = "top",
      labels = function(x) {
        vapply(strsplit(x, "\n"), function(z) {
          paste0(z[1], "\n", z[2])
        }, character(1))
      }
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      caption = "Tile colour = Pr(BRR > 1); label = Pr(BRR > 1) (median BRR)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey30"),
      plot.caption = ggplot2::element_text(size = 8.5, colour = "grey40"),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9, lineheight = 0.95),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.95),
      legend.position = "right"
    )
}

# ---- Run ---------------------------------------------------------------------
heatmap_data <- prep_brr_heatmap_data(
  draw_df          = draw_level_xy_serostatus,
  outcome          = OUTCOME,
  coverage         = COVERAGE,
  age_levels       = AGE_LEVELS,
  setting_levels   = SETTING_LEVELS,
  setting_labels   = SETTING_LABELS
)

p_brr_heatmap <- plot_brr_outbreak_heatmap(
  heat_df   = heatmap_data,
  title     = "Outbreak response immunisation",
  subtitle  = paste0(OUTCOME, " · ", COVERAGE, " · RR_seropos = 0")
)

print(p_brr_heatmap)

if (!dir.exists("06_Results")) dir.create("06_Results", recursive = TRUE)
ggplot2::ggsave(
  SAVE_PATH,
  plot = p_brr_heatmap,
  width = PLOT_WIDTH,
  height = PLOT_HEIGHT,
  dpi = 300,
  bg = "white"
)
message("[05_heatmap] Saved: ", SAVE_PATH)
