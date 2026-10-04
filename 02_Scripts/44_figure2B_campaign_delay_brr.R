# =============================================================================
# 44_figure2B_campaign_delay_brr.R
#
# Figure 2B: Benefit-risk ratio and outbreak response immunisation campaign
# delay.
#
#   x = campaign start week
#   y = median DALY BRR (log scale)
#   columns = Low / Moderate / High setting
#   rows = 18-64 / >=65 years
#   colour = protection mechanism
#   linetype = base versus serostatus-adjusted risk assumption
#
# Thick lines are interpolated to 50% coverage. Interpolation is only within
# the observed week x coverage grid; no extrapolation is used.
# =============================================================================

if (!file.exists("06_Results/brr_weeksweep_summary_finite.xlsx") ||
    !file.exists("06_Results/brr_weeksweep_summary_finite_ve0.xlsx")) {
  stop("ORI summary Excel files are missing from 06_Results/")
}

suppressMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

summary_data <- bind_rows(
  read_excel("06_Results/brr_weeksweep_summary_finite.xlsx"),
  read_excel("06_Results/brr_weeksweep_summary_finite_ve0.xlsx")
) %>%
  filter(outcome == "DALY", RR_seropos == 0,
         AgeCat %in% c("18-64", "65+")) %>%
  mutate(
    age_group = factor(AgeCat, levels = c("18-64", "65+"),
                       labels = c("18-64", "≥65")),
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    mechanism = factor(VE_label,
      levels = c("Disease blocking only", "Disease and infection blocking")),
    week = as.numeric(week),
    coverage = as.numeric(sub("cov", "", Coverage))
  )

# Interpolate over the observed coverage range at each campaign start week.
# The endpoints are the actual 10% and 90% coverage values; 50% is the fixed
# reference coverage for the thick central line.
interpolate_coverage <- function(data, value_col) {
  data %>%
    group_by(age_group, setting, mechanism, week) %>%
    group_modify(~ {
      x <- .x$coverage
      y <- .x[[value_col]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 2) {
        return(tibble(coverage = c(10, 50, 90), BRR = NA_real_))
      }
      o <- order(x[ok])
      x <- x[ok][o]
      y <- y[ok][o]
      tibble(
        coverage = c(10, 50, 90),
        BRR = approx(x, y, xout = c(10, 50, 90), rule = 1)$y
      )
    }) %>%
    ungroup()
}

base_grid <- interpolate_coverage(summary_data, "brr_base_med") %>%
  mutate(risk_assumption = "Base risk")
adj_grid <- interpolate_coverage(summary_data, "brr_adj_med") %>%
  mutate(risk_assumption = "Serostatus-adjusted risk")

coverage_grid <- bind_rows(base_grid, adj_grid) %>%
  mutate(
    risk_assumption = factor(risk_assumption,
      levels = c("Base risk", "Serostatus-adjusted risk")),
    BRR = ifelse(is.finite(BRR) & BRR > 0, BRR, NA_real_)
  )

line_data <- coverage_grid %>% filter(coverage == 50)

find_crossing <- function(week, brr, target = 1) {
  ok <- is.finite(week) & is.finite(brr) & brr > 0
  week <- week[ok]
  brr <- brr[ok]
  o <- order(week)
  week <- week[o]
  brr <- brr[o]
  delta <- log(brr) - log(target)
  exact <- which(delta == 0)
  if (length(exact)) return(week[[exact[[1]]]])
  crossing <- which(delta[-length(delta)] * delta[-1] < 0)
  if (!length(crossing)) return(NA_real_)
  i <- crossing[[1]]
  week[[i]] + (0 - delta[[i]]) * (week[[i + 1]] - week[[i]]) /
    (delta[[i + 1]] - delta[[i]])
}

crossover_data <- line_data %>%
  filter(risk_assumption == "Base risk") %>%
  group_by(age_group, setting, mechanism) %>%
  summarise(
    week = find_crossing(week, BRR),
    BRR = 1,
    .groups = "drop"
  ) %>%
  filter(is.finite(week))

mechanism_cols <- c(
  "Disease blocking only" = "#2a9d8f",
  "Disease and infection blocking" = "#c44e52"
)

## Build an effect-size style "determinants of BRR" plot ---------------------------------

# Reference scenario: Low transmission, age 18-64, disease-only, week 2, 50% coverage, base risk
ref_row <- coverage_grid %>%
  filter(
    coverage == 50,
    week == 1,
    setting == "Low",
    age_group == "18-64",
    mechanism == "Disease blocking only",
    risk_assumption == "Base risk"
  )
ref_brr <- if (nrow(ref_row)) ref_row$BRR[[1]] else NA_real_
if (!is.finite(ref_brr) || ref_brr <= 0) stop("Reference BRR is missing or non-positive")

get_brr <- function(coverage = 50, week = 1, setting = "Low",
                    age_group = "18-64", mechanism = "Disease blocking only",
                    risk_assumption = "Base risk") {
  row <- coverage_grid %>%
    filter(
      coverage == !!coverage,
      week == !!week,
      setting == !!setting,
      age_group == !!age_group,
      mechanism == !!mechanism,
      risk_assumption == !!risk_assumption
    )
  if (!nrow(row)) return(NA_real_)
  row$BRR[[1]]
}

# Comparisons to show (labelled as "Reference → Alternative")
comparisons <- tibble::tibble(
  group = c(
    "Transmission", "Transmission",
    "Age group",
    "Protection mechanism",
    "Campaign timing", "Campaign timing", "Campaign timing",
    "Coverage", "Coverage",
    "Risk assumption"
  ),
    label = c(
    "Low → Moderate", "Low → High",
    "18–64 → ≥65",
    "Disease only → D+I",
    "Week 1 → Week 8", "Week 1 → Week 16", "Week 1 → Week 24",
    "50% → 10%", "50% → 90%",
    "Base → serostatus-adjusted"
  )
)

# Compute alternative BRR for each comparison, keeping other factors at reference
comparisons <- comparisons %>%
  rowwise() %>%
  mutate(
    alt_brr = case_when(
      label == "Low → Moderate" ~ get_brr(coverage = 50, week = 2, setting = "Moderate"),
      label == "Low → High" ~ get_brr(coverage = 50, week = 2, setting = "High"),
      label == "18–64 → ≥65" ~ get_brr(coverage = 50, week = 1, setting = "Low", age_group = "≥65"),
      label == "Disease only → D+I" ~ get_brr(coverage = 50, week = 2, setting = "Low", mechanism = "Disease and infection blocking"),
      label == "Week 1 → Week 8" ~ get_brr(coverage = 50, week = 8),
      label == "Week 1 → Week 16" ~ get_brr(coverage = 50, week = 16),
      label == "Week 1 → Week 24" ~ get_brr(coverage = 50, week = 24),
      label == "50% → 10%" ~ get_brr(coverage = 10, week = 2),
      label == "50% → 90%" ~ get_brr(coverage = 90, week = 2),
      label == "Base → serostatus-adjusted" ~ get_brr(coverage = 50, week = 2, risk_assumption = "Serostatus-adjusted risk"),
      TRUE ~ NA_real_
    ),
    ratio = ifelse(is.finite(alt_brr) & alt_brr > 0, alt_brr / ref_brr, NA_real_)
  ) %>%
  ungroup()

# Prepare plotting positions (top-down order)
comparisons <- comparisons %>%
  mutate(y = rev(seq_len(n()))) %>%
  arrange(desc(y))

# Prepare segments: from 1 (reference) to ratio (alternative)
seg_data <- comparisons %>%
  mutate(x1 = 1, x2 = ratio)

# Colour by direction (increase vs decrease)
seg_data <- seg_data %>%
  mutate(direction = case_when(is.na(ratio) ~ "missing",
                               ratio > 1 ~ "higher",
                               ratio < 1 ~ "lower",
                               TRUE ~ "equal"))

cols_dir <- c(higher = "#c44e52", lower = "#2a9d8f", equal = "#6c6c6c", missing = "#d9d9d9")

p <- ggplot(seg_data) +
  geom_segment(aes(x = pmin(x1, x2, na.rm = TRUE), xend = pmax(x1, x2, na.rm = TRUE),
                   y = y, yend = y, colour = direction), linewidth = 1.2, na.rm = TRUE) +
  geom_point(aes(x = x1, y = y), shape = 21, fill = "white", size = 2.2, colour = "#333333") +
  geom_point(aes(x = x2, y = y, fill = direction), shape = 21, size = 3, colour = "#333333", na.rm = TRUE) +
  scale_colour_manual(values = cols_dir, guide = "none") +
  scale_fill_manual(values = cols_dir, guide = "none") +
  scale_x_log10(
    limits = c(0.08, 12),
    breaks = c(0.1, 0.3, 1, 3, 10),
    labels = c("0.1", "0.3", "1", "3", "10")
  ) +
  scale_y_continuous(expand = expansion(add = c(0.5, 0.5)), breaks = seg_data$y,
                     labels = paste0(seg_data$group, ": ", seg_data$label)) +
  labs(
    title = "B. Determinants of BRR: fold-change relative to reference",
    x = "Fold-change in median DALY BRR relative to reference",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0, face = "plain", size = 11),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.margin = margin(8, 12, 8, 8)
  )

# Save outputs
out_files <- c(
  "06_Results/figure2B_determinants_brr.png",
  "07_Final_Results/Main_Figures/figure2B_determinants_brr.png"
)
dir.create("07_Final_Results/Main_Figures", showWarnings = FALSE, recursive = TRUE)
for (out_file in out_files) {
  ggsave(out_file, p, width = 9, height = 6.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}

## B. Window for favourable outbreak-response vaccination (median BRR > 1)

# Use interpolated 50% coverage grid
line50 <- coverage_grid %>% filter(coverage == 50)

# Compute latest week where BRR > 1 for each combination
window_df <- line50 %>%
  group_by(age_group, setting, mechanism, risk_assumption) %>%
  summarise(
    latest_week = if (any(is.finite(BRR) & BRR > 1)) max(week[is.finite(BRR) & BRR > 1], na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

# Prepare base and adjusted sets
base_win <- window_df %>% filter(risk_assumption == "Base risk")
adj_win <- window_df %>% filter(risk_assumption == "Serostatus-adjusted risk")

# Order y levels: for each setting (Low, Moderate, High) list Disease-only then D+I
y_levels <- c(
  "Low|Disease blocking only", "Low|Disease and infection blocking",
  "Moderate|Disease blocking only", "Moderate|Disease and infection blocking",
  "High|Disease blocking only", "High|Disease and infection blocking"
)

prep_y <- function(df) {
  df %>% mutate(y = paste(setting, mechanism, sep = "|"),
                y = factor(y, levels = y_levels))
}

base_win <- prep_y(base_win)
adj_win <- prep_y(adj_win)

## Compute windows for median BRR>1 and Pr(BRR>1) >= 0.9
# Interpolate probability grid at 50% coverage
prob_base_grid <- interpolate_coverage(summary_data, "brr_base_prob_gt1") %>% mutate(risk_assumption = "Base risk")
prob_adj_grid <- interpolate_coverage(summary_data, "brr_adj_prob_gt1") %>% mutate(risk_assumption = "Serostatus-adjusted risk")
prob_grid <- bind_rows(prob_base_grid, prob_adj_grid) %>%
  mutate(prob = ifelse(is.finite(BRR), BRR, NA_real_))

win_df_med <- coverage_grid %>% filter(coverage == 50) %>%
  group_by(age_group, setting, mechanism, risk_assumption) %>%
  summarise(latest_med = if (any(is.finite(BRR) & BRR > 1)) max(week[is.finite(BRR) & BRR > 1]) else NA_real_, .groups = "drop")
win_df_prob <- prob_grid %>% filter(coverage == 50) %>%
  group_by(age_group, setting, mechanism, risk_assumption) %>%
  summarise(latest_p90 = if (any(is.finite(prob) & prob >= 0.9)) max(week[is.finite(prob) & prob >= 0.9]) else NA_real_, .groups = "drop")

win_df <- full_join(win_df_med, win_df_prob, by = c("age_group","setting","mechanism","risk_assumption"))

# Build plotting y-order: for each setting (Low, Moderate, High), within each mechanism (DB then D+I),
# show Base (upper) then Adjusted (lower).
settings <- c("Low","Moderate","High")
mechs <- c("Disease blocking only","Disease and infection blocking")
risk_levels <- c("Base risk","Serostatus-adjusted risk")

plot_rows <- expand.grid(setting = settings, mechanism = mechs, risk_assumption = risk_levels, stringsAsFactors = FALSE)
plot_rows <- plot_rows[order(match(plot_rows$setting, settings), match(plot_rows$mechanism, mechs), match(plot_rows$risk_assumption, risk_levels)), ]
plot_rows$y <- seq_len(nrow(plot_rows))
plot_rows <- plot_rows %>% mutate(y = rev(y))

plot_df <- plot_rows %>%
  left_join(win_df, by = c("setting","mechanism","risk_assumption")) %>%
  mutate(mech_col = mechanism)

# Plot: pale (median) bar behind, saturated (Pr>=90%) bar on top (shorter), colored by mechanism.
p_window <- ggplot(plot_df) +
  geom_segment(data = filter(plot_df, !is.na(latest_med)),
               aes(x = 1, xend = latest_med, y = y, yend = y, colour = mech_col),
               linewidth = 10, alpha = 0.35, lineend = "butt") +
  geom_segment(data = filter(plot_df, !is.na(latest_p90)),
               aes(x = 1, xend = latest_p90, y = y, yend = y, colour = mech_col),
               linewidth = 6, alpha = 1.0, lineend = "butt") +
  scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
  scale_x_continuous(name = "Vaccination campaign start week",
                     breaks = c(1,8,16,24,32,42,52), limits = c(1,52)) +
  scale_y_continuous(breaks = plot_df$y, labels = paste0(plot_df$setting, "  ", ifelse(plot_df$mechanism=="Disease blocking only","DB","D+I"), "  (", ifelse(plot_df$risk_assumption=="Base risk","Base","Adj"), ")"), expand = expansion(add = c(0.6,0.6))) +
  facet_grid(age_group ~ ., scales = "free_y", space = "free_y") +
  labs(title = "B. Window for favourable outbreak-response vaccination (median BRR > 1)",
       subtitle = "Reference: 50% coverage; pale = median BRR>1 window, saturated = Pr(BRR>1) ≥ 90%",
       y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_text(hjust = 0, face = "plain", size = 11),
        axis.title.x = element_text(face = "bold"), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11),
        plot.margin = margin(8,12,8,8))

# Save window plot
out_files_w <- c(
  "06_Results/figure2B_window_favourable.png",
  "07_Final_Results/Main_Figures/figure2B_window_favourable.png"
)
for (out_file in out_files_w) {
  ggsave(out_file, p_window, width = 10, height = 6, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}

## Optional: Time-series BRR with coverage points at selected weeks (1,13,26,39)
point_weeks <- c(1, 13, 26, 39)
point_data <- coverage_grid %>%
  filter(coverage %in% c(10, 50, 90), week %in% point_weeks)

shapes_map <- c(`10` = 1, `50` = 16, `90` = 2) # 1=open circle,16=filled circle,2=open triangle
sizes_map <- c(`10` = 2, `50` = 3.5, `90` = 3)

p_time <- ggplot() +
  geom_line(
    data = line_data,
    aes(x = week, y = BRR, colour = mechanism, linetype = risk_assumption,
        group = interaction(mechanism, risk_assumption)),
    linewidth = 1.05, na.rm = TRUE
  ) +
  geom_point(
    data = point_data,
    aes(x = week, y = BRR, colour = mechanism, shape = factor(coverage)),
    size = 3, stroke = 0.9, na.rm = TRUE
  ) +
  scale_shape_manual(values = shapes_map, name = "Coverage") +
  scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
  scale_linetype_manual(
    values = c("Base risk" = "solid", "Serostatus-adjusted risk" = "dashed"),
    name = "Risk assumption"
  ) +
  scale_x_continuous(
    name = "Vaccination campaign start week",
    breaks = c(1, 13, 26, 39), limits = c(1, 52), expand = c(0, 0)
  ) +
  scale_y_log10(
    name = "Median DALY benefit-risk ratio (BRR; log scale)",
    breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"),
    limits = c(0.01, NA), expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "B. BRR over campaign start week — coverage points at selected weeks",
    subtitle = "Points: ○=10%, ●=50%, △=90% (colours = protection mechanism)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(colour = "#e1e0d9", linewidth = 0.3),
    panel.border = element_rect(colour = "#c8c8c8", fill = NA, linewidth = 0.35),
    panel.spacing = unit(10, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14, colour = "#1a1a1a"),
    axis.title = element_text(face = "bold", size = 13, colour = "#1a1a1a"),
    axis.text = element_text(size = 11.5, colour = "#1a1a1a"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(colour = "#5a5a5a", size = 12),
    plot.margin = margin(10, 12, 8, 8)
  ) +
  facet_grid(age_group ~ setting)

# Save time-series plot
out_files_time <- c(
  "06_Results/figure2B_timeseries_points.png",
  "07_Final_Results/Main_Figures/figure2B_timeseries_points.png"
)
for (out_file in out_files_time) {
  ggsave(out_file, p_time, width = 12, height = 9.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}

## C. Effect of coverage on per-vaccinee BRR (inset-style plot)
# Interpolate BRR across coverage 10-90 for selected representative weeks.
interpolate_coverage_seq <- function(data, value_col, covs = seq(10, 90, by = 1)) {
  data %>%
    group_by(age_group, setting, mechanism, week) %>%
    group_modify(~ {
      x <- as.numeric(sub("cov", "", .x$Coverage))
      y <- .x[[value_col]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 2) return(tibble(coverage = covs, BRR = NA_real_, week = unique(.x$week)))
      o <- order(x[ok])
      x <- x[ok][o]
      y <- y[ok][o]
      tibble(coverage = covs, BRR = approx(x, y, xout = covs, rule = 1)$y, week = unique(.x$week))
    }) %>%
    ungroup()
}

# Generate fine coverage grid for base risk
base_seq <- interpolate_coverage_seq(summary_data, "brr_base_med") %>% mutate(risk_assumption = "Base risk")

# Choose requested representative weeks and map to nearest available weeks in data
requested_weeks <- c(2, 13, 26)
available_weeks <- sort(unique(line_data$week))
map_nearest <- function(req) available_weeks[which.min(abs(available_weeks - req))]
actual_weeks <- vapply(requested_weeks, map_nearest, numeric(1))

cov_plot_df <- base_seq %>%
  filter(week %in% actual_weeks, age_group == "18-64", setting == "Moderate") %>%
  group_by(age_group, setting, mechanism, week) %>%
  mutate(brr50 = BRR[coverage == 50]) %>%
  ungroup() %>%
  mutate(ratio = ifelse(is.finite(BRR) & is.finite(brr50) & brr50 > 0, BRR / brr50, NA_real_)) %>%
  mutate(week_label = paste0("wk", week))

p_coverage <- ggplot(cov_plot_df, aes(x = coverage, y = ratio, colour = mechanism, linetype = factor(week_label))) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "#b8b8b3") +
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  scale_x_continuous(name = "Coverage (%)", breaks = c(10,30,50,70,90), limits = c(10,90)) +
  scale_y_continuous(name = "BRR relative to 50% coverage", trans = "identity") +
  scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
  labs(title = "C. Effect of coverage on per-vaccinee BRR",
       subtitle = paste0("Representative weeks (requested → actual): ", paste(requested_weeks, "→", actual_weeks, collapse = ", ")))+
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Save coverage inset plot
out_files_cov <- c(
  "06_Results/figure2B_coverage_inset.png",
  "07_Final_Results/Main_Figures/figure2B_coverage_inset.png"
)
for (out_file in out_files_cov) {
  ggsave(out_file, p_coverage, width = 5.5, height = 4, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}

## Combined 3x2 panel: ribbons for coverage effect and lines for 50% (base solid, adjusted dashed)
cov_wide <- coverage_grid %>%
  filter(coverage %in% c(10,50,90)) %>%
  select(age_group, setting, mechanism, week, coverage, BRR, risk_assumption) %>%
  tidyr::pivot_wider(names_from = coverage, values_from = BRR, names_prefix = "cov_")

# Separate base and adjusted mid lines
cov_base <- cov_wide %>% filter(risk_assumption == "Base risk")
cov_adj  <- cov_wide %>% filter(risk_assumption == "Serostatus-adjusted risk")

# Prepare edge lines for cov_10 and cov_90 so we can show which edge is high/low
cov_edges <- cov_base %>%
  tidyr::pivot_longer(cols = tidyselect::starts_with("cov_"), names_to = "cov_lbl", values_to = "val") %>%
  filter(cov_lbl %in% c("cov_10", "cov_90")) %>%
  mutate(coverage_label = ifelse(cov_lbl == "cov_10", "10%", "90%"))

p_combined <- ggplot() +
  geom_ribbon(data = cov_base, aes(x = week, ymin = cov_10, ymax = cov_90, fill = mechanism),
              alpha = 0.22, inherit.aes = FALSE) +
  # thin edge lines to show which side corresponds to 10% vs 90%
  geom_line(data = cov_edges, aes(x = week, y = val, colour = mechanism, linetype = coverage_label), size = 0.6, inherit.aes = FALSE) +
  geom_line(data = cov_base, aes(x = week, y = cov_50, colour = mechanism), linewidth = 1.1) +
  geom_line(data = cov_adj, aes(x = week, y = cov_50, colour = mechanism), linewidth = 0.9, linetype = "dashed") +
  scale_fill_manual(values = mechanism_cols, guide = "none") +
  scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
  scale_linetype_manual(values = c("10%" = "dotted", "90%" = "dotdash"), name = "Edge coverage") +
  scale_x_continuous(name = "Vaccination campaign start week", breaks = c(1,13,26,39), limits = c(1,52)) +
  scale_y_log10(name = "Median DALY BRR (log)", breaks = c(0.1,1,10,100), labels = c("0.1","1","10","100"), limits = c(0.01, NA)) +
  facet_grid(age_group ~ setting) +
  labs(title = "B. BRR trajectories with coverage effect (ribbon = 10–90%, line = 50%)",
       subtitle = "Solid = base risk (50%); dashed = serostatus-adjusted (50%)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(), panel.grid.major.y = element_line(colour = "#e1e0d9", linewidth = 0.3),
        strip.text = element_text(face = "bold"), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16))

# Save combined panel
out_files_comb <- c(
  "06_Results/figure2B_combined_ribbon.png",
  "07_Final_Results/Main_Figures/figure2B_combined_ribbon.png"
)
for (out_file in out_files_comb) {
  ggsave(out_file, p_combined, width = 12, height = 9.5, dpi = 300, bg = "#fcfcfb")
  message("Saved: ", out_file)
}

## D. ORI-style: BRR versus remaining epidemic infection risk (remaining attack rate)
# x-axis = remaining symptomatic attack rate (percent) in the absence of vaccination
# colour = protection mechanism; linetype = base vs serostatus-adjusted; shape = setting (Low/Moderate/High)

ori_rdata_candidates <- c(
  "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep.RData",
  "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep_ve0.RData",
  "00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep_ve989.RData",
  # fall back to 01_Data if processed files live there
  "01_Data/draw_level_xy_serostatus_finite_weeksweep.RData",
  "01_Data/draw_level_xy_serostatus_finite_weeksweep_ve0.RData",
  "01_Data/draw_level_xy_serostatus_finite_weeksweep_ve989.RData",
  # complete draw-level object also present in 01_Data
  "01_Data/draw_level_xy_serostatus_finite.RData",
  "01_Data/draw_level_xy_serostatus.RData"
)
found_rdata <- ori_rdata_candidates[file.exists(ori_rdata_candidates)]

if (length(found_rdata) == 0) {
  message("[ORI plot] draw-level processed data not found in 00_Data/0_2_Processed/; skipping ORI remaining-risk plot.")
} else {
  load(found_rdata[[1]]) # expects object draw_level_xy_serostatus or similar
  # tolerant: accept either draw_level_xy_serostatus or all_weeks_brr / draw-level object
  if (exists("draw_level_xy_serostatus")) ori_draws <- draw_level_xy_serostatus else if (exists("all_weeks_brr")) ori_draws <- all_weeks_brr else ori_draws <- NULL

  if (is.null(ori_draws)) {
    message("[ORI plot] expected object 'draw_level_xy_serostatus' or 'all_weeks_brr' not found in RData; skipping.")
  } else {
    # choose the symptomatic attack-rate metric if present (symp_10k or inf_10k)
    if ("symp_10k" %in% names(ori_draws)) {
      ori_draws <- ori_draws %>% rename(symp_per10k = symp_10k)
    } else if ("inf_10k" %in% names(ori_draws)) {
      ori_draws <- ori_draws %>% rename(symp_per10k = inf_10k)
    } else if ("x_10k_base" %in% names(ori_draws)) {
      # x_10k_base is the outcome-specific x denominator (DALY/SAE/death) per 10k;
      # when outcome == "DALY" this is DALY-per-10k not symptomatic attack; fall back but warn.
      ori_draws <- ori_draws %>% rename(symp_per10k = x_10k_base)
      warning("[ORI plot] using 'x_10k_base' as a proxy for remaining infection risk (may not be symptomatic attack rate)")
    } else {
      message("[ORI plot] no suitable per-10k infection/symptomatic column found; skipping ORI remaining-risk plot.")
      ori_draws <- NULL
    }
  }

  if (!is.null(ori_draws)) {
    # detect best candidate column for remaining infection risk (per-10k)
    cols <- names(ori_draws)
    prefer_cols <- c("symp_nv_10k", "symp_10k", "inf_10k", "AR_travel", "symp_per10k", "x_10k_base")
    chosen <- intersect(prefer_cols, cols)[1]
    if (is.na(chosen) || is.null(chosen)) chosen <- cols[grepl("symp|inf|attack|ar|AR", cols, ignore.case = TRUE)][1]

    if (is.na(chosen) || is.null(chosen)) {
      message("[ORI plot] no suitable remaining-risk column found in draw-level object; skipping ORI plot.")
    } else {
      if (chosen == "x_10k_base") message("[ORI plot] using 'x_10k_base' as a proxy for remaining infection risk (not ideal)")

      ori_summary <- ori_draws %>%
        filter(RR_seropos == 0, outcome == "DALY", AgeCat %in% c("18-64","65+")) %>%
        mutate(age_group = factor(AgeCat, levels = c("18-64","65+"), labels = c("18-64","≥65")),
               setting = factor(setting, levels = c("Low","Moderate","High")),
               week = as.numeric(week),
               coverage = as.numeric(ifelse(!is.null(Coverage), sub("cov", "", Coverage), NA))) %>%
        group_by(age_group, setting, week) %>%
        summarise(remaining_per10k = median(.data[[chosen]], na.rm = TRUE), .groups = "drop") %>%
        mutate(remaining_pct = remaining_per10k / 100) # per10k -> percent

      # BRR at 50% coverage (base and adjusted)
      brr50 <- coverage_grid %>% filter(coverage == 50) %>%
        select(age_group, setting, mechanism, week, risk_assumption, BRR) %>%
        rename(risk = risk_assumption, brr = BRR)

      # join and remove rows without remaining-risk
      ori_plot_df <- brr50 %>%
        left_join(ori_summary %>% select(age_group, setting, week, remaining_pct), by = c("age_group","setting","week")) %>%
        filter(!is.na(remaining_pct)) %>%
        mutate(setting = factor(setting, levels = c("Low","Moderate","High")))

      # use alpha to denote Low/Moderate/High (paler = Low)
      alpha_map <- c(Low = 0.5, Moderate = 0.75, High = 1)

      p_ori <- ggplot(ori_plot_df, aes(x = remaining_pct, y = brr, colour = mechanism, linetype = risk, group = interaction(mechanism, risk))) +
        geom_line(linewidth = 1.05, na.rm = TRUE) +
        geom_point(aes(alpha = setting), size = 2.6, shape = 21, fill = "white", colour = "black", na.rm = TRUE) +
        scale_alpha_manual(values = alpha_map, name = "2022 intensity") +
        scale_colour_manual(values = mechanism_cols, name = "Protection mechanism") +
        scale_linetype_manual(values = c("Base risk" = "solid", "Serostatus-adjusted risk" = "dashed"), name = "Risk assumption") +
        scale_x_continuous(name = "Remaining infection risk after campaign initiation (%)", limits = c(0, 6), breaks = seq(0,6,by=1)) +
        scale_y_log10(name = "Median DALY BRR (log)", breaks = c(0.1,1,10,100), labels = c("0.1","1","10","100"), limits = c(0.01, NA)) +
        geom_hline(yintercept = 1, linetype = "dotted", colour = "#444444") +
        facet_grid(. ~ age_group) +
        labs(title = "B. Outbreak-response benefit–risk by remaining epidemic infection risk",
             subtitle = "Coverage fixed at 50%; colour = mechanism, solid/dashed = base/serostatus-adjusted") +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")

      out_files_ori <- c(
        "06_Results/figure2B_remaining_risk_vs_brr.png",
        "07_Final_Results/Main_Figures/figure2B_remaining_risk_vs_brr.png"
      )
      for (of in out_files_ori) {
        ggsave(of, p_ori, width = 10, height = 5, dpi = 300, bg = "#fcfcfb")
        message("Saved: ", of)
      }
    }
  }
}
