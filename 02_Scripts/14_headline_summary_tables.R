# =============================================================================
# 14_headline_summary_tables.R
#
# Compact headline summary: Outcome x Setting x Age group -> median Benefit
# (per 10,000 vaccinated), median Risk (per 10,000), median BRR, Pr(BRR>1).
# No duration/week breakdown -- see the two heatmaps
# (finite_weeksweep_heatmap.R, 13_travel_duration_entryday_heatmap_setting.R)
# for how BRR varies across timing/coverage/duration; this table is the
# single-number-per-stratum companion, not a replacement.
#
#   Table A (outbreak response): base case = fixed week-2/50%-coverage
#   scenario (seropositive vaccinees assumed zero vaccine-attributable risk),
#   median + Pr(BRR>1) taken across posterior draws at that one scenario --
#   the same convention headline_benefit_risk_plane_*.png already uses.
#
#   Table B (travel vaccination): duration (7d/14d/30d/90d) pivoted into
#   column groups -- the same compression already used for mechanism on the
#   ori table -- so each outcome's table stays at Setting x Age group (6
#   rows) while still showing all 4 simulated durations. Entry timing is
#   pooled across the within-draw sampled entry days. Previously this pooled
#   across the WHOLE duration(7-180d) x entry-week(1-52) surface, but that
#   surface's per-cell values include cells extrapolated past what was
#   actually simulated for some duration/entry-week combinations (see
#   figure2a_travel_infection_risk_surface.png's investigation: only
#   duration<=90d x entry week<=40 is fully computed for every state without
#   extrapolation) -- silently mixing real and extrapolated values into one
#   "pooled" number. Using travel_long (draw-level psa_df) at the 4
#   actually-simulated durations avoids that entirely.
#
# Output:
#   - 06_Results/headline_summary_outbreak.docx
#   - 06_Results/headline_summary_travel.docx (DALY, main)
#   - 06_Results/headline_summary_travel_Death.docx
#   - 06_Results/headline_summary_travel_SAE.docx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
library(dplyr)
library(flextable)
library(officer)
library(writexl)

age_levels <- c("18-64", "65+")
setting_levels <- c("Low", "Moderate", "High")
outcome_levels <- c("SAE", "Death", "DALY")

fmt_num  <- function(x) sprintf("%.1f", x)  # ori (Table A) -- 1 decimal place
fmt_num2 <- function(x) sprintf("%.2f", x)  # travel (Table B) -- 2 decimal places
                                             # (BRR values are small, ~0.01-0.05;
                                             # 1 decimal collapsed most to "0.0")
fmt_pct <- function(p) sprintf("%.0f%%", 100 * p)
fmt_ci  <- function(med, lo, hi) sprintf("%.1f (%.1f–%.1f)", med, lo, hi)   # ori -- 1 decimal
fmt_ci2 <- function(med, lo, hi) sprintf("%.2f (%.2f–%.2f)", med, lo, hi)   # travel -- 2 decimals

build_docx <- function(df, out_file, caption_text) {
  merge_cols <- intersect(c("Outcome", "Setting", "Age group"), names(df))
  ft <- flextable(df) %>%
    set_caption(caption = caption_text) %>%
    merge_v(j = merge_cols) %>%
    valign(j = merge_cols, valign = "top") %>%
    autofit() %>%
    theme_vanilla()
  doc <- read_docx() %>% body_add_flextable(ft)
  print(doc, target = out_file)
  cat("Saved:", out_file, "\n")

  xlsx_file <- sub("\\.docx$", ".xlsx", out_file)
  writexl::write_xlsx(df, xlsx_file)
  cat("Saved:", xlsx_file, "\n")
}

# =============================================================================
# Table A -- Outbreak response (base case: week 2, 50% coverage, seroadjusted)
# =============================================================================

mechanism_levels <- c("Disease blocking only", "Disease and infection blocking")

load("01_Data/draw_level_xy_serostatus_finite.RData")  # -> draw_level_xy_serostatus

outbreak_summary <- draw_level_xy_serostatus %>%
  filter(
    RR_seropos == 0,
    AgeCat %in% age_levels,
    !(outcome == "Death" & AgeCat == "18-64")  # same exclusion as 10_streamline_mechanism_figure.R
  ) %>%
  group_by(setting, AgeCat, outcome, VE_label) %>%
  summarise(
    Benefit_med = median(averted_10k, na.rm = TRUE),
    Benefit_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
    Benefit_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),
    Risk_med    = median(x_10k_adj, na.rm = TRUE),
    Risk_lo     = quantile(x_10k_adj, 0.025, na.rm = TRUE),
    Risk_hi     = quantile(x_10k_adj, 0.975, na.rm = TRUE),
    BRR_med     = median(brr_adj, na.rm = TRUE),
    BRR_lo      = quantile(brr_adj, 0.025, na.rm = TRUE),
    BRR_hi      = quantile(brr_adj, 0.975, na.rm = TRUE),
    prob_gt1    = mean(brr_adj > 1, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  rename(age_group = AgeCat, Mechanism = VE_label) %>%
  mutate(
    setting   = factor(setting, levels = setting_levels),
    age_group = factor(age_group, levels = age_levels),
    outcome   = factor(outcome, levels = outcome_levels),
    Mechanism = factor(Mechanism, levels = mechanism_levels)
  ) %>%
  arrange(outcome, setting, age_group, Mechanism) %>%
  transmute(
    Outcome = outcome, Setting = setting, `Age group` = age_group, Mechanism = Mechanism,
    `Benefit per 10,000` = fmt_ci(Benefit_med, Benefit_lo, Benefit_hi),
    `Risk per 10,000`    = fmt_ci(Risk_med, Risk_lo, Risk_hi),
    BRR                  = fmt_ci(BRR_med, BRR_lo, BRR_hi),
    `Pr(BRR>1)`          = fmt_pct(prob_gt1)
  )

print(outbreak_summary, n = 50)

build_docx(
  outbreak_summary, "06_Results/headline_summary_outbreak.docx",
  paste(
    "Outbreak-response vaccination: headline benefit-risk summary, by mechanism (disease-blocking only vs.",
    "disease and infection blocking).",
    "Base case = fixed week-2/50%-coverage scenario, seroadjusted (seropositive vaccinees assumed zero",
    "vaccine-attributable risk). Benefit/Risk/BRR = median across 1000 posterior draws; Pr(BRR>1) = draw-level probability."
  )
)


# =============================================================================
# Table B -- Travel vaccination (base case: 30-day trip, entry timing pooled)
#
# No Mechanism column here (unlike Table A): the travel PSA has only ever
# simulated the "disease and infection blocking" mechanism -- there is no
# VE0/"disease blocking only" sweep for travel to report, so adding a
# Mechanism axis would be silently misleading (implying a comparison that
# was never run). If a VE0 travel scenario is added later, extend this
# table the same way Table A was extended.
# =============================================================================

load("01_Data/psa_df_bra_travel_finite.RData")  # -> psa_df (draw-level)
load("01_Data/setting_key.RData")               # -> setting_key

psa_df <- psa_df %>%
  mutate(setting = unname(setting_key[state])) %>%
  filter(!is.na(setting)) %>%
  mutate(
    setting   = factor(setting, levels = setting_levels),
    days      = factor(days, levels = c("7d", "14d", "30d", "90d")),
    age_group = factor(age_group, levels = age_levels)
  )

# One row per draw x state x age x days x outcome (same construction as
# 11_travel_headline_figures.R's travel_long).
travel_long <- bind_rows(
  psa_df %>% transmute(draw, state, setting, age_group, days,
                        outcome = "DALY", benefit = daly_averted, risk = daly_sae, brr = brr_daly),
  psa_df %>% transmute(draw, state, setting, age_group, days,
                        outcome = "SAE", benefit = averted_10k_sae, risk = excess_10k_sae, brr = brr_sae),
  psa_df %>% transmute(draw, state, setting, age_group, days,
                        outcome = "Death", benefit = averted_10k_death, risk = excess_10k_death, brr = brr_death)
)

duration_levels <- c("7d", "14d", "30d", "90d")

travel_all <- travel_long %>%
  group_by(outcome, setting, age_group, days) %>%
  summarise(
    Benefit_med = median(benefit, na.rm = TRUE),
    Benefit_lo  = quantile(benefit, 0.025, na.rm = TRUE),
    Benefit_hi  = quantile(benefit, 0.975, na.rm = TRUE),
    Risk_med    = median(risk, na.rm = TRUE),
    Risk_lo     = quantile(risk, 0.025, na.rm = TRUE),
    Risk_hi     = quantile(risk, 0.975, na.rm = TRUE),
    BRR_med     = median(brr, na.rm = TRUE),
    BRR_lo      = quantile(brr, 0.025, na.rm = TRUE),
    BRR_hi      = quantile(brr, 0.975, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  mutate(outcome = factor(outcome, levels = outcome_levels))

# Pr(BRR>1) dropped from the table (too dense, same reasoning as the ori table).
travel_metric_cols   <- c("Benefit", "Risk", "BRR")
travel_metric_labels <- c(
  Benefit  = "Benefit\n(per 10,000)",
  Risk     = "Risk\n(per 10,000)",
  BRR      = "BRR"
)

build_duration_block <- function(outcome_filter, days_filter) {
  travel_all %>%
    filter(outcome == outcome_filter, days == days_filter) %>%
    select(setting, age_group,
           Benefit_med, Benefit_lo, Benefit_hi,
           Risk_med, Risk_lo, Risk_hi,
           BRR_med, BRR_lo, BRR_hi) %>%
    rename_with(~ paste0(.x, "_", days_filter), -c(setting, age_group))
}

# Outcome-split (DALY main, Death/SAE companion files), duration pivoted to
# columns within each -- Setting x Age group stays at 6 rows per file.
build_travel_outcome_table <- function(outcome_filter) {
  blocks <- lapply(duration_levels, function(d) build_duration_block(outcome_filter, d))
  df <- Reduce(function(a, b) dplyr::left_join(a, b, by = c("setting", "age_group")), blocks) %>%
    mutate(setting = factor(setting, levels = setting_levels)) %>%
    arrange(setting, age_group)

  out <- df %>% select(setting, age_group)
  for (d in duration_levels) {
    out[[paste0("Benefit_", d)]] <- fmt_ci2(df[[paste0("Benefit_med_", d)]], df[[paste0("Benefit_lo_", d)]], df[[paste0("Benefit_hi_", d)]])
    out[[paste0("Risk_",    d)]] <- fmt_ci2(df[[paste0("Risk_med_",    d)]], df[[paste0("Risk_lo_",    d)]], df[[paste0("Risk_hi_",    d)]])
    out[[paste0("BRR_",     d)]] <- fmt_ci2(df[[paste0("BRR_med_",     d)]], df[[paste0("BRR_lo_",     d)]], df[[paste0("BRR_hi_",     d)]])
  }

  out %>% dplyr::rename(Setting = setting, `Age group` = age_group)
}

build_travel_docx <- function(df, out_file, caption_text) {
  header_labels <- c(
    Setting = "Setting", `Age group` = "Age group",
    unlist(lapply(duration_levels, function(d) setNames(travel_metric_labels, paste0(names(travel_metric_labels), "_", d))))
  )

  ft <- flextable::flextable(df) %>%
    flextable::set_header_labels(values = as.list(header_labels)) %>%
    flextable::add_header_row(
      values = c("", "", duration_levels),
      colwidths = c(1, 1, rep(length(travel_metric_cols), length(duration_levels)))
    ) %>%
    flextable::theme_booktabs() %>%
    flextable::bold(part = "header") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::align(j = 1:2, align = "left", part = "all") %>%
    flextable::merge_v(j = c("Setting", "Age group")) %>%
    flextable::valign(j = c("Setting", "Age group"), valign = "top") %>%
    flextable::fontsize(size = 8, part = "all") %>%
    flextable::autofit()

  doc <- officer::read_docx() %>%
    officer::body_add_par(caption_text, style = "heading 2") %>%
    flextable::body_add_flextable(ft)

  print(doc, target = out_file)
  cat("Saved:", out_file, "\n")

  xlsx_file <- sub("\\.docx$", ".xlsx", out_file)
  writexl::write_xlsx(df, xlsx_file)
  cat("Saved:", xlsx_file, "\n")
}

travel_summary_daly <- build_travel_outcome_table("DALY")
print(travel_summary_daly, n = 50)

build_travel_docx(
  travel_summary_daly, "06_Results/headline_summary_travel.docx",
  "Travel vaccination: headline benefit-risk summary by duration -- DALY (see headline_summary_travel_{Death,SAE}.docx for the other outcomes). Entry timing pooled across the within-draw sampled entry days. Benefit/Risk/BRR = median across 1000 posterior draws; Pr(BRR>1) = draw-level probability."
)

for (oc in c("Death", "SAE")) {
  tbl <- build_travel_outcome_table(oc)
  print(tbl, n = 50)
  build_travel_docx(
    tbl, sprintf("06_Results/headline_summary_travel_%s.docx", oc),
    sprintf(
      "Travel vaccination: headline benefit-risk summary by duration -- %s. Entry timing pooled across the within-draw sampled entry days. Benefit/Risk/BRR = median across 1000 posterior draws; Pr(BRR>1) = draw-level probability.",
      oc
    )
  )
}
