# =============================================================================
# 27_main_table2_travel_death.R
#
# Main Table 2 (Death companion): Absolute death benefits and vaccine-
# attributable risks for travel vaccination, by trip duration. Same structure
# as 20_main_table2_travel_daly.R, with two differences:
#   1) 18-64 is EXCLUDED -- vaccine-attributable death risk (excess_10k_death)
#      is identically 0 for 18-64 in psa_df, making brr_death NA for that age
#      group (same reasoning as 25_main_table1_death.R's exclusion for the
#      outbreak-response side). Only 65+ rows are shown.
#   2) 2 decimal places (already used here, same as the DALY/SAE companions)
#      -- kept as-is since Death magnitudes (benefit ~0.01, risk ~0.6) need at
#      least 2 decimals to stay legible.
#
# Output: 06_Results/Main_Table_2_Travel_Death.docx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(writexl)
})

setting_levels  <- c("Low", "Moderate", "High")
duration_levels <- c("7d", "14d", "30d", "90d")

fmt_ci2 <- function(med, lo, hi) sprintf("%.2f (%.2f\u2013%.2f)", med, lo, hi)

load("01_Data/psa_df_bra_travel_finite.RData")  # -> psa_df (draw-level)
load("01_Data/setting_key.RData")               # -> setting_key

psa_df <- psa_df %>%
  mutate(setting = unname(setting_key[state])) %>%
  filter(!is.na(setting), age_group == "65+") %>%
  mutate(
    setting   = factor(setting, levels = setting_levels),
    days      = factor(days, levels = duration_levels),
    age_group = factor(age_group, levels = "65+")
  )

# Death only, 65+ only -- same construction as travel_long in
# 11_travel_headline_figures.R / 14_headline_summary_tables.R, filtered to
# just the Death outcome.
travel_death <- psa_df %>%
  transmute(draw, state, setting, age_group, days,
            benefit = averted_10k_death, risk = excess_10k_death,
            brr = ifelse(is.finite(excess_10k_death) & excess_10k_death > 0, averted_10k_death / excess_10k_death, NA_real_))

# ---- Setting x Duration (65+ only): Benefit, Risk ---------------------------
group_summary <- travel_death %>%
  group_by(setting, age_group, days) %>%
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
    .groups = "drop"
  ) %>%
  mutate(
    Benefit = fmt_ci2(Benefit_med, Benefit_lo, Benefit_hi),
    Risk    = fmt_ci2(Risk_med, Risk_lo, Risk_hi)
    ,BRR     = fmt_ci2(BRR_med, BRR_lo, BRR_hi)
  )

sub_cols <- c("Benefit", "Risk", "BRR")

build_block <- function(d) {
  group_summary %>%
    filter(days == d) %>%
    select(setting, age_group, dplyr::all_of(sub_cols)) %>%
    rename_with(~ paste0(.x, "___", d), dplyr::all_of(sub_cols))
}

blocks <- lapply(duration_levels, build_block)
wide <- Reduce(function(a, b) dplyr::left_join(a, b, by = c("setting", "age_group")), blocks) %>%
  arrange(setting, age_group) %>%
  rename(Transmission = setting, `Age group` = age_group)

print(wide, n = 10, width = Inf)

# ---- Assemble two-row-header flextable --------------------------------------
header_labels <- c(Transmission = "Transmission", `Age group` = "Age group")
for (d in duration_levels) {
  for (sc in sub_cols) {
    header_labels[paste0(sc, "___", d)] <- sc
  }
}

ft <- flextable::flextable(wide) %>%
  flextable::set_header_labels(values = as.list(header_labels)) %>%
  flextable::add_header_row(
    values = c("", "", duration_levels),
    colwidths = c(1, 1, rep(length(sub_cols), length(duration_levels)))
  ) %>%
  flextable::set_caption(caption = "Absolute death benefits and vaccine-attributable risks for travel vaccination, by trip duration (65+ only)") %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Transmission", "Age group")) %>%
  flextable::valign(j = c("Transmission", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par("Main Table 2 (Death)", style = "heading 1") %>%
  officer::body_add_par(
    "Absolute death benefits and vaccine-attributable risks for travel vaccination, by trip duration",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par(
    paste(
      "Values are median (95% uncertainty interval) across 1000 posterior draws, per 10,000 vaccinated.",
      "18-64 is omitted: vaccine-attributable death risk is 0 for this age group, leaving the death BRR",
      "undefined. Benefit = deaths averted; Risk = vaccine-attributable deaths. Risk is identical across all",
      "four durations within an age group (vaccine safety risk depends on age, not trip length) -- only Benefit",
      "changes with duration. Entry timing pooled across the within-draw sampled entry days. Travel vaccination",
      "has only ever been simulated under the disease and infection blocking mechanism (no VE0/disease-blocking-",
      "only comparison, unlike Main Table 1's outbreak-response results) and has no serostatus-adjustment split",
      "(single risk value, not base/adjusted)."
    ),
    style = "Normal"
  )

print(doc, target = "06_Results/Main_Table_2_Travel_Death.docx")
message("Saved: 06_Results/Main_Table_2_Travel_Death.docx")

writexl::write_xlsx(wide, "06_Results/Main_Table_2_Travel_Death.xlsx")
message("Saved: 06_Results/Main_Table_2_Travel_Death.xlsx")
