# =============================================================================
# 25_main_table1_death.R
#
# Main Table 1 (Death companion): Absolute death benefits and vaccine-
# attributable risks for outbreak-response immunisation under the reference
# scenario (week 2, 50% coverage). Same structure as 19_main_table1_daly.R,
# with two differences:
#   1) 18-64 is EXCLUDED -- vaccine-attributable death risk (x_10k_base) is
#      identically 0 for 18-64 in this dataset, making brr_base/brr_adj
#      NaN/undefined for that age group (same exclusion convention as
#      14_headline_summary_tables.R's Table A and 10_streamline_mechanism_
#      figure.R). Only 65+ rows are shown.
#   2) 2 decimal places (not 1) -- Death benefit/risk/BRR values are all
#      O(0.01-2), so 1 decimal would collapse most values to 0.0-2.0 and lose
#      the distinctions that matter near the BRR=1 threshold.
#
# Output: 06_Results/Main_Table_1_Death.docx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(writexl)
})

setting_levels <- c("Low", "Moderate", "High")
mechanism_abbrev <- c("Disease blocking only" = "DB", "Disease and infection blocking" = "D+I")

fmt_ci <- function(med, lo, hi) sprintf("%.2f (%.2f\u2013%.2f)", med, lo, hi)

load("01_Data/draw_level_xy_serostatus_finite.RData")  # -> draw_level_xy_serostatus

d <- draw_level_xy_serostatus %>%
  filter(outcome == "Death", RR_seropos == 0, AgeCat == "65+")

# ---- Setting x Mechanism (65+ only): Benefit, Risk (base/adj), BRR (base/adj) ----
group_summary <- d %>%
  group_by(setting, AgeCat, VE_label) %>%
  summarise(
    Benefit_med = median(averted_10k, na.rm = TRUE),
    Benefit_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
    Benefit_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),
    RiskBase_med = median(x_10k_base, na.rm = TRUE),
    RiskBase_lo  = quantile(x_10k_base, 0.025, na.rm = TRUE),
    RiskBase_hi  = quantile(x_10k_base, 0.975, na.rm = TRUE),
    RiskAdj_med  = median(x_10k_adj, na.rm = TRUE),
    RiskAdj_lo   = quantile(x_10k_adj, 0.025, na.rm = TRUE),
    RiskAdj_hi   = quantile(x_10k_adj, 0.975, na.rm = TRUE),
    BRRBase_med = median(brr_base, na.rm = TRUE),
    BRRBase_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
    BRRBase_hi  = quantile(brr_base, 0.975, na.rm = TRUE),
    BRRAdj_med  = median(brr_adj, na.rm = TRUE),
    BRRAdj_lo   = quantile(brr_adj, 0.025, na.rm = TRUE),
    BRRAdj_hi   = quantile(brr_adj, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    setting = factor(setting, levels = setting_levels),
    AgeCat  = factor(AgeCat, levels = "65+"),
    mech    = unname(mechanism_abbrev[as.character(VE_label)]),
    Benefit  = fmt_ci(Benefit_med, Benefit_lo, Benefit_hi),
    `Risk (base)` = fmt_ci(RiskBase_med, RiskBase_lo, RiskBase_hi),
    `Risk (adj.)` = fmt_ci(RiskAdj_med, RiskAdj_lo, RiskAdj_hi),
    `BRR (base)` = fmt_ci(BRRBase_med, BRRBase_lo, BRRBase_hi),
    `BRR (adj.)` = fmt_ci(BRRAdj_med, BRRAdj_lo, BRRAdj_hi)
  )

mech_order  <- c("DB", "D+I")
sub_cols    <- c("Benefit", "Risk (base)", "Risk (adj.)", "BRR (base)", "BRR (adj.)")

build_block <- function(m) {
  group_summary %>%
    filter(mech == m) %>%
    select(setting, AgeCat, dplyr::all_of(sub_cols)) %>%
    rename_with(~ paste0(.x, "___", m), dplyr::all_of(sub_cols))
}

blocks <- lapply(mech_order, build_block)
wide <- Reduce(function(a, b) dplyr::left_join(a, b, by = c("setting", "AgeCat")), blocks) %>%
  arrange(setting, AgeCat) %>%
  rename(Transmission = setting, `Age group` = AgeCat)

print(wide, n = 10, width = Inf)

# ---- Assemble two-row-header flextable --------------------------------------
header_labels <- c(Transmission = "Transmission", `Age group` = "Age group")
for (m in mech_order) {
  for (sc in sub_cols) {
    header_labels[paste0(sc, "___", m)] <- sc
  }
}
mech_group_labels <- c("Disease blocking only (DB)", "Disease and infection blocking (D+I)")

ft <- flextable::flextable(wide) %>%
  flextable::set_header_labels(values = as.list(header_labels)) %>%
  flextable::add_header_row(
    values = c("", "", mech_group_labels),
    colwidths = c(1, 1, rep(length(sub_cols), length(mech_order)))
  ) %>%
  flextable::set_caption(caption = "Absolute death benefits and vaccine-attributable risks for outbreak-response immunisation under the reference scenario (65+ only)") %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Transmission", "Age group")) %>%
  flextable::valign(j = c("Transmission", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par("Main Table 1 (Death)", style = "heading 1") %>%
  officer::body_add_par(
    "Absolute death benefits and vaccine-attributable risks for outbreak-response immunisation under the reference scenario",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par(
    paste(
      "DB = disease blocking only; D+I = disease and infection blocking. Reference scenario = outbreak-response",
      "vaccination campaign starting week 2 at 50% coverage. Values are median (95% uncertainty interval) across",
      "1000 posterior draws, per 10,000 vaccinated. 18-64 is omitted: vaccine-attributable death risk is 0 for",
      "this age group, leaving the death BRR undefined. Benefit = deaths averted. Risk (base) = vaccine-",
      "attributable deaths unadjusted for recipient serostatus (identical for DB and D+I -- vaccine safety does",
      "not depend on outbreak intensity or protection mechanism). Risk (adj.) = serostatus-adjusted (seropositive",
      "vaccinees assumed zero vaccine-attributable risk beyond this scenario's baseline); this DOES differ by",
      "mechanism. BRR (base)/(adj.) = benefit-risk ratio (deaths averted / deaths caused) using the",
      "correspondingly-adjusted risk denominator; summarised across draws directly (not the ratio of the two",
      "columns' medians)."
    ),
    style = "Normal"
  )

print(doc, target = "06_Results/Main_Table_1_Death.docx")
message("Saved: 06_Results/Main_Table_1_Death.docx")

writexl::write_xlsx(wide, "06_Results/Main_Table_1_Death.xlsx")
message("Saved: 06_Results/Main_Table_1_Death.xlsx")
