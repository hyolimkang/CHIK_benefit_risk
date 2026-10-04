# =============================================================================
# 24_main_table1_sae.R
#
# Main Table 1 (SAE companion): Absolute SAE benefits and vaccine-attributable
# risks for outbreak-response immunisation under the reference scenario
# (week 2, 50% coverage). Same structure as 19_main_table1_daly.R (rows =
# Transmission x Age group, columns = Mechanism x {Benefit, Risk (base),
# Risk (adj.), BRR (base), BRR (adj.)}) -- only the outcome filter and output
# filenames differ. Unlike Death, SAE risk is nonzero for both age groups, so
# no row exclusion is needed here (contrast with 25_main_table1_death.R).
#
# Output: 06_Results/Main_Table_1_SAE.docx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(writexl)
})

age_levels     <- c("18-64", "65+")
setting_levels <- c("Low", "Moderate", "High")
mechanism_abbrev <- c("Disease blocking only" = "DB", "Disease and infection blocking" = "D+I")

fmt_ci <- function(med, lo, hi) sprintf("%.1f (%.1f\u2013%.1f)", med, lo, hi)

load("01_Data/draw_level_xy_serostatus_finite.RData")  # -> draw_level_xy_serostatus

d <- draw_level_xy_serostatus %>%
  filter(outcome == "SAE", RR_seropos == 0, AgeCat %in% age_levels)

# ---- Setting x Age group x Mechanism: Benefit, Risk (base/adj), BRR (base/adj) ----
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
    AgeCat  = factor(AgeCat, levels = age_levels),
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
  flextable::set_caption(caption = "Absolute SAE benefits and vaccine-attributable risks for outbreak-response immunisation under the reference scenario") %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:2, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Transmission", "Age group")) %>%
  flextable::valign(j = c("Transmission", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par("Main Table 1 (SAE)", style = "heading 1") %>%
  officer::body_add_par(
    "Absolute SAE benefits and vaccine-attributable risks for outbreak-response immunisation under the reference scenario",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par(
    paste(
      "DB = disease blocking only; D+I = disease and infection blocking. Reference scenario = outbreak-response",
      "vaccination campaign starting week 2 at 50% coverage. Values are median (95% uncertainty interval) across",
      "1000 posterior draws, per 10,000 vaccinated. Benefit = SAEs averted. Risk (base) = vaccine-attributable",
      "SAEs unadjusted for recipient serostatus (identical for DB and D+I within an age group -- vaccine safety",
      "does not depend on outbreak intensity or protection mechanism). Risk (adj.) = serostatus-adjusted",
      "(seropositive vaccinees assumed zero vaccine-attributable risk beyond this scenario's baseline); this DOES",
      "differ by mechanism. BRR (base)/(adj.) = benefit-risk ratio (SAEs averted / SAEs caused) using the",
      "correspondingly-adjusted risk denominator; summarised across draws directly (not the ratio of the two",
      "columns' medians)."
    ),
    style = "Normal"
  )

print(doc, target = "06_Results/Main_Table_1_SAE.docx")
message("Saved: 06_Results/Main_Table_1_SAE.docx")

writexl::write_xlsx(wide, "06_Results/Main_Table_1_SAE.xlsx")
message("Saved: 06_Results/Main_Table_1_SAE.xlsx")
