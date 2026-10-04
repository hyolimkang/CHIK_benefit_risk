# =============================================================================
# 18_rebuild_ori_brr_table.R
#
# Quick rebuild of BRR_table_ori_setting.docx (+ its new supplementary
# companion) WITHOUT re-running 03_brazil_all_draws_ori_v3.R's expensive
# draw-level simulation (SECTIONS 01-07). draw_level_xy_serostatus_finite.RData
# already has every column that section needs (VE_label, setting, x_10k_base/
# adj, brr_base/adj -- verified by inspection), so this script just loads it
# and re-runs the reporting-table assembly (03's SECTION 08B/08C, verbatim)
# on top.
#
# Use this whenever ONLY the table-formatting/layout code changes (e.g. the
# mechanism-pivoted-to-columns compression) -- rerun 03 in full only if the
# underlying draw-level data itself needs to change (new posterior, new rho
# fix, etc.), since only then would draw_level_xy_serostatus_finite.RData
# itself be stale.
#
# Output:
#   - 06_Results/BRR_table_ori_setting_finite.docx (main, compact)
#   - 06_Results/BRR_table_ori_setting_supplementary_finite.docx (full, 36-row)
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(flextable)
  library(officer)
  library(writexl)
})

DATA_MODE   <- "flat"
MODE_SUFFIX <- "_finite"

load("01_Data/draw_level_xy_serostatus_finite.RData")  # -> draw_level_xy_serostatus
stopifnot(all(c("VE_label", "setting", "x_10k_base", "x_10k_adj", "brr_base", "brr_adj") %in%
                colnames(draw_level_xy_serostatus)))

# ---- Helper functions (verbatim from 03_brazil_all_draws_ori_v3.R SECTION 08A) --
make_pr_gt1_wide <- function(ceac_ob, brr_type_filter) {
  ceac_ob %>%
    filter(
      abs(log10(threshold)) < 1e-12,
      brr_type == brr_type_filter,
      RR_seropos == 0
    ) %>%
    transmute(
      Outcome     = outcome,
      Setting     = setting,
      `Age group` = AgeCat,
      VE_col      = VE_label,
      pr_fmt      = sprintf("%.0f%%", 100 * p_accept)
    ) %>%
    pivot_wider(
      names_from  = VE_col,
      values_from = pr_fmt,
      names_glue  = "{VE_col} Pr(BRR>1)"
    )
}

fmt_ci <- function(med, lo, hi) {
  sprintf("%.1f (%.1f–%.1f)", med, lo, hi)
}

make_brr_ceac_outbreak <- function(brr_long,
                                   thresholds = 10^seq(-1, 1, by = 0.02),
                                   group_vars = c("setting","VE_label","AgeCat","outcome")) {
  brr_long %>%
    tidyr::crossing(threshold = thresholds) %>%
    group_by(across(all_of(group_vars)), threshold) %>%
    summarise(
      p_accept = mean(brr > threshold),
      .groups = "drop"
    )
}

# ---- SECTION 08B/08C, verbatim from 03_brazil_all_draws_ori_v3.R ----------
brr_draw_summary_true <- draw_level_xy_serostatus %>%
  dplyr::mutate(
    brr_base = ifelse(is.infinite(brr_base), NA, brr_base),
    brr_adj  = ifelse(is.infinite(brr_adj),  NA, brr_adj),
    setting  = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  dplyr::group_by(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting) %>%
  dplyr::summarise(
    brr_base_med = quantile(brr_base, 0.50,  na.rm = TRUE),
    brr_base_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
    brr_base_hi  = quantile(brr_base, 0.975, na.rm = TRUE),

    brr_adj_med  = quantile(brr_adj,  0.50,  na.rm = TRUE),
    brr_adj_lo   = quantile(brr_adj,  0.025, na.rm = TRUE),
    brr_adj_hi   = quantile(brr_adj,  0.975, na.rm = TRUE),

    av_med = quantile(averted_10k, 0.50,  na.rm = TRUE),
    av_lo  = quantile(averted_10k, 0.025, na.rm = TRUE),
    av_hi  = quantile(averted_10k, 0.975, na.rm = TRUE),

    ca_sae_base_med   = quantile(sae_10k_base,   0.50,  na.rm = TRUE),
    ca_sae_base_lo    = quantile(sae_10k_base,   0.025, na.rm = TRUE),
    ca_sae_base_hi    = quantile(sae_10k_base,   0.975, na.rm = TRUE),

    ca_death_base_med = quantile(death_10k_base, 0.50,  na.rm = TRUE),
    ca_death_base_lo  = quantile(death_10k_base, 0.025, na.rm = TRUE),
    ca_death_base_hi  = quantile(death_10k_base, 0.975, na.rm = TRUE),

    ca_daly_base_med  = quantile(daly_10k_base,  0.50,  na.rm = TRUE),
    ca_daly_base_lo   = quantile(daly_10k_base,  0.025, na.rm = TRUE),
    ca_daly_base_hi   = quantile(daly_10k_base,  0.975, na.rm = TRUE),

    ca_sae_adj_med    = quantile(sae_10k_adj,    0.50,  na.rm = TRUE),
    ca_sae_adj_lo     = quantile(sae_10k_adj,    0.025, na.rm = TRUE),
    ca_sae_adj_hi     = quantile(sae_10k_adj,    0.975, na.rm = TRUE),

    ca_death_adj_med  = quantile(death_10k_adj,  0.50,  na.rm = TRUE),
    ca_death_adj_lo   = quantile(death_10k_adj,  0.025, na.rm = TRUE),
    ca_death_adj_hi   = quantile(death_10k_adj,  0.975, na.rm = TRUE),

    ca_daly_adj_med   = quantile(daly_10k_adj,   0.50,  na.rm = TRUE),
    ca_daly_adj_lo    = quantile(daly_10k_adj,   0.025, na.rm = TRUE),
    ca_daly_adj_hi    = quantile(daly_10k_adj,   0.975, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    scenario  = Scenario,
    age_group = AgeCat
  ) %>%
  dplyr::arrange(outcome, Scenario, setting, AgeCat, VE_label, RR_seropos)

brr_draw_summary_true_filtered <- brr_draw_summary_true %>%
  dplyr::filter(
    AgeCat %in% c("18-64", "65+"),
    RR_seropos == 0.0
  ) %>%
  dplyr::select(-scenario, -age_group)

ca_summary <- brr_draw_summary_true_filtered %>%
  filter(outcome == "DALY") %>%
  select(Scenario, AgeCat, VE_label, RR_seropos, setting,
         starts_with("ca_"))

brr_long_serostatus <- draw_level_xy_serostatus %>%
  mutate(
    outcome = recode(outcome, "sae" = "SAE", "death" = "Death", "daly" = "DALY"),
    setting = factor(setting, levels = c("Low", "Moderate", "High")),
    AgeCat  = factor(AgeCat,  levels = c("18-64", "65+"))
  ) %>%
  filter(
    is.finite(brr_base), brr_base > 0,
    is.finite(brr_adj),  brr_adj  > 0
  ) %>%
  transmute(
    Region, setting, VE_label, AgeCat, outcome, RR_seropos,
    brr_base,
    brr_adj
  )

brr_long_serostatus_long <- brr_long_serostatus %>%
  filter(!is.na(AgeCat)) %>%
  pivot_longer(
    cols      = c(brr_base, brr_adj),
    names_to  = "brr_type",
    values_to = "brr"
  ) %>%
  filter(is.finite(brr), brr > 0)

brr_range <- brr_long_serostatus_long %>%
  filter(RR_seropos == 0) %>%
  summarise(
    min_brr = min(brr, na.rm = TRUE),
    max_brr = max(brr, na.rm = TRUE)
  )

brr_min <- brr_range$min_brr
brr_max <- brr_range$max_brr

thresholds_auto <- 10^seq(
  floor(log10(brr_min)),
  ceiling(log10(brr_max)),
  by = 0.02
)
ceac_ob <- make_brr_ceac_outbreak(
  brr_long_serostatus_long,
  group_vars = c("setting", "VE_label", "AgeCat", "outcome", "RR_seropos", "brr_type")
)

pr_gt1_wide_setting_base <- make_pr_gt1_wide(ceac_ob, "brr_base")
pr_gt1_wide_setting_adj  <- make_pr_gt1_wide(ceac_ob, "brr_adj")

ca_summary <- brr_draw_summary_true_filtered %>%
  filter(RR_seropos == 0, outcome == "DALY") %>%
  select(Scenario, AgeCat, VE_label, setting,
         ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi,
         ca_death_base_med, ca_death_base_lo, ca_death_base_hi,
         ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi,
         ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi,
         ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi,
         ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)

brr_draw_for_wide <- brr_draw_summary_true_filtered %>%
  filter(RR_seropos == 0) %>%
  select(outcome, Scenario, AgeCat, VE_label, setting,
         av_med,       av_lo,       av_hi,
         brr_base_med, brr_base_lo, brr_base_hi,
         brr_adj_med,  brr_adj_lo,  brr_adj_hi) %>%
  left_join(ca_summary, by = c("Scenario", "AgeCat", "VE_label", "setting"))

brr_table_wide_serostatus <- brr_draw_for_wide %>%
  mutate(
    setting = factor(setting, levels = c("Low", "Moderate", "High")),

    av_formatted            = fmt_ci(av_med,       av_lo,       av_hi),
    brr_base_formatted      = fmt_ci(brr_base_med, brr_base_lo, brr_base_hi),
    brr_adj_formatted       = fmt_ci(brr_adj_med,  brr_adj_lo,  brr_adj_hi),
    ca_sae_base_formatted   = fmt_ci(ca_sae_base_med,   ca_sae_base_lo,   ca_sae_base_hi),
    ca_death_base_formatted = fmt_ci(ca_death_base_med, ca_death_base_lo, ca_death_base_hi),
    ca_daly_base_formatted  = fmt_ci(ca_daly_base_med,  ca_daly_base_lo,  ca_daly_base_hi),
    ca_sae_adj_formatted    = fmt_ci(ca_sae_adj_med,    ca_sae_adj_lo,    ca_sae_adj_hi),
    ca_death_adj_formatted  = fmt_ci(ca_death_adj_med,  ca_death_adj_lo,  ca_death_adj_hi),
    ca_daly_adj_formatted   = fmt_ci(ca_daly_adj_med,   ca_daly_adj_lo,   ca_daly_adj_hi)
  ) %>%
  select(outcome, Scenario, AgeCat, setting, VE_label,
         av_formatted, brr_base_formatted, brr_adj_formatted,
         ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
         ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted) %>%
  pivot_wider(
    names_from  = VE_label,
    values_from = c(av_formatted, brr_base_formatted, brr_adj_formatted,
                    ca_sae_base_formatted, ca_death_base_formatted, ca_daly_base_formatted,
                    ca_sae_adj_formatted,  ca_death_adj_formatted,  ca_daly_adj_formatted),
    names_glue  = "{VE_label}_{.value}"
  ) %>%
  arrange(outcome, Scenario, setting, AgeCat) %>%
  dplyr::rename(
    Outcome     = outcome,
    `Age group` = AgeCat,
    Setting     = setting
  )

brr_table_wide_serostatus <- brr_table_wide_serostatus %>%
  dplyr::rename(
    `DB (Averted)`           = `Disease blocking only_av_formatted`,
    `DB (BRR base)`          = `Disease blocking only_brr_base_formatted`,
    `DB (BRR adj)`           = `Disease blocking only_brr_adj_formatted`,
    `DB (SAE caused base)`   = `Disease blocking only_ca_sae_base_formatted`,
    `DB (SAE caused adj)`    = `Disease blocking only_ca_sae_adj_formatted`,
    `DB (Death caused base)` = `Disease blocking only_ca_death_base_formatted`,
    `DB (Death caused adj)`  = `Disease blocking only_ca_death_adj_formatted`,
    `DB (DALY caused base)`  = `Disease blocking only_ca_daly_base_formatted`,
    `DB (DALY caused adj)`   = `Disease blocking only_ca_daly_adj_formatted`,
    `DIB (Averted)`           = `Disease and infection blocking_av_formatted`,
    `DIB (BRR base)`          = `Disease and infection blocking_brr_base_formatted`,
    `DIB (BRR adj)`           = `Disease and infection blocking_brr_adj_formatted`,
    `DIB (SAE caused base)`   = `Disease and infection blocking_ca_sae_base_formatted`,
    `DIB (SAE caused adj)`    = `Disease and infection blocking_ca_sae_adj_formatted`,
    `DIB (Death caused base)` = `Disease and infection blocking_ca_death_base_formatted`,
    `DIB (Death caused adj)`  = `Disease and infection blocking_ca_death_adj_formatted`,
    `DIB (DALY caused base)`  = `Disease and infection blocking_ca_daly_base_formatted`,
    `DIB (DALY caused adj)`   = `Disease and infection blocking_ca_daly_adj_formatted`
  )

pr_gt1_wide_serostatus_base <- pr_gt1_wide_setting_base %>%
  dplyr::rename_with(
    ~ gsub("Pr(BRR>1)", "Pr(BRR_base>1)", .x, fixed = TRUE),
    .cols = -c(Outcome, Setting, `Age group`)
  )

pr_gt1_wide_serostatus_adj <- pr_gt1_wide_setting_adj %>%
  dplyr::rename_with(
    ~ gsub("Pr(BRR>1)", "Pr(BRR_adj>1)", .x, fixed = TRUE),
    .cols = -c(Outcome, Setting, `Age group`)
  )

brr_table_wide_serostatus2 <- brr_table_wide_serostatus %>%
  left_join(pr_gt1_wide_serostatus_base, by = c("Outcome", "Setting", "Age group")) %>%
  left_join(pr_gt1_wide_serostatus_adj,  by = c("Outcome", "Setting", "Age group"))

brr_table_wide_serostatus2 <- brr_table_wide_serostatus2 %>%
  relocate(`DB (Averted)`,                                  .after = `Age group`) %>%
  relocate(`DB (BRR base)`,                                 .after = `DB (Averted)`) %>%
  relocate(`DB (BRR adj)`,                                  .after = `DB (BRR base)`) %>%
  relocate(`Disease blocking only Pr(BRR_base>1)`,          .after = `DB (BRR adj)`) %>%
  relocate(`Disease blocking only Pr(BRR_adj>1)`,           .after = `Disease blocking only Pr(BRR_base>1)`) %>%
  relocate(`DB (SAE caused base)`,                          .after = `Disease blocking only Pr(BRR_adj>1)`) %>%
  relocate(`DB (SAE caused adj)`,                           .after = `DB (SAE caused base)`) %>%
  relocate(`DB (Death caused base)`,                        .after = `DB (SAE caused adj)`) %>%
  relocate(`DB (Death caused adj)`,                         .after = `DB (Death caused base)`) %>%
  relocate(`DB (DALY caused base)`,                         .after = `DB (Death caused adj)`) %>%
  relocate(`DB (DALY caused adj)`,                          .after = `DB (DALY caused base)`) %>%
  relocate(`DIB (Averted)`,                                 .after = `DB (DALY caused adj)`) %>%
  relocate(`DIB (BRR base)`,                                .after = `DIB (Averted)`) %>%
  relocate(`DIB (BRR adj)`,                                 .after = `DIB (BRR base)`) %>%
  relocate(`Disease and infection blocking Pr(BRR_base>1)`, .after = `DIB (BRR adj)`) %>%
  relocate(`Disease and infection blocking Pr(BRR_adj>1)`,  .after = `Disease and infection blocking Pr(BRR_base>1)`) %>%
  relocate(`DIB (SAE caused base)`,                         .after = `Disease and infection blocking Pr(BRR_adj>1)`) %>%
  relocate(`DIB (SAE caused adj)`,                          .after = `DIB (SAE caused base)`) %>%
  relocate(`DIB (Death caused base)`,                       .after = `DIB (SAE caused adj)`) %>%
  relocate(`DIB (Death caused adj)`,                        .after = `DIB (Death caused base)`) %>%
  relocate(`DIB (DALY caused base)`,                        .after = `DIB (Death caused adj)`) %>%
  relocate(`DIB (DALY caused adj)`,                         .after = `DIB (DALY caused base)`) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      ~ tidyr::replace_na(as.character(.x), "beneficial")
    )
  )

part1 <- brr_table_wide_serostatus2 %>%
  dplyr::mutate(
    Risk_base = dplyr::case_when(
      Outcome == "SAE"   ~ `DB (SAE caused base)`,
      Outcome == "Death" ~ `DB (Death caused base)`,
      Outcome == "DALY"  ~ `DB (DALY caused base)`,
      TRUE ~ NA_character_
    ),
    Risk_adj = dplyr::case_when(
      Outcome == "SAE"   ~ `DB (SAE caused adj)`,
      Outcome == "Death" ~ `DB (Death caused adj)`,
      Outcome == "DALY"  ~ `DB (DALY caused adj)`,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `DB (Averted)`,
    Risk_base,
    Risk_adj,
    BRR_base = `DB (BRR base)`,
    BRR_adj  = `DB (BRR adj)`,
    prob_base = `Disease blocking only Pr(BRR_base>1)`,
    prob_adj  = `Disease blocking only Pr(BRR_adj>1)`
  ) %>%
  dplyr::mutate(mechanism = "Disease blocking only")

part2 <- brr_table_wide_serostatus2 %>%
  dplyr::mutate(
    Risk_base = dplyr::case_when(
      Outcome == "SAE"   ~ `DIB (SAE caused base)`,
      Outcome == "Death" ~ `DIB (Death caused base)`,
      Outcome == "DALY"  ~ `DIB (DALY caused base)`,
      TRUE ~ NA_character_
    ),
    Risk_adj = dplyr::case_when(
      Outcome == "SAE"   ~ `DIB (SAE caused adj)`,
      Outcome == "Death" ~ `DIB (Death caused adj)`,
      Outcome == "DALY"  ~ `DIB (DALY caused adj)`,
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(
    Outcome, Setting, `Age group`,
    Benefit = `DIB (Averted)`,
    Risk_base,
    Risk_adj,
    BRR_base = `DIB (BRR base)`,
    BRR_adj  = `DIB (BRR adj)`,
    prob_base = `Disease and infection blocking Pr(BRR_base>1)`,
    prob_adj  = `Disease and infection blocking Pr(BRR_adj>1)`
  ) %>%
  dplyr::mutate(mechanism = "Disease and infection blocking")

brr_table_final_long <- dplyr::bind_rows(part1, part2) %>%
  dplyr::mutate(
    Setting = factor(Setting, levels = c("High", "Moderate", "Low"))
  ) %>%
  dplyr::arrange(Outcome, Setting, `Age group`) %>%
  dplyr::select(
    Outcome, Setting, `Age group`, mechanism,
    Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj, prob_base, prob_adj
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(Benefit, Risk_base, Risk_adj, BRR_base, BRR_adj),
      ~ {
        x <- gsub(" \\(", "\n(", .x, fixed = FALSE)
        dplyr::case_when(
          x %in% c(
            "NA (NA-NA)", "NA (NA–NA)",
            "NA\n(NA-NA)", "NA\n(NA–NA)",
            "NA (NA NA)", "NA\n(NA NA)"
          ) ~ "beneficial",
          TRUE ~ x
        )
      }
    )
  )

# ---- Main table: mechanism pivoted to columns, not rows -------------------
# Pr(BRR>1) dropped from the main table (too dense) -- still in the
# supplementary long table (brr_table_final_long) above for full detail.
metric_cols <- c("Benefit", "Risk_base", "Risk_adj", "BRR_base", "BRR_adj")

brr_table_final_wide <- part1 %>%
  dplyr::select(Outcome, Setting, `Age group`, dplyr::all_of(metric_cols)) %>%
  dplyr::rename_with(~ paste0(.x, "_DB"), dplyr::all_of(metric_cols)) %>%
  dplyr::left_join(
    part2 %>%
      dplyr::select(Outcome, Setting, `Age group`, dplyr::all_of(metric_cols)) %>%
      dplyr::rename_with(~ paste0(.x, "_DIB"), dplyr::all_of(metric_cols)),
    by = c("Outcome", "Setting", "Age group")
  ) %>%
  dplyr::mutate(Setting = factor(Setting, levels = c("High", "Moderate", "Low"))) %>%
  dplyr::arrange(Outcome, Setting, `Age group`) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::ends_with(c("_DB", "_DIB")) & dplyr::matches("Benefit|Risk|BRR"),
      ~ {
        x <- gsub(" \\(", "\n(", .x, fixed = FALSE)
        dplyr::case_when(
          x %in% c("NA (NA-NA)", "NA (NA–NA)", "NA\n(NA-NA)", "NA\n(NA–NA)",
                    "NA (NA NA)", "NA\n(NA NA)") ~ "beneficial",
          TRUE ~ x
        )
      }
    )
  )

metric_labels <- c(
  Benefit   = "Benefit\n(per 10,000)",
  Risk_base = "Risk (base)\n(per 10,000)",
  Risk_adj  = "Risk (adj.)\n(per 10,000)",
  BRR_base  = "BRR (base)",
  BRR_adj   = "BRR (adj.)"
)
header_labels_wide <- c(
  Setting = "Setting", `Age group` = "Age group",
  setNames(metric_labels, paste0(names(metric_labels), "_DB")),
  setNames(metric_labels, paste0(names(metric_labels), "_DIB"))
)

# One table PER outcome (Setting x Age group, 6 rows each) instead of one
# combined 18-row table -- DALY is the main-text table; Death and SAE get
# their own files. Outcome column dropped within each file since it's now
# constant (redundant).
build_outcome_table <- function(outcome_filter, out_file, caption_text) {
  df <- brr_table_final_wide %>%
    dplyr::filter(Outcome == outcome_filter) %>%
    dplyr::select(-Outcome)

  ft <- flextable::flextable(df) %>%
    flextable::set_header_labels(values = as.list(header_labels_wide)) %>%
    flextable::add_header_row(
      values = c("", "", "Disease blocking only", "Disease and infection blocking"),
      colwidths = c(1, 1, length(metric_cols), length(metric_cols))
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
  message("Saved: ", out_file)

  xlsx_file <- sub("\\.docx$", ".xlsx", out_file)
  writexl::write_xlsx(df, xlsx_file)
  message("Saved: ", xlsx_file)
}

# Main text table: DALY only.
build_outcome_table(
  "DALY", paste0("06_Results/BRR_table_ori_setting", MODE_SUFFIX, ".docx"),
  "Benefit-Risk Ratio (BRR) by Setting and Age group -- DALY"
)
# Companion single-outcome tables (Death, SAE) -- own files, not in the main table.
build_outcome_table(
  "Death", paste0("06_Results/BRR_table_ori_setting_Death", MODE_SUFFIX, ".docx"),
  "Benefit-Risk Ratio (BRR) by Setting and Age group -- Death"
)
build_outcome_table(
  "SAE", paste0("06_Results/BRR_table_ori_setting_SAE", MODE_SUFFIX, ".docx"),
  "Benefit-Risk Ratio (BRR) by Setting and Age group -- SAE"
)

# ---- Supplementary table: full long format (mechanism as rows, 36 rows) ---
ft_brr_long <- flextable::flextable(brr_table_final_long) %>%
  flextable::set_header_labels(
    Outcome     = "Outcome",
    Setting     = "Setting",
    `Age group` = "Age group",
    mechanism   = "Vaccine protection\nmechanism",
    Benefit     = "Benefit:\nOutcomes averted\n(per 10,000)",
    Risk_base   = "Risk (base):\nOutcomes attributable\n(per 10,000)",
    Risk_adj    = "Risk (adjusted):\nOutcomes attributable\n(per 10,000)",
    BRR_base    = "BRR (base):\n(Prevented per 1 caused)",
    BRR_adj     = "BRR (adjusted):\n(Prevented per 1 caused)",
    prob_base   = "Probability\n(BRR_base > 1)\n(%)",
    prob_adj    = "Probability\n(BRR_adj > 1)\n(%)"
  ) %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = 1:4, align = "left", part = "all") %>%
  flextable::merge_v(j = c("Outcome", "Setting", "Age group")) %>%
  flextable::valign(j = c("Outcome", "Setting", "Age group"), valign = "top") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc_supp <- officer::read_docx() %>%
  officer::body_add_par(
    "Supplementary Table: Benefit-Risk Ratio (BRR) by Outcome, Setting, Age group, and Mechanism (full breakdown)",
    style = "heading 2"
  ) %>%
  flextable::body_add_flextable(ft_brr_long)

print(doc_supp, target = paste0("06_Results/BRR_table_ori_setting_supplementary", MODE_SUFFIX, ".docx"))
message("Saved: 06_Results/BRR_table_ori_setting_supplementary", MODE_SUFFIX, ".docx")

writexl::write_xlsx(brr_table_final_long, paste0("06_Results/BRR_table_ori_setting_supplementary", MODE_SUFFIX, ".xlsx"))
message("Saved: 06_Results/BRR_table_ori_setting_supplementary", MODE_SUFFIX, ".xlsx")
