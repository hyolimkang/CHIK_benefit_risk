# =============================================================================
# 32_figure1_state_summary_table.R
#
# Results-section background table: per-state summary of the reconstructed
# 2022 chikungunya outbreak underlying Figure 1 (16_figure1_epi_map_ar.R) --
# transmission setting classification, population, peak weekly case rate
# (the classification metric itself), and posterior cumulative infection
# attack rate (median + 95% UI, from the Stan fit's lambda draws). Same data,
# same cache, same row order as Figure 1 panels B/C -- this is just the
# numbers behind those two panels laid out as a table for the manuscript
# Results section (background paragraph before the benefit-risk results).
#
# Output: 06_Results/Table_S_state_outbreak_summary.docx / .xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(writexl)
})

load("01_Data/setting_key.RData")  # -> setting_key

state_abbrev_map <- c(
  ce = "Ceará", bh = "Bahia", pa = "Paraíba", pn = "Pernambuco",
  rg = "Rio Grande do Norte", pi = "Piauí", ag = "Alagoas",
  tc = "Tocantins", mg = "Minas Gerais", se = "Sergipe", go = "Goiás"
)

# ---- Peak weekly cases per million (classification metric) -----------------
load("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_2_Processed/observed_2022.RData")
load("01_Data/pop_by_state.RData")  # -> pop_by_state (region, tot_pop)

obj_names <- c("observed_ag", "observed_bh", "observed_ce", "observed_go", "observed_mg",
               "observed_pa", "observed_pi", "observed_pn", "observed_rg", "observed_se", "observed_tc")
peak_rate_summary <- dplyr::bind_rows(lapply(obj_names, function(nm) {
  df <- get(nm)
  data.frame(region = unique(df$region), peak_weekly = max(df$Observed, na.rm = TRUE))
})) %>%
  dplyr::left_join(pop_by_state, by = "region") %>%
  dplyr::mutate(peak_per_million = peak_weekly / (tot_pop / 1e6))

# ---- Posterior cumulative attack rate (same cache as Figure 1) -------------
ar_cache_file <- "01_Data/state_ar_hazard_summary_finite.RData"
if (!file.exists(ar_cache_file)) {
  message("Cache not found -- loading 01_Data/posterior_finite_all.RData (~1.3GB, this is slow)...")
  load("01_Data/posterior_finite_all.RData")
  ar_summary <- dplyr::bind_rows(lapply(names(state_abbrev_map), function(ab) {
    post <- get(paste0("posterior_", ab))
    cum_hazard <- rowSums(post$lambda)
    ar_draw <- 1 - exp(-cum_hazard)
    data.frame(region = state_abbrev_map[[ab]],
               ar_med = median(ar_draw),
               ar_lo  = quantile(ar_draw, 0.025, names = FALSE),
               ar_hi  = quantile(ar_draw, 0.975, names = FALSE))
  }))
  save(ar_summary, file = ar_cache_file)
} else {
  load(ar_cache_file)  # -> ar_summary
}

# ---- Assemble table, same sort as Figure 1 (tier, then AR within tier) -----
setting_levels <- c("Low", "Moderate", "High")
setting_rank <- c(Low = 1, Moderate = 2, High = 3)

tbl <- ar_summary %>%
  left_join(peak_rate_summary, by = "region") %>%
  mutate(
    setting = factor(setting_key[region], levels = setting_levels)
  ) %>%
  arrange(setting_rank[as.character(setting)], ar_med) %>%
  transmute(
    Setting = setting,
    State = region,
    `Population` = scales::comma(round(tot_pop)),
    `Peak weekly cases per million` = sprintf("%.0f", peak_per_million),
    `Cumulative infection attack rate, % (95% UI)` = sprintf(
      "%.1f (%.1f\u2013%.1f)", 100 * ar_med, 100 * ar_lo, 100 * ar_hi
    )
  )

print(as.data.frame(tbl))

# ---- Word table --------------------------------------------------------------
ft <- flextable::flextable(tbl) %>%
  flextable::merge_v(j = "Setting") %>%
  flextable::valign(j = "Setting", valign = "top") %>%
  flextable::set_caption(caption = "Reconstructed 2022 chikungunya outbreak, by state") %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Setting", "State"), align = "left", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par("Reconstructed chikungunya outbreak: state-level summary", style = "heading 2") %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par(
    paste(
      "Transmission setting classification, population, peak weekly reported case rate (the metric",
      "used to define the Low/Moderate/High tiers; thresholds at 100 and 200 cases per million per",
      "week), and posterior cumulative infection attack rate (1 - exp(-sum of the fitted weekly hazard),",
      "median and 95% uncertainty interval across posterior draws) for each of the 11 states included",
      "in this analysis. States are ordered by transmission setting, then by attack rate within each",
      "setting; this ordering, and the underlying data, match Figure 1 panels B and C."
    ),
    style = "Normal"
  )

print(doc, target = "06_Results/Table_S_state_outbreak_summary.docx")
message("Saved: 06_Results/Table_S_state_outbreak_summary.docx")

writexl::write_xlsx(tbl, "06_Results/Table_S_state_outbreak_summary.xlsx")
message("Saved: 06_Results/Table_S_state_outbreak_summary.xlsx")
