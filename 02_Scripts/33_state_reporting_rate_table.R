# =============================================================================
# 33_state_reporting_rate_table.R
#
# Results-section background table: per-state PREDICTED REPORTING RATE (rho),
# the Stan-fitted parameter that converts observed/reported case counts into
# the true (unreported-inclusive) incidence used to derive the cumulative
# infection attack rate in Table_S_state_outbreak_summary /
# 16_figure1_epi_map_ar.R panel B. rho comes directly from each state's Stan
# posterior (posterior_<abbrev>$rho, 1000 draws) -- see lhs_samples_sir_finite.R's
# note that rho/gamma/sigma/beta_observed/I0 are all taken directly from the
# fitted posterior (not independently resampled), so this rho is the exact
# quantity used everywhere downstream, not a separate estimate.
#
# Output: 06_Results/Table_S_state_reporting_rate.docx / .xlsx
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(flextable)
  library(officer)
  library(writexl)
})

load("01_Data/setting_key.RData")  # -> setting_key
setting_levels <- c("Low", "Moderate", "High")
setting_rank <- c(Low = 1, Moderate = 2, High = 3)

state_abbrev_map <- c(
  ce = "Ceará", bh = "Bahia", pa = "Paraíba", pn = "Pernambuco",
  rg = "Rio Grande do Norte", pi = "Piauí", ag = "Alagoas",
  tc = "Tocantins", mg = "Minas Gerais", se = "Sergipe", go = "Goiás"
)

rho_cache_file <- "01_Data/state_rho_summary_finite.RData"
if (!file.exists(rho_cache_file)) {
  message("Cache not found -- loading 01_Data/posterior_finite_all.RData (~1.3GB, this is slow)...")
  load("01_Data/posterior_finite_all.RData")
  rho_summary <- dplyr::bind_rows(lapply(names(state_abbrev_map), function(ab) {
    post <- get(paste0("posterior_", ab))
    data.frame(
      region = state_abbrev_map[[ab]],
      rho_med = median(post$rho, na.rm = TRUE),
      rho_lo  = quantile(post$rho, 0.025, na.rm = TRUE, names = FALSE),
      rho_hi  = quantile(post$rho, 0.975, na.rm = TRUE, names = FALSE)
    )
  }))
  save(rho_summary, file = rho_cache_file)
  message("Saved: ", rho_cache_file)
} else {
  load(rho_cache_file)  # -> rho_summary
}

# Reuse the same cached AR summary as Table_S_state_outbreak_summary / Figure 1,
# so this table's row order and the AR it feeds into line up exactly.
load("01_Data/state_ar_hazard_summary_finite.RData")  # -> ar_summary

tbl <- rho_summary %>%
  left_join(ar_summary %>% select(region, ar_med, ar_lo, ar_hi), by = "region") %>%
  mutate(setting = factor(setting_key[region], levels = setting_levels)) %>%
  arrange(setting_rank[as.character(setting)], ar_med) %>%
  transmute(
    Setting = setting,
    State = region,
    `Predicted reporting rate, % (95% UI)` = sprintf(
      "%.1f (%.1f\u2013%.1f)", 100 * rho_med, 100 * rho_lo, 100 * rho_hi
    ),
    `Resulting cumulative infection attack rate, % (95% UI)` = sprintf(
      "%.1f (%.1f\u2013%.1f)", 100 * ar_med, 100 * ar_lo, 100 * ar_hi
    )
  )

print(as.data.frame(tbl))

ft <- flextable::flextable(tbl) %>%
  flextable::merge_v(j = "Setting") %>%
  flextable::valign(j = "Setting", valign = "top") %>%
  flextable::set_caption(caption = "Predicted reporting rate by state, and the resulting inferred attack rate") %>%
  flextable::theme_booktabs() %>%
  flextable::bold(part = "header") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Setting", "State"), align = "left", part = "all") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par("Predicted reporting rate by state", style = "heading 2") %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par(
    paste(
      "Reporting rate (rho) is a state-specific parameter fitted jointly with the transmission model",
      "in Stan, representing the fraction of true infections captured by the reported case surveillance",
      "data; median and 95% uncertainty interval across 1000 posterior draws. True (reporting-adjusted)",
      "incidence = reported cases / rho, and the cumulative infection attack rate reported here and in",
      "Figure 1 / Table_S_state_outbreak_summary is derived from this reporting-adjusted incidence",
      "(1 - exp(-cumulative fitted weekly hazard)). States are ordered by transmission setting, then by",
      "attack rate within each setting, matching Figure 1 and Table_S_state_outbreak_summary."
    ),
    style = "Normal"
  )

print(doc, target = "06_Results/Table_S_state_reporting_rate.docx")
message("Saved: 06_Results/Table_S_state_reporting_rate.docx")

writexl::write_xlsx(tbl, "06_Results/Table_S_state_reporting_rate.xlsx")
message("Saved: 06_Results/Table_S_state_reporting_rate.xlsx")
