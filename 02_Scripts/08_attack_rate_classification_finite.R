# =============================================================================
# 08_attack_rate_classification_finite.R
#
# Recomputes the state-level transmission-setting classification (Low /
# Moderate / High) using the SAME method as the companion paper (Kang et al.,
# EClinicalMedicine 2025; 90:103690), for consistency across both manuscripts:
#
#   1) Exclude states with insufficient epidemic signal: peak weekly
#      symptomatic cases < 20 per million population, OR total symptomatic
#      cases < 500 (stochastic fluctuation, not an identifiable outbreak).
#   2) Among the remaining states, k-means clustering on (peak cases per
#      million, total cases) is used to examine the empirical clustering
#      pattern (exploratory/validation step -- not the operative rule).
#   3) States are categorised by peak symptomatic cases per million
#      population: <100 = Low, 100-200 = Moderate, >200 = High.
#
# This REPLACES the previously-active case-rate-per-100,000 method (rank-based
# 4 High / 4 Moderate / 3 Low), which was found to misclassify Rio Grande do
# Norte and Sergipe as Moderate (both are Low under the legacy/published
# rule -- confirmed against the published Table S3 classification, see chat
# record 2026-08-19). That case-rate method was itself a replacement for an
# even earlier phi-based method (vaccination-contaminated FOI) -- see git
# history for that superseded version of this script.
#
# Verified: applying this method to observed_2022.RData (raw weekly reported
# cases, all 11 states) + pop_by_state.RData reproduces the published
# classification exactly (11/11 states) -- High: Ceara, Piaui, Paraiba,
# Alagoas; Moderate: Tocantins, Pernambuco; Low: Bahia, Rio Grande do Norte,
# Minas Gerais, Sergipe, Goias.
#
# Output:
#   - 01_Data/setting_key.RData        (Low/Moderate/High membership)
#   - 06_Results/case_rate_by_state.docx (peak/total case-rate table)
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages(library(dplyr))

# ---- 0) Load weekly observed case data (all 11 states, 52 epi weeks) -------
load("c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact/00_Data/0_2_Processed/observed_2022.RData")
load("01_Data/pop_by_state.RData")  # -> pop_by_state (region, tot_pop)

obj_names <- c("observed_ag", "observed_bh", "observed_ce", "observed_go", "observed_mg",
               "observed_pa", "observed_pi", "observed_pn", "observed_rg", "observed_se", "observed_tc")

state_summary <- dplyr::bind_rows(lapply(obj_names, function(nm) {
  df <- get(nm)
  data.frame(
    region      = unique(df$region),
    peak_weekly = max(df$Observed, na.rm = TRUE),
    total_cases = sum(df$Observed, na.rm = TRUE)
  )
}))

# ---- 1) Peak cases per million population -----------------------------------
rate_summary <- state_summary %>%
  dplyr::left_join(pop_by_state, by = "region") %>%
  dplyr::mutate(peak_per_million = peak_weekly / (tot_pop / 1e6))

# ---- 2) Exclusion: insufficient epidemic signal ------------------------------
rate_summary <- rate_summary %>%
  dplyr::mutate(excluded = peak_per_million < 20 | total_cases < 500)

if (any(rate_summary$excluded)) {
  cat("States excluded (insufficient epidemic signal):\n")
  print(rate_summary %>% dplyr::filter(excluded) %>% dplyr::select(region, peak_per_million, total_cases))
}
rate_summary <- rate_summary %>% dplyr::filter(!excluded)

# ---- 3) k-means clustering (exploratory validation of the fixed thresholds) --
set.seed(123)
km_scaled <- scale(rate_summary %>% dplyr::select(peak_per_million, total_cases))
km_fit <- kmeans(km_scaled, centers = 3, nstart = 25)
rate_summary$km_cluster <- km_fit$cluster
cluster_order <- rate_summary %>% dplyr::group_by(km_cluster) %>%
  dplyr::summarise(m = mean(peak_per_million), .groups = "drop") %>% dplyr::arrange(m)
label_map <- setNames(c("Low", "Moderate", "High"), cluster_order$km_cluster)
rate_summary$setting_kmeans <- label_map[as.character(rate_summary$km_cluster)]

# ---- 4) Operative rule: fixed thresholds on peak cases per million ----------
rate_summary <- rate_summary %>%
  dplyr::mutate(setting = dplyr::case_when(
    peak_per_million < 100 ~ "Low",
    peak_per_million <= 200 ~ "Moderate",
    TRUE ~ "High"
  ))

cat("\n=== Transmission-setting classification (peak cases per million; legacy/published rule) ===\n")
print(rate_summary %>% dplyr::arrange(dplyr::desc(peak_per_million)) %>%
        dplyr::select(region, peak_weekly, total_cases, peak_per_million, setting, setting_kmeans))

if (any(rate_summary$setting != rate_summary$setting_kmeans)) {
  cat("\nNOTE: fixed-threshold tier and k-means cluster disagree for:\n")
  print(rate_summary %>% dplyr::filter(setting != setting_kmeans) %>%
          dplyr::select(region, peak_per_million, setting, setting_kmeans))
}

# ---- 5) Save --------------------------------------------------------------
if (file.exists("01_Data/setting_key.RData") && !file.exists("01_Data/setting_key_backup_caserate.RData")) {
  file.copy("01_Data/setting_key.RData", "01_Data/setting_key_backup_caserate.RData")
}

setting_key <- setNames(rate_summary$setting, rate_summary$region)
stopifnot(length(setting_key) == 11)
save(setting_key, file = "01_Data/setting_key.RData")
cat("\nSaved: 01_Data/setting_key.RData (peak-cases-per-million method)\n")
print(setting_key)

# ---- 6) Region-level case-rate table (for supplementary reporting) ---------
rate_table <- rate_summary %>%
  dplyr::arrange(dplyr::desc(peak_per_million)) %>%
  dplyr::transmute(
    Region = region,
    `Peak weekly cases` = peak_weekly,
    `Peak per million pop.` = sprintf("%.1f", peak_per_million),
    `Total cases (2022)` = total_cases,
    Setting = setting
  )

ft <- flextable::flextable(rate_table) %>%
  flextable::set_caption(caption = paste0(
    "State-level transmission-setting classification, following the criteria applied in the companion ",
    "modelling study (Kang et al., EClinicalMedicine 2025;90:103690). States with peak weekly symptomatic ",
    "cases <20 per million population or total symptomatic cases <500 were excluded as insufficient epidemic ",
    "signal (none excluded among these 11 states). Remaining states were categorised by peak symptomatic ",
    "cases per million population: <100 = Low, 100-200 = Moderate, >200 = High."
  )) %>%
  flextable::autofit() %>%
  flextable::theme_vanilla()

doc <- officer::read_docx() %>% flextable::body_add_flextable(ft)
print(doc, target = "06_Results/case_rate_by_state.docx")
cat("Saved: 06_Results/case_rate_by_state.docx\n")
