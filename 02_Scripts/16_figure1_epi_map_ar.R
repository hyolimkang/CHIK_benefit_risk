# =============================================================================
# 16_figure1_epi_map_ar.R
#
# Figure 1. Epidemiological reconstruction and analytical framework.
#   Panel A: Map of the 11 Brazilian states included in this analysis,
#            coloured by transmission setting (Low/Moderate/High, from
#            01_Data/setting_key.RData -- peak-cases-per-million classification,
#            08_attack_rate_classification_finite.R), zoomed to the focal
#            region (not full Brazil) so the other 16 states are visible only
#            as unshaded context, not the analytical footprint.
#   Panel B: Posterior cumulative infection attack rate by state,
#            1 - exp(-sum_t(lambda_t)), median + 95% UI, points coloured by
#            the same Low/Moderate/High setting -- sorted, so the reader can
#            see the classification is consistent with the underlying
#            posterior hazard, not just the case-rate proxy used to derive it.
#
# Data sources:
#   - 01_Data/setting_key.RData (Low/Moderate/High membership, 11 states)
#   - 01_Data/posterior_finite_all.RData (state-level Stan posteriors,
#     posterior_<abbrev>$lambda: 1000 draws x 54 weeks daily/weekly hazard)
#     -- LARGE (~1.3GB); the cumulative-hazard summary is cached to
#     01_Data/state_ar_hazard_summary_finite.RData after first run so this
#     script doesn't reload the full posterior on every rerun.
#   - Brazil state boundaries: rnaturalearth::ne_states() (geobr's metadata
#     endpoint was failing at the time this was written -- see comment below
#     if switching back).
#
# Output: 06_Results/figure1_epi_map_ar.png
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(patchwork)
  library(scales)
})

load("01_Data/setting_key.RData")  # -> setting_key (region -> Low/Moderate/High)

pal_setting <- c(Low = "#2A9D8F", Moderate = "#F4A261", High = "#C44E52")

theme_fig1 <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
}

# ---- A) State abbreviation -> full name map (matches posterior_finite_all.RData
# object names and setting_key's names exactly) ------------------------------
state_abbrev_map <- c(
  ce = "Ceará", bh = "Bahia", pa = "Paraíba", pn = "Pernambuco",
  rg = "Rio Grande do Norte", pi = "Piauí", ag = "Alagoas",
  tc = "Tocantins", mg = "Minas Gerais", se = "Sergipe", go = "Goiás"
)
stopifnot(setequal(unname(state_abbrev_map), names(setting_key)))

# ---- A2) Peak weekly cases per million -- the actual classification metric -
# (same logic as 08_attack_rate_classification_finite.R). Shown alongside the
# posterior cumulative AR (panel B) because they measure different things --
# peak intensity vs whole-season cumulative burden -- and don't always rank
# the same way (e.g. Sergipe vs Pernambuco).
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

# ---- B) Posterior cumulative attack rate per state (cached) ----------------
ar_cache_file <- "01_Data/state_ar_hazard_summary_finite.RData"

if (!file.exists(ar_cache_file)) {
  message("Cache not found -- loading 01_Data/posterior_finite_all.RData (~1.3GB, this is slow)...")
  load("01_Data/posterior_finite_all.RData")

  ar_summary <- dplyr::bind_rows(lapply(names(state_abbrev_map), function(ab) {
    post <- get(paste0("posterior_", ab))
    cum_hazard <- rowSums(post$lambda)          # per-draw sum_t(lambda_t)
    ar_draw <- 1 - exp(-cum_hazard)             # per-draw cumulative attack rate
    data.frame(
      region = state_abbrev_map[[ab]],
      ar_med = median(ar_draw),
      ar_lo  = quantile(ar_draw, 0.025, names = FALSE),
      ar_hi  = quantile(ar_draw, 0.975, names = FALSE)
    )
  }))

  save(ar_summary, file = ar_cache_file)
  message("Saved: ", ar_cache_file)
} else {
  load(ar_cache_file)  # -> ar_summary
}


# Sort by classification tier first (Low -> Moderate -> High, matching panel
# A's legend order), then by AR within each tier -- NOT by raw AR value alone.
# The peak-cases-per-million classification (setting_key) and the posterior
# cumulative attack rate (1-exp(-sum(lambda))) measure different things --
# peak intensity vs whole-season cumulative burden -- so their rankings don't
# perfectly agree (e.g. Sergipe has a lower peak but a longer/more sustained
# epidemic than Pernambuco, giving it a higher cumulative AR despite being
# classified Low vs Pernambuco's Moderate). Sorting by tier keeps that
# genuine divergence visible (Sergipe's point still sits right of
# Pernambuco's) while still grouping the panel into clean Low/Moderate/High
# blocks that read consistently with panel A.
setting_rank <- c(Low = 1, Moderate = 2, High = 3)
ar_summary <- ar_summary %>%
  mutate(setting = setting_key[region],
         setting = factor(setting, levels = c("Low", "Moderate", "High"))) %>%
  arrange(setting_rank[as.character(setting)], ar_med) %>%
  mutate(region = factor(region, levels = region))

# Same row order as ar_summary (panel B), so panels B and C line up state-
# for-state -- not peak_per_million's own order, which (by construction) IS
# monotonic within tier and would otherwise silently hide the divergence.
peak_rate_summary <- peak_rate_summary %>%
  mutate(setting = setting_key[region],
         setting = factor(setting, levels = c("Low", "Moderate", "High")),
         region = factor(region, levels = levels(ar_summary$region)))

# ---- C) Brazil state boundaries (all 27, for context) -----------------------
# geobr::read_state() is the more standard source for Brazilian IBGE state
# boundaries, but its metadata endpoint (download_metadata2()) was erroring
# out (object 'check_con' not found) when this was written -- rnaturalearth's
# admin-1 boundaries are simpler and sufficient for a state-level choropleth.
bra_states <- rnaturalearth::ne_states(country = "Brazil", returnclass = "sf")

bra_states <- bra_states %>%
  mutate(setting = setting_key[name],
         setting = factor(setting, levels = c("Low", "Moderate", "High")),
         focal = !is.na(setting))

focal_states <- bra_states %>% filter(focal)

# Zoom to the focal region: bbox of the 11 states + ~6% padding, rather than
# the full country (the whole point is to show these are 11 states WITHIN a
# much larger country, not to show all of Brazil).
bbox <- sf::st_bbox(focal_states)
pad_x <- (bbox["xmax"] - bbox["xmin"]) * 0.08
pad_y <- (bbox["ymax"] - bbox["ymin"]) * 0.08
xlim <- c(bbox["xmin"] - pad_x, bbox["xmax"] + pad_x)
ylim <- c(bbox["ymin"] - pad_y, bbox["ymax"] + pad_y)

# Label points: point-on-surface (not centroid) so labels for oddly-shaped/
# concave states still fall inside the polygon.
focal_labels <- focal_states %>%
  mutate(lbl_geom = sf::st_point_on_surface(geometry)) %>%
  sf::st_set_geometry("lbl_geom")

# ---- D) Panel A -- map -------------------------------------------------------
p_map <- ggplot() +
  geom_sf(data = bra_states, fill = "#F5F1E8", colour = "grey70", linewidth = 0.25) +
  geom_sf(data = focal_states, aes(fill = setting), colour = "grey30", linewidth = 0.3) +
  geom_sf_text(data = focal_labels, aes(label = name), size = 2.6, colour = "grey10",
               fontface = "bold", check_overlap = TRUE) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  scale_fill_manual(values = pal_setting, name = "Transmission setting", na.translate = FALSE) +
  labs(title = "A. States included (n = 11)") +
  theme_fig1() +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5))

# ---- E) Panel B -- posterior cumulative attack rate, sorted forest plot ----
p_forest <- ggplot(ar_summary, aes(x = ar_med, y = region, colour = setting)) +
  geom_errorbarh(aes(xmin = ar_lo, xmax = ar_hi), height = 0.2, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = pal_setting, name = "Transmission setting") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = expression(paste("Cumulative infection attack rate")),
       y = NULL,
       title = "B. Posterior infection attack rate by state") +
  guides(colour = "none") +  # dropped at the ggplot level -- theme(legend.position=...)
                              # alone isn't enough because plot_annotation()'s `&
                              # theme(legend.position = "bottom")` below re-applies to
                              # every panel and would silently re-enable it otherwise
  theme_fig1() +
  theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_text(colour = "black", size = 14, face = "bold"),
        plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0))

# ---- E2) Panel C -- peak weekly cases per million (the classification -----
# metric itself), same state order as panel B, with the two threshold lines
# (100 / 200 per million) that actually define the Low/Moderate/High cutoffs
# -- so it's visually obvious why each state landed in its tier, in contrast
# to panel B's different (posterior, whole-season) metric.
p_peak <- ggplot(peak_rate_summary, aes(x = peak_per_million, y = region, colour = setting)) +
  geom_vline(xintercept = c(100, 200), linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(size = 3) +
  scale_colour_manual(values = pal_setting, name = "Transmission setting") +
  labs(x = "Peak weekly cases per million",
       y = NULL,
       title = "C. Peak case rate by state") +
  guides(colour = "none") +
  theme_fig1() +
  theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
        axis.text.y = element_blank())

# ---- F) Combine --------------------------------------------------------------
p_fig1 <- (p_map | p_forest | p_peak) +
  plot_layout(widths = c(1.1, 1, 1), guides = "collect") +
  plot_annotation(
    #title = "Figure 1. Epidemiological reconstruction and analytical framework",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  ) &
  theme(legend.position = "bottom")

ggsave("06_Results/figure1_epi_map_ar.png", p_fig1, width = 17, height = 7, dpi = 300, bg = "white")
message("Saved: 06_Results/figure1_epi_map_ar.png")
