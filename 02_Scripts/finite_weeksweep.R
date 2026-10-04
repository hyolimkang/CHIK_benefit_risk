### -----------------------------------------------------------------------
### Vaccination-timing x coverage sweep, finite-history (long-term FOI +
### 8-year exposure cap) model.
###
### Point-estimate only (posterior/LHS medians) -- NOT full LHS-UI. A
### 52-week start-week x 10-level coverage x 11-region x 4-scenario grid at
### full 1000-draw UI would be ~5.7M individual SIR runs; medians make it a
### ~23k-run grid instead (see design discussion this was built from).
###
### Self-contained: reuses fit_prevacc_*_finite (via posterior_list_finite)
### and lhs_combined_finite from the finite-history pipeline
### (age_struc_fitting_region_finite_2022.R / lhs_samples_sir_finite.R /
### model_0803_finite_history.R), but does NOT load, modify, or depend on
### any of the *session-global* state those scripts set (posterior_idx_list,
### lhs_idx_list, regions, coverage_set, etc). Calls
### sirv_sim_coverageSwitch_ixchiq() directly instead of going through
### run_simulation_scenarios_ui_ixchiq(), because that wrapper (a) hardcodes
### delay = 2 and (b) sources total_coverage from the region_coverage global
### lookup rather than its own total_coverage argument -- neither of which
### can be swept without changing the wrapper. All object/file names below
### are new; nothing produced by the existing finite-history scripts is
### overwritten.
### -----------------------------------------------------------------------

library(purrr)
library(dplyr)
library(ggplot2)
library(writexl)

source("01_Script/1_1_Functions/sim_functions_final.R")   # sirv_sim_coverageSwitch_ixchiq()

load("00_Data/0_2_Processed/lhs_combined_finite.RData")   # posterior_list_finite, lhs_combined_finite
load("00_Data/0_2_Processed/bra_pop_2022_cleaned.RData")  # N_bahia, N_ceara, N_mg, N_pemam, N_pa,
                                                            # N_rg, N_pi, N_ag, N_tc, N_se, N_go

age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
                mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
                mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
                mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))

# Matches make_finite_prevacc_data()'s default (introduction_year = 2014L,
# fitting_year = 2022L for all 11 states) in age_struc_fitting_region_finite_2022.R.
years_since_introduction <- 2022L - 2014L

N_by_region <- list(
  "Bahia"               = as.numeric(N_bahia$Bahia),
  "Ceará"               = as.numeric(N_ceara$Ceará),
  "Minas Gerais"        = as.numeric(N_mg$`Minas Gerais`),
  "Pernambuco"          = as.numeric(N_pemam$Pernambuco),
  "Paraíba"             = as.numeric(N_pa$Paraíba),
  "Rio Grande do Norte" = as.numeric(N_rg$`Rio Grande do Norte`),
  "Piauí"               = as.numeric(N_pi$Piauí),
  "Alagoas"             = as.numeric(N_ag$Alagoas),
  "Tocantins"           = as.numeric(N_tc$Tocantins),
  "Sergipe"             = as.numeric(N_se$Sergipe),
  "Goiás"               = as.numeric(N_go$Goiás)
)

scenario_labels <- c("1-11y", "12-17y", "18-64y", "65y+")

# Same four target_age vectors as model_0803_finite_history.R's target_age_list.
target_age_list <- list(
  c(0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0),
  c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1)
)

#' Sweep vaccination start week x total coverage for one region/scenario,
#' using posterior/LHS *medians* (point estimate -- no UI).
#'
#' base_beta/I0 come from posterior_list_finite[[region]] (the finite-history
#' Stan fit). gamma/rho/sigma/foi/ve_ix/wd come from lhs_combined_finite[[region]]
#' -- matching exactly which object each parameter is drawn from inside
#' run_simulation_scenarios_ui_ixchiq() (age_struc_fitting_region_func_updated.R
#' L2874-2882), just median instead of per-draw. R0 uses the same finite-history
#' cap as that function's finite_history = TRUE branch (L2909-2910): 1 - exp(-foi
#' * pmin(age_groups, years_since_introduction)).
simulate_weeksweep_finite <- function(
    region_name, target_age, scenario_label,
    posterior, lhs_sample, N,
    age_groups                = age_groups,
    years_since_introduction  = 8,
    week_range                = 1:52,
    coverage_range             = seq(0.1, 1, by = 0.1),
    T                          = 52L,
    A                          = 20L
) {

  base_beta_med <- apply(posterior$beta_observed, 2, median)
  I0_med        <- apply(posterior$I0,            2, median)
  gamma_med     <- median(lhs_sample$gamma)
  rho_med       <- median(lhs_sample$rho)
  sigma_med     <- median(lhs_sample$sigma)
  foi_med       <- median(lhs_sample$foi)
  ve_med        <- median(lhs_sample$ve_ix)
  wd_med        <- median(lhs_sample$wd)

  R0_med <- 1 - exp(-foi_med * pmin(age_groups, years_since_introduction))

  run_one <- function(delay_wk, coverage, target) {
    sim <- sirv_sim_coverageSwitch_ixchiq(
      T = T, A = A, N = N, r = rep(0, A),
      base_beta = base_beta_med, I0_draw = I0_med, R0 = R0_med,
      rho = rho_med, gamma = gamma_med, sigma = sigma_med,
      delay = delay_wk, target_age = target, total_coverage = coverage,
      weekly_delivery_speed = wd_med,
      VE_block = ve_med, VE_inf = ve_med, coverage_threshold = 1
    )
    sum(sim$true_symptomatic)
  }

  # No-vaccination baseline: delay past the simulation horizon so the vaccine
  # branch never fires, coverage/target moot but set to 0 for clarity.
  baseline_cases <- run_one(delay_wk = T + 1L, coverage = 0, target = rep(0, A))

  grid <- expand.grid(week = week_range, coverage = coverage_range)
  grid$total_cases_vacc <- map2_dbl(
    grid$week, grid$coverage,
    ~ run_one(delay_wk = .x, coverage = .y, target = target_age)
  )

  grid$region                <- region_name
  grid$scenario               <- scenario_label
  grid$total_cases_baseline   <- baseline_cases
  grid$cases_averted          <- baseline_cases - grid$total_cases_vacc
  grid$pct_averted            <- 100 * grid$cases_averted / baseline_cases

  grid[, c("region", "scenario", "week", "coverage",
           "total_cases_baseline", "total_cases_vacc",
           "cases_averted", "pct_averted")]
}

## Full grid: 11 regions x 4 scenarios x 52 weeks x 10 coverage levels --------

weeksweep_grid_finite <- map_dfr(names(N_by_region), function(region_name) {
  map_dfr(seq_along(target_age_list), function(s) {
    simulate_weeksweep_finite(
      region_name              = region_name,
      target_age                = target_age_list[[s]],
      scenario_label             = scenario_labels[s],
      posterior                  = posterior_list_finite[[region_name]],
      lhs_sample                 = lhs_combined_finite[[region_name]],
      N                          = N_by_region[[region_name]],
      age_groups                 = age_groups,
      years_since_introduction   = years_since_introduction
    )
  })
})

dir.create("02_Outputs/2_2_Tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("02_Outputs/2_1_Figures", showWarnings = FALSE, recursive = TRUE)

save(weeksweep_grid_finite, file = "00_Data/0_2_Processed/finite_weeksweep_grid.RData")
write_xlsx(weeksweep_grid_finite, path = "02_Outputs/2_2_Tables/finite_weeksweep_grid.xlsx")


## Heatmap: % of symptomatic cases averted, by start week x coverage --------
## Sequential metric (magnitude of impact) -> single hue, light-to-dark.

national_weeksweep <- weeksweep_grid_finite %>%
  group_by(scenario, week, coverage) %>%
  summarise(
    total_cases_baseline = sum(total_cases_baseline),
    cases_averted         = sum(cases_averted),
    .groups = "drop"
  ) %>%
  mutate(
    pct_averted = 100 * cases_averted / total_cases_baseline,
    scenario    = factor(scenario, levels = scenario_labels)
  )

p_weeksweep_heatmap <- ggplot(
  national_weeksweep,
  aes(x = week, y = factor(coverage * 100), fill = pct_averted)
) +
  geom_tile() +
  facet_wrap(~scenario, ncol = 2) +
  scale_fill_distiller(
    palette   = "Blues",
    direction = 1,
    name      = "% symptomatic\ncases averted"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(8), expand = c(0, 0)) +
  labs(
    title    = "Impact of vaccination timing and coverage, finite-history baseline",
    subtitle = "National total, by target-age scenario. Point estimate (posterior/LHS medians), no uncertainty interval.",
    x        = "Vaccination campaign start week",
    y        = "Total coverage (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title.position = "plot",
    strip.text          = element_text(face = "bold"),
    panel.grid           = element_blank()
  )

ggsave(
  "02_Outputs/2_1_Figures/finite_weeksweep_heatmap.png",
  plot = p_weeksweep_heatmap, width = 10, height = 7, dpi = 400, bg = "white"
)

print(p_weeksweep_heatmap)
