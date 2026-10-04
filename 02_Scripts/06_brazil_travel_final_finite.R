# =============================================================================
# 06_brazil_travel_final_finite.R
#
# Rebuilds psa_df for the TRAVEL vaccination scenario under the finite-history
# baseline-immunity assumption (DATA_MODE = "flat" throughout the rest of
# this project). Adapted from the final PSA block of 06_brazil_travel_final.R
# (the "using 1000 random draws of attack rates by each state" section,
# ~line 553-825 of that file -- the file has two earlier, superseded PSA
# loops above it that are NOT run here, since this block is the one whose
# output (psa_df) actually gets saved and consumed downstream by
# 07_brazil_travel_final_update.R).
#
# ARCHITECTURE (rewritten 2026-08, see chat record):
# Local transmission risk for a trip now comes DIRECTLY from the m3-refit
# Stan posterior's own `phi_pred` (weekly force-of-infection generated
# quantity), converted to a daily hazard and windowed by entry_day/duration
# via a rolling cumsum. Each PSA draw `d` is independently paired to a Stan
# posterior draw index (posterior_idx_by_state), so transmission uncertainty
# = Stan posterior, VE/natural-history/safety uncertainty = LHS, entry timing
# uncertainty = sampled entry day -- three separate uncertainty sources, not
# mixed.
#
# REMOVED: the previous "target AR" architecture --
#   ar_target (lhs_sample$ar_<state>) -> Ht = -log(1-ar_target) ->
#   m = Ht/H0 -> foi_daily_scaled = foi0*m
# -- which took only the WEEKLY SHAPE of a single resimulated VE0/cov50 phi
# trajectory and rescaled its total magnitude to hit a separately-derived
# per-state AR target (itself built via 02b_setup_ar_by_state.R from
# all_draws_ix_true_finite.RData / total_S0, i.e. a full forward-resimulated
# "pre-vaccination baseline" scenario -- a different quantity from the
# historical-fit AR by design, not a bug). 09_diagnose_ar_incidence_vs_hazard.R
# confirmed phi_pred is already internally self-consistent with the Stan fit's
# own new_exposed-based AR (ratio = 1 to machine precision, all 11 states),
# so no external rescaling is needed -- phi_pred is used as-is.
#
# symp_prop_d (p_sym) is still drawn per-draw from lhs_sample$symp_overall
# (uncertainty propagation kept, NOT fixed to the 0.524 point estimate used
# in the primary transmission refit -- deliberate choice, see chat record).
#
# Output:
#   - 01_Data/psa_df_bra_travel_finite.RData
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")
source("02_Scripts/01_setup.R")
source("02_Scripts/02_setup_age_props.R")
library(tictoc)

load("01_Data/rho_df.RData")

# posterior_list (keyed by full region name) -- previously obtained as a
# side effect of sourcing 02b_setup_ar_by_state.R (which built it only to
# then derive the now-removed ar_<state>/ar_target columns). Loaded directly
# here instead, since phi_pred is all this script needs from it.
load("01_Data/posterior_finite_all.RData")  # -> posterior_ce, posterior_bh, ...
posterior_list <- list(
  "Ceará"               = posterior_ce,
  "Bahia"                = posterior_bh,
  "Paraíba"              = posterior_pa,
  "Pernambuco"           = posterior_pn,
  "Rio Grande do Norte"  = posterior_rg,
  "Piauí"                = posterior_pi,
  "Alagoas"              = posterior_ag,
  "Tocantins"            = posterior_tc,
  "Minas Gerais"         = posterior_mg,
  "Sergipe"              = posterior_se,
  "Goiás"                = posterior_go
)
names(posterior_list) <- enc2utf8(names(posterior_list))

# ---- FOI shape+magnitude straight from the corrected Stan posterior --------
# NOTE (2026-08): previously built from the R-resimulator's own VE0/cov50
# rep-1 `sim_out$phi` (a single deterministic trajectory per state), then
# rescaled in magnitude to hit a separately-derived `lhs_sample$ar_<state>`
# target (target AR -> Ht -> multiplier m -> scaled FOI). Confirmed via
# 09_diagnose_ar_incidence_vs_hazard.R that AR_hazard = 1-exp(-sum(phi_pred))
# and AR_incidence (new_exposed/S_start) are IDENTICAL within the Stan fit
# itself (ratio = 1 to machine precision, all 11 states) -- so phi_pred is
# already the model's own self-consistent transmission trajectory and needs
# no external rescaling. Now: phi_pred is drawn per-PSA-draw directly from
# the m3-refit posterior (posterior_list, loaded directly from
# posterior_finite_all.RData above), independently paired to the LHS draw
# index via posterior_idx_by_state -- see chat record, 2026-08.
weekly_to_daily_foi <- function(phi_weekly) {
  phi_weekly <- as.numeric(phi_weekly)
  rep(phi_weekly / 7, each = 7)  # 52*7 = 364 days
}

states_to_run <- enc2utf8(names(posterior_list))

set.seed(123)
posterior_idx_by_state <- lapply(states_to_run, function(st) {
  sample(seq_len(nrow(posterior_list[[st]]$phi_pred)),
         size = nrow(lhs_sample), replace = TRUE)
})
names(posterior_idx_by_state) <- states_to_run

# ---- PSA loop: 1000 draws x 11 states x 2 age groups x 4 durations ---------
n_entry_samples <- 100

all_risk <- all_risk %>% dplyr::filter(!is.na(SAE))

psa_out_list <- list()
tic("PSA Total Run (finite-history)")

for (d in seq_len(nrow(lhs_sample))) {

  draw_pars <- lhs_sample[d, ]
  ve_d      <- lhs_sample$ve[d]

  p_sae_vacc_u65   <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65    <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65 <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65  <- lhs_sample$p_death_vacc_65[d]

  p_sae_nat_64 <- lhs_sample$p_sae_nat_64[d]
  p_sae_nat_65 <- lhs_sample$p_sae_nat_65[d]

  p_death_nat_64 <- lhs_sample$p_death_nat_64[d]
  p_death_nat_65 <- lhs_sample$p_death_nat_65[d]

  travel_days <- list(
    "7d"  = lhs_sample$trav_7d[d],
    "14d" = lhs_sample$trav_14d[d],
    "30d" = lhs_sample$trav_30d[d],
    "90d" = lhs_sample$trav_90d[d]
  )

  symp_prop_d <- as.numeric(draw_pars$symp_overall)

  for (st in states_to_run) {

    idx_phi <- posterior_idx_by_state[[st]][d]
    phi_weekly <- posterior_list[[st]]$phi_pred[idx_phi, ]
    foi_daily <- weekly_to_daily_foi(phi_weekly)

    H_total_state  <- sum(foi_daily, na.rm = TRUE)
    AR_total_state <- 1 - exp(-H_total_state)

    for (i in seq_len(nrow(all_risk))) {

      age <- as.character(all_risk$age_group[i])

      if (age == "65+") {
        p_sae_vacc   <- p_sae_vacc_65
        p_death_vacc <- p_death_vacc_65
      } else {
        p_sae_vacc   <- p_sae_vacc_u65
        p_death_vacc <- p_death_vacc_u65
      }

      # all_risk is pre-filtered to !is.na(SAE), which leaves only "18-64" and
      # "65+" -- manuscript scope excludes ages 1-17, same as the
      # outbreak-response pipeline.
      if (age == "18-64") {
        p_sae_nat   <- p_sae_nat_64
        p_death_nat <- p_death_nat_64
      } else if (age == "65+") {
        p_sae_nat   <- p_sae_nat_65
        p_death_nat <- p_death_nat_65
      } else {
        stop("Unexpected age group in all_risk (expected only 18-64/65+): ", age)
      }

      for (days_label in names(travel_days)) {

        D <- max(1L, round(travel_days[[days_label]]))
        L_days <- length(foi_daily)
        max_entry <- max(1L, L_days - D + 1L)
        entry_days <- sample.int(max_entry, size = n_entry_samples, replace = TRUE)

        # ---- Vectorised over all n_entry_samples entry days at once --------
        # Mathematically identical to calling compute_ar_travel() once per
        # entry day in a loop (verified numerically against the original
        # loop-based code with a fixed seed -- see chat record): a rolling
        # window sum via cumsum() replaces compute_ar_travel()'s per-call
        # slice-and-sum, and compute_outcome()/compute_daly_one_age_specific()
        # already operate elementwise (only +, -, *, ifelse()/case_when()),
        # so calling them once with length-n_entry_samples vector arguments
        # returns the same per-entry-day values as n_entry_samples separate
        # scalar calls -- just without paying R's per-call overhead 100x.
        if (D >= L_days) {
          H_trip <- rep(H_total_state, length(entry_days))
        } else {
          cs <- c(0, cumsum(foi_daily))
          H_trip <- cs[entry_days + D] - cs[entry_days]
        }
        AR_travel <- 1 - exp(-H_trip)

        res <- compute_outcome(
          AR           = AR_travel,
          p_hosp       = symp_prop_d * p_sae_nat,
          p_death      = symp_prop_d * p_death_nat,
          p_sae_vacc   = p_sae_vacc,
          p_death_vacc = p_death_vacc,
          VE_hosp      = ve_d,
          VE_death     = ve_d
        )

        # compute_outcome()'s internal ifelse(excess==0, NA, averted/excess)
        # takes its result LENGTH from the *test* argument. excess_10k_sae/
        # excess_10k_death are scalars (not AR-dependent) while averted_10k_*
        # are now length-100 vectors (AR_travel is vectorized), so
        # ifelse() silently truncates brr_sae/brr_death to length 1.
        # Recompute locally without ifelse() to keep full length + correct values.
        res$brr_sae   <- res$averted_10k_sae   / res$excess_10k_sae
        res$brr_sae[res$excess_10k_sae == 0]     <- NA_real_
        res$brr_death <- res$averted_10k_death / res$excess_10k_death
        res$brr_death[res$excess_10k_death == 0] <- NA_real_

        symp_nv_10k <- 1e4 * AR_travel * symp_prop_d
        symp_v_10k  <- symp_nv_10k * (1 - ve_d)

        nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - res$risk_nv_hosp)
        nonhosp_symp_v_10k  <- pmax(0, symp_v_10k  - res$risk_v_hosp)

        dz_nv <- compute_daly_one_age_specific(
          age_group = age,
          deaths_10k = res$risk_nv_death,
          hosp_10k = res$risk_nv_hosp,
          nonhosp_symp_10k = nonhosp_symp_nv_10k,
          symp_10k = symp_nv_10k,
          sae_10k = 0, deaths_sae_10k = 0,
          draw_pars = draw_pars
        )

        dz_v <- compute_daly_one_age_specific(
          age_group = age,
          deaths_10k = res$risk_v_death,
          hosp_10k = res$risk_v_hosp,
          nonhosp_symp_10k = nonhosp_symp_v_10k,
          symp_10k = symp_v_10k,
          sae_10k = 0, deaths_sae_10k = 0,
          draw_pars = draw_pars
        )

        # sae doesn't depend on entry day / AR at all -- compute once
        # (previously recomputed identically inside the 100x entry-day loop).
        sae <- compute_daly_one_age_specific(
          age_group = age,
          deaths_10k = 0, hosp_10k = 0, nonhosp_symp_10k = 0, symp_10k = 0,
          sae_10k = 1e4 * p_sae_vacc,
          deaths_sae_10k = 1e4 * p_death_vacc,
          draw_pars = draw_pars
        )

        daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
        daly_sae     <- sae$daly_sae  # scalar
        # same ifelse()-length bug as brr_sae/brr_death above: daly_sae is
        # scalar but daly_averted is now length-100, so avoid ifelse() here.
        brr_daly <- daly_averted / daly_sae
        brr_daly[daly_sae <= 0] <- NA_real_

        risk_nv_10k_sae <- res$risk_nv_hosp + res$risk_nv_death
        risk_v_10k_sae  <- (res$risk_v_hosp + res$risk_v_death) + res$excess_10k_sae

        psa_out_list[[length(psa_out_list) + 1]] <- data.frame(
          draw         = d,
          state        = st,
          H_total      = H_total_state,
          AR_total     = AR_total_state,
          AR_total_pct = AR_total_state * 100,
          days         = days_label,
          age_group    = age,
          entry_day    = median(entry_days),

          H_trip       = median(H_trip, na.rm = TRUE),
          AR           = median(AR_travel, na.rm = TRUE),

          risk_nv_10k_sae = median(risk_nv_10k_sae, na.rm = TRUE),
          risk_v_10k_sae  = median(risk_v_10k_sae,  na.rm = TRUE),
          averted_10k_sae = median(res$averted_10k_sae, na.rm = TRUE),
          excess_10k_sae  = median(res$excess_10k_sae,  na.rm = TRUE),
          brr_sae         = median(res$brr_sae,         na.rm = TRUE),

          risk_nv_10k_death = median(res$risk_nv_death, na.rm = TRUE),
          risk_v_10k_death  = median(res$risk_v_death + res$excess_10k_death, na.rm = TRUE),
          averted_10k_death = median(res$averted_10k_death, na.rm = TRUE),
          excess_10k_death  = median(res$excess_10k_death,  na.rm = TRUE),
          brr_death         = median(res$brr_death,         na.rm = TRUE),

          daly_nv      = median(dz_nv$daly_dz,   na.rm = TRUE),
          daly_v       = median(dz_v$daly_dz,    na.rm = TRUE),
          daly_averted = median(daly_averted,    na.rm = TRUE),
          daly_sae     = median(daly_sae,        na.rm = TRUE),
          brr_daly     = median(brr_daly,        na.rm = TRUE),

          yll_nv = median(dz_nv$yll_dz, na.rm = TRUE),
          yll_v  = median(dz_v$yll_dz,  na.rm = TRUE),
          yld_nv = median(dz_nv$yld_dz, na.rm = TRUE),
          yld_v  = median(dz_v$yld_dz,  na.rm = TRUE),

          row.names = NULL
        )
      }
    }
  }

  if (d %% 50 == 0) {
    cat("Completed", d, "of", nrow(lhs_sample), "PSA draws\n")
  }
}

toc()

psa_df <- dplyr::bind_rows(psa_out_list)
save(psa_df, file = "01_Data/psa_df_bra_travel_finite.RData")
cat("Saved: 01_Data/psa_df_bra_travel_finite.RData (", nrow(psa_df), "rows )\n")

# =============================================================================
# Fine entry_day x duration 2D grid, for a full heatmap (state/age/mechanism)
# -----------------------------------------------------------------------------
# Deliberately a SEPARATE loop from the discrete-duration PSA above, not an
# extension of it: the 7/14/30/90-day durations above are each drawn from
# lhs_sample$trav_7d/14d/30d/90d, which carry +-10% PSA jitter from
# 01_setup.R's core randomLHS(k = 51) call (line ~143-146). Adding grid
# duration points (e.g. 21d, 45d, ...) the same way would mean growing k,
# which regenerates the ENTIRE correlated LHS matrix and silently reshuffles
# every already-validated parameter draw in the pipeline -- not acceptable.
# So this grid uses deterministic (unjittered) duration values instead; it is
# a genuinely different quantity from the discrete table above, not a superset
# of it, even though both draw on the same underlying FOI/risk machinery.
#
# entry_day is also NOT Monte-Carlo-sampled-then-median-collapsed here (as it
# is above) -- it's a fixed weekly grid, kept as its own output axis, because
# the whole point of this block is to show how risk varies BY entry timing.
#
# Output is a compact per-cell summary across the 1000 PSA draws (median +
# Pr(BRR>1)), not a draw-level long table: storing state x age x duration x
# entry_day x draw individually would be ~31M rows. Values are instead
# accumulated into plain numeric arrays (draw as the last dimension) and
# collapsed with apply(..., median/mean) once the draw loop finishes -- this
# avoids ever materialising a giant long-format data.frame.
# =============================================================================

duration_grid  <- sort(unique(c(seq(7, 180, by = 7), 30, 90, 180)))  # 28 pts; 30/90 forced in (not on the 7d cadence) so they line up exactly with the discrete table's 30d/90d for a direct sanity check
entry_day_grid <- seq(1, 364, by = 7)                          # 52 points, weekly
age_levels     <- c("18-64", "65+")

n_state <- length(states_to_run)
n_age   <- length(age_levels)
n_dur   <- length(duration_grid)
n_ent   <- length(entry_day_grid)
n_draws <- nrow(lhs_sample)

arr_dimnames <- list(states_to_run, age_levels,
                      as.character(duration_grid), as.character(entry_day_grid), NULL)
new_grid_arr <- function() {
  array(NA_real_, dim = c(n_state, n_age, n_dur, n_ent, n_draws), dimnames = arr_dimnames)
}
grid_AR        <- new_grid_arr()
grid_brr_sae   <- new_grid_arr()
grid_brr_death <- new_grid_arr()
grid_brr_daly  <- new_grid_arr()
# Benefit/risk components (per 10,000), not just their ratio -- needed for
# 14_headline_summary_tables.R's Benefit/Risk/BRR summary, which pools
# across the whole duration x entry_day surface instead of picking one
# discrete duration.
grid_averted_sae   <- new_grid_arr()
grid_excess_sae    <- new_grid_arr()
grid_averted_death <- new_grid_arr()
grid_excess_death  <- new_grid_arr()
grid_daly_averted  <- new_grid_arr()
grid_daly_sae      <- new_grid_arr()

tic("Fine entry_day x duration grid (finite-history)")

for (d in seq_len(n_draws)) {

  draw_pars <- lhs_sample[d, ]
  ve_d      <- lhs_sample$ve[d]

  p_sae_vacc_u65   <- lhs_sample$p_sae_vacc_u65[d]
  p_sae_vacc_65    <- lhs_sample$p_sae_vacc_65[d]
  p_death_vacc_u65 <- lhs_sample$p_death_vacc_u65[d]
  p_death_vacc_65  <- lhs_sample$p_death_vacc_65[d]

  p_sae_nat_64 <- lhs_sample$p_sae_nat_64[d]
  p_sae_nat_65 <- lhs_sample$p_sae_nat_65[d]

  p_death_nat_64 <- lhs_sample$p_death_nat_64[d]
  p_death_nat_65 <- lhs_sample$p_death_nat_65[d]

  symp_prop_d <- as.numeric(draw_pars$symp_overall)

  for (s_idx in seq_along(states_to_run)) {

    st  <- states_to_run[s_idx]

    idx_phi <- posterior_idx_by_state[[st]][d]
    phi_weekly <- posterior_list[[st]]$phi_pred[idx_phi, ]
    foi_daily <- weekly_to_daily_foi(phi_weekly)

    L_days <- length(foi_daily)
    cs <- c(0, cumsum(foi_daily))
    H_full <- sum(foi_daily, na.rm = TRUE)

    for (dur_idx in seq_along(duration_grid)) {

      D <- duration_grid[dur_idx]
      max_entry  <- max(1L, L_days - D + 1L)
      valid_mask <- entry_day_grid <= max_entry
      ent        <- entry_day_grid[valid_mask]
      if (length(ent) == 0L) next

      if (D >= L_days) {
        AR_travel <- rep(1 - exp(-H_full), length(ent))
      } else {
        AR_travel <- 1 - exp(-(cs[ent + D] - cs[ent]))
      }

      for (age_idx in seq_along(age_levels)) {

        age <- age_levels[age_idx]

        if (age == "65+") {
          p_sae_vacc   <- p_sae_vacc_65
          p_death_vacc <- p_death_vacc_65
          p_sae_nat    <- p_sae_nat_65
          p_death_nat  <- p_death_nat_65
        } else {
          p_sae_vacc   <- p_sae_vacc_u65
          p_death_vacc <- p_death_vacc_u65
          p_sae_nat    <- p_sae_nat_64
          p_death_nat  <- p_death_nat_64
        }

        res <- compute_outcome(
          AR           = AR_travel,
          p_hosp       = symp_prop_d * p_sae_nat,
          p_death      = symp_prop_d * p_death_nat,
          p_sae_vacc   = p_sae_vacc,
          p_death_vacc = p_death_vacc,
          VE_hosp      = ve_d,
          VE_death     = ve_d
        )
        # same ifelse()-length fix as the discrete-duration loop above.
        res$brr_sae   <- res$averted_10k_sae   / res$excess_10k_sae
        res$brr_sae[res$excess_10k_sae == 0]     <- NA_real_
        res$brr_death <- res$averted_10k_death / res$excess_10k_death
        res$brr_death[res$excess_10k_death == 0] <- NA_real_

        symp_nv_10k <- 1e4 * AR_travel * symp_prop_d
        symp_v_10k  <- symp_nv_10k * (1 - ve_d)
        nonhosp_symp_nv_10k <- pmax(0, symp_nv_10k - res$risk_nv_hosp)
        nonhosp_symp_v_10k  <- pmax(0, symp_v_10k  - res$risk_v_hosp)

        dz_nv <- compute_daly_one_age_specific(
          age_group = age, deaths_10k = res$risk_nv_death, hosp_10k = res$risk_nv_hosp,
          nonhosp_symp_10k = nonhosp_symp_nv_10k, symp_10k = symp_nv_10k,
          sae_10k = 0, deaths_sae_10k = 0, draw_pars = draw_pars
        )
        dz_v <- compute_daly_one_age_specific(
          age_group = age, deaths_10k = res$risk_v_death, hosp_10k = res$risk_v_hosp,
          nonhosp_symp_10k = nonhosp_symp_v_10k, symp_10k = symp_v_10k,
          sae_10k = 0, deaths_sae_10k = 0, draw_pars = draw_pars
        )
        sae <- compute_daly_one_age_specific(
          age_group = age, deaths_10k = 0, hosp_10k = 0, nonhosp_symp_10k = 0, symp_10k = 0,
          sae_10k = 1e4 * p_sae_vacc, deaths_sae_10k = 1e4 * p_death_vacc, draw_pars = draw_pars
        )

        daly_averted <- dz_nv$daly_dz - dz_v$daly_dz
        daly_sae     <- sae$daly_sae
        brr_daly <- daly_averted / daly_sae
        brr_daly[daly_sae <= 0] <- NA_real_

        grid_AR[s_idx, age_idx, dur_idx, valid_mask, d]        <- AR_travel
        grid_brr_sae[s_idx, age_idx, dur_idx, valid_mask, d]   <- res$brr_sae
        grid_brr_death[s_idx, age_idx, dur_idx, valid_mask, d] <- res$brr_death
        grid_brr_daly[s_idx, age_idx, dur_idx, valid_mask, d]  <- brr_daly

        grid_averted_sae[s_idx, age_idx, dur_idx, valid_mask, d]   <- res$averted_10k_sae
        grid_excess_sae[s_idx, age_idx, dur_idx, valid_mask, d]    <- res$excess_10k_sae
        grid_averted_death[s_idx, age_idx, dur_idx, valid_mask, d] <- res$averted_10k_death
        grid_excess_death[s_idx, age_idx, dur_idx, valid_mask, d]  <- res$excess_10k_death
        grid_daly_averted[s_idx, age_idx, dur_idx, valid_mask, d]  <- daly_averted
        grid_daly_sae[s_idx, age_idx, dur_idx, valid_mask, d]      <- daly_sae
      }
    }
  }

  if (d %% 50 == 0) cat("Grid: completed", d, "of", n_draws, "PSA draws\n")
}

toc()

# ---- Collapse across draws into a compact per-cell summary -----------------
# 95% UI (2.5%/97.5% quantiles) added alongside median (chat record
# 2026-08-27): the raw per-draw array is already fully in memory at this
# point, so quantile() costs essentially nothing extra on top of median() --
# both are sort-based; this does NOT require re-running the simulation loop,
# only changes what's extracted from data already computed.
summarise_grid <- function(arr, prefix) {
  med    <- apply(arr, c(1, 2, 3, 4), median, na.rm = TRUE)
  lo     <- apply(arr, c(1, 2, 3, 4), quantile, probs = 0.025, na.rm = TRUE)
  hi     <- apply(arr, c(1, 2, 3, 4), quantile, probs = 0.975, na.rm = TRUE)
  pr_gt1 <- apply(arr, c(1, 2, 3, 4), function(x) mean(x > 1, na.rm = TRUE))
  med_df <- as.data.frame.table(med, responseName = paste0(prefix, "_median"))
  lo_df  <- as.data.frame.table(lo,  responseName = paste0(prefix, "_lo"))
  hi_df  <- as.data.frame.table(hi,  responseName = paste0(prefix, "_hi"))
  pr_df  <- as.data.frame.table(pr_gt1, responseName = paste0(prefix, "_prob_gt1"))
  med_df %>%
    dplyr::left_join(lo_df, by = c("Var1", "Var2", "Var3", "Var4")) %>%
    dplyr::left_join(hi_df, by = c("Var1", "Var2", "Var3", "Var4")) %>%
    dplyr::left_join(pr_df, by = c("Var1", "Var2", "Var3", "Var4"))
}

grid_AR_df <- {
  med <- apply(grid_AR, c(1, 2, 3, 4), median, na.rm = TRUE)
  lo  <- apply(grid_AR, c(1, 2, 3, 4), quantile, probs = 0.025, na.rm = TRUE)
  hi  <- apply(grid_AR, c(1, 2, 3, 4), quantile, probs = 0.975, na.rm = TRUE)
  as.data.frame.table(med, responseName = "AR_median") %>%
    dplyr::left_join(as.data.frame.table(lo, responseName = "AR_lo"), by = c("Var1","Var2","Var3","Var4")) %>%
    dplyr::left_join(as.data.frame.table(hi, responseName = "AR_hi"), by = c("Var1","Var2","Var3","Var4"))
}
grid_sae_df   <- summarise_grid(grid_brr_sae,   "brr_sae")
grid_death_df <- summarise_grid(grid_brr_death, "brr_death")
grid_daly_df  <- summarise_grid(grid_brr_daly,  "brr_daly")

# Benefit/risk components -- absolute per-10,000 quantities, not ratios, so
# no Pr(>1) is meaningful, but median + 95% UI still are.
median_only_df <- function(arr, name) {
  base_name <- sub("_median$", "", name)
  med <- apply(arr, c(1, 2, 3, 4), median, na.rm = TRUE)
  lo  <- apply(arr, c(1, 2, 3, 4), quantile, probs = 0.025, na.rm = TRUE)
  hi  <- apply(arr, c(1, 2, 3, 4), quantile, probs = 0.975, na.rm = TRUE)
  out <- as.data.frame.table(med, responseName = name) %>%
    dplyr::left_join(as.data.frame.table(lo, responseName = "lo_tmp"), by = c("Var1","Var2","Var3","Var4")) %>%
    dplyr::left_join(as.data.frame.table(hi, responseName = "hi_tmp"), by = c("Var1","Var2","Var3","Var4"))
  names(out)[names(out) == "lo_tmp"] <- paste0(base_name, "_lo")
  names(out)[names(out) == "hi_tmp"] <- paste0(base_name, "_hi")
  out
}
grid_averted_sae_df   <- median_only_df(grid_averted_sae,   "averted_sae_median")
grid_excess_sae_df    <- median_only_df(grid_excess_sae,    "excess_sae_median")
grid_averted_death_df <- median_only_df(grid_averted_death, "averted_death_median")
grid_excess_death_df  <- median_only_df(grid_excess_death,  "excess_death_median")
grid_daly_averted_df  <- median_only_df(grid_daly_averted,  "daly_averted_median")
grid_daly_sae_df      <- median_only_df(grid_daly_sae,      "daly_sae_median")

psa_grid_travel <- grid_AR_df %>%
  dplyr::left_join(grid_sae_df,   by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_death_df, by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_daly_df,  by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_averted_sae_df,   by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_excess_sae_df,    by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_averted_death_df, by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_excess_death_df,  by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_daly_averted_df,  by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::left_join(grid_daly_sae_df,      by = c("Var1", "Var2", "Var3", "Var4")) %>%
  dplyr::rename(state = Var1, age_group = Var2, duration = Var3, entry_day = Var4) %>%
  dplyr::mutate(
    duration  = as.integer(as.character(duration)),
    entry_day = as.integer(as.character(entry_day))
  )

save(psa_grid_travel, file = "01_Data/psa_grid_bra_travel_finite.RData")
cat("Saved: 01_Data/psa_grid_bra_travel_finite.RData (", nrow(psa_grid_travel), "rows )\n")

# ---- Setting-level (Low/Moderate/High) pooled summary -----------------------
# Same convention as 11_travel_headline_figures.R's build_travel_synthesis_data:
# pool states WITHIN a setting together with draws (not median-of-state-medians)
# before summarising, using setting_key.RData (case-rate-based classification,
# 08_attack_rate_classification_finite.R). Only feasible here because the raw
# per-state x per-draw grid arrays are still in memory at this point in the
# script -- they're discarded after this block, which is why this summary
# can't be reconstructed later from psa_grid_bra_travel_finite.RData alone.
load("01_Data/setting_key.RData")  # -> setting_key (region -> Low/Moderate/High)

# 95% UI (2.5%/97.5% quantiles) added alongside median (chat record
# 2026-08-27) -- `merged` (all draws x states-in-setting) is already fully
# built in memory right here before collapsing, so quantile() is essentially
# free on top of median(): no new simulation, just extracting 2 more numbers
# per cell from data that already exists at this point.
summarise_grid_by_setting <- function(arr, prefix, has_prob = TRUE) {
  out_list <- list()
  for (s in c("Low", "Moderate", "High")) {
    st_in_setting <- states_to_run[setting_key[states_to_run] == s]
    idx <- match(st_in_setting, states_to_run)
    sub  <- arr[idx, , , , , drop = FALSE]           # state_sub x age x dur x ent x draws
    perm <- aperm(sub, c(2, 3, 4, 1, 5))              # age x dur x ent x state_sub x draws
    d <- dim(perm)
    merged <- array(perm, dim = c(d[1], d[2], d[3], d[4] * d[5]),
                     dimnames = list(dimnames(perm)[[1]], dimnames(perm)[[2]], dimnames(perm)[[3]], NULL))
    med <- apply(merged, c(1, 2, 3), median, na.rm = TRUE)
    lo  <- apply(merged, c(1, 2, 3), quantile, probs = 0.025, na.rm = TRUE)
    hi  <- apply(merged, c(1, 2, 3), quantile, probs = 0.975, na.rm = TRUE)
    med_df <- as.data.frame.table(med, responseName = paste0(prefix, "_median"))
    lo_df  <- as.data.frame.table(lo,  responseName = paste0(prefix, "_lo"))
    hi_df  <- as.data.frame.table(hi,  responseName = paste0(prefix, "_hi"))
    med_df <- med_df %>%
      dplyr::left_join(lo_df, by = c("Var1", "Var2", "Var3")) %>%
      dplyr::left_join(hi_df, by = c("Var1", "Var2", "Var3"))
    if (has_prob) {
      pr_gt1 <- apply(merged, c(1, 2, 3), function(x) mean(x > 1, na.rm = TRUE))
      pr_df  <- as.data.frame.table(pr_gt1, responseName = paste0(prefix, "_prob_gt1"))
      med_df <- dplyr::left_join(med_df, pr_df, by = c("Var1", "Var2", "Var3"))
    }
    med_df$setting <- s
    out_list[[s]] <- med_df
  }
  dplyr::bind_rows(out_list)
}

grid_AR_setting_df    <- summarise_grid_by_setting(grid_AR,        "AR",        has_prob = FALSE)
grid_sae_setting_df   <- summarise_grid_by_setting(grid_brr_sae,   "brr_sae")
grid_death_setting_df <- summarise_grid_by_setting(grid_brr_death, "brr_death")
grid_daly_setting_df  <- summarise_grid_by_setting(grid_brr_daly,  "brr_daly")

grid_averted_sae_setting_df   <- summarise_grid_by_setting(grid_averted_sae,   "averted_sae",   has_prob = FALSE)
grid_excess_sae_setting_df    <- summarise_grid_by_setting(grid_excess_sae,    "excess_sae",    has_prob = FALSE)
grid_averted_death_setting_df <- summarise_grid_by_setting(grid_averted_death, "averted_death", has_prob = FALSE)
grid_excess_death_setting_df  <- summarise_grid_by_setting(grid_excess_death,  "excess_death",  has_prob = FALSE)
grid_daly_averted_setting_df  <- summarise_grid_by_setting(grid_daly_averted,  "daly_averted",  has_prob = FALSE)
grid_daly_sae_setting_df      <- summarise_grid_by_setting(grid_daly_sae,      "daly_sae",      has_prob = FALSE)

psa_grid_travel_setting <- grid_AR_setting_df %>%
  dplyr::left_join(grid_sae_setting_df,   by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_death_setting_df, by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_daly_setting_df,  by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_averted_sae_setting_df,   by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_excess_sae_setting_df,    by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_averted_death_setting_df, by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_excess_death_setting_df,  by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_daly_averted_setting_df,  by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::left_join(grid_daly_sae_setting_df,      by = c("Var1", "Var2", "Var3", "setting")) %>%
  dplyr::rename(age_group = Var1, duration = Var2, entry_day = Var3) %>%
  dplyr::mutate(
    duration  = as.integer(as.character(duration)),
    entry_day = as.integer(as.character(entry_day)),
    setting   = factor(setting, levels = c("Low", "Moderate", "High"))
  )

save(psa_grid_travel_setting, file = "01_Data/psa_grid_bra_travel_finite_setting.RData")
cat("Saved: 01_Data/psa_grid_bra_travel_finite_setting.RData (", nrow(psa_grid_travel_setting), "rows )\n")
