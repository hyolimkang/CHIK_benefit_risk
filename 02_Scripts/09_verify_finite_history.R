# =============================================================================
# 09_verify_finite_history.R
#
# Fast consistency check for the finite-history / short-term FOI fix before
# trusting the regenerated ft_brr / brr_table_final_long for the manuscript
# revision (see revision plan, Part A). Not a formal parameter-recovery
# exercise -- just sanity checks that nothing looks obviously wrong.
#
# Sections 1-2 run against the companion simulation repo
# (CHIK_vaccine_impact/CHIK_ORV_impact) and are self-contained (only load
# small/medium cached RData files, no re-fitting).
#
# Section 3 requires objects that only exist in an interactive session that
# has already run 03_brazil_all_draws_ori_v3.R through Section 08
# (brr_table_final_long / ft_brr are not saved to disk anywhere in the
# pipeline) -- run this section in that session, not standalone.
#
# Sections 4-5 run against this repo's saved 01_Data/all_draws_*_true.RData
# (already produced by the current DATA_MODE = "flat" run).
# =============================================================================

## ---- 1) FOI self-consistency: production draws vs bra_state_short_term_foi ----
verify_foi_self_consistency <- function(
    sim_repo = "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"
) {
  owd <- getwd(); on.exit(setwd(owd))
  setwd(sim_repo)

  load("00_Data/0_2_Processed/bra_state_short_term_foi.RData")  # -> results_df
  st_foi <- results_df
  load("00_Data/0_2_Processed/lhs_combined_finite.RData")       # -> lhs_combined_finite

  regions <- names(lhs_combined_finite)
  out <- do.call(rbind, lapply(regions, function(reg) {
    row <- st_foi[st_foi$state == reg, ]
    q <- quantile(lhs_combined_finite[[reg]]$foi, c(0.5, 0.025, 0.975), na.rm = TRUE)
    data.frame(region = reg, H_median = row$H_median, H_lo = row$H_lo, H_hi = row$H_hi,
               foi_med = q[1], foi_lo = q[2], foi_hi = q[3])
  }))
  out$median_delta_pct <- 100 * (out$foi_med - out$H_median) / out$H_median
  attr(out, "max_abs_delta_pct") <- max(abs(out$median_delta_pct))
  out
}

## ---- 2) Old (long-term avg FOI) vs new (short-term FOI) baseline immunity ----
verify_s0_shift <- function(
    sim_repo = "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact",
    age_cap_years = 2022 - 2014
) {
  owd <- getwd(); on.exit(setwd(owd))
  setwd(sim_repo)

  age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
                   mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
                   mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
                   mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))

  load("00_Data/0_2_Processed/lhs_combined_finite.RData")   # new (short-term) draws
  load("00_Data/0_2_Processed/bra_foi_states.RData")        # old (long-term, per-pixel) draws
  regions <- names(lhs_combined_finite)

  bra_foi_states_df <- if (inherits(bra_foi_states, "sf")) {
    sf::st_drop_geometry(bra_foi_states)
  } else bra_foi_states
  foi_cols <- paste0("foi", 1:100)

  set.seed(123)
  foi_draws_list <- lapply(regions, function(reg) {
    row <- bra_foi_states_df[bra_foi_states_df$NAME_1 == reg, ]
    vals <- unlist(row[, foi_cols], use.names = FALSE)
    sample(vals[!is.na(vals)], size = 1000, replace = TRUE)
  })
  names(foi_draws_list) <- regions

  out <- do.call(rbind, lapply(regions, function(reg) {
    foi_old_med <- median(foi_draws_list[[reg]])
    foi_new_med <- median(lhs_combined_finite[[reg]]$foi)
    S0_old <- mean(1 - exp(-foi_old_med * pmin(age_groups, age_cap_years)))
    S0_new <- mean(1 - exp(-foi_new_med * pmin(age_groups, age_cap_years)))
    data.frame(region = reg, foi_old_longterm = foi_old_med, foi_new_shortterm = foi_new_med,
               mean_S0_old = S0_old, mean_S0_new = S0_new,
               foi_ratio_new_over_old = foi_new_med / foi_old_med)
  }))
  out[order(-out$foi_ratio_new_over_old), ]
}

## ---- 3) Headline number: does the 65+ / mechanism finding survive? ----------
## RUN THIS IN THE LIVE SESSION where brr_table_final_long already exists
## (after 03_brazil_all_draws_ori_v3.R has completed Section 08B with
## DATA_MODE = "flat"). Re-run with DATA_MODE = "original" to get the
## comparison slice, or save brr_table_final_long to disk from each run:
##   save(brr_table_final_long, file = "01_Data/brr_table_final_long_flat.RData")
##   save(brr_table_final_long, file = "01_Data/brr_table_final_long_original.RData")
verify_65plus_mechanism_headline <- function(brr_table_final_long) {
  brr_table_final_long[brr_table_final_long$`Age group` == "65+", ]
}

## ---- 4) Data hygiene: NA / Inf / negative values in saved draw objects ------
verify_draw_hygiene <- function(
    data_dir = "c:/Users/user/OneDrive/CHIK_benefit_risk/01_Data",
    objs = c("all_draws_ix_true", "all_draws_hosp_true", "all_draws_fatal_true",
             "all_draws_daly_true", "all_draws_sae_true")
) {
  do.call(rbind, lapply(objs, function(f) {
    e <- new.env()
    load(file.path(data_dir, paste0(f, ".RData")), envir = e)
    df <- get(f, envir = e)
    data.frame(object = f, nrow = nrow(df),
               na_pre = sum(is.na(df$total_pre)), na_post = sum(is.na(df$total_post)),
               neg_pre = sum(df$total_pre < 0, na.rm = TRUE),
               neg_post = sum(df$total_post < 0, na.rm = TRUE),
               n_regions = length(unique(df$Region)))
  }))
}

## ---- 5) Burden-level plausibility across regions (flat mode) ----------------
verify_burden_ranking <- function(
    data_dir = "c:/Users/user/OneDrive/CHIK_benefit_risk/01_Data",
    object = "all_draws_daly_true"
) {
  e <- new.env()
  load(file.path(data_dir, paste0(object, ".RData")), envir = e)
  df <- get(object, envir = e)
  agg <- aggregate(total_pre ~ Region, data = df, FUN = median)
  agg[order(-agg$total_pre), ]
}

## ---- 6) Travel scenario scope check (informational only, no re-run needed) --
## compute_ar()/compute_outcome() in 01_setup.R use the fitted time-varying
## daily FOI from the actual 2022 outbreak (already state/outbreak-specific),
## not the catalytic-model S0 that changed under this fix -- travel-scenario
## BRR numbers should be UNCHANGED by this fix. Reviewer #1's flagged
## inconsistency (line 315: text says "risks outweigh benefits" for
## 90-day/18-64/low-transmission while BRR > 1, 91.5% probability of benefit)
## is a wording bug in the manuscript text, not a numbers bug -- fix in prose,
## no pipeline change needed.

if (sys.nframe() == 0) {
  cat("=== 1) FOI self-consistency ===\n")
  foi_check <- verify_foi_self_consistency()
  print(foi_check, digits = 4, row.names = FALSE)
  cat(sprintf("\nMax |median delta %%|: %.2f%% (small = consistent by construction)\n\n",
              attr(foi_check, "max_abs_delta_pct")))

  cat("=== 2) Old (long-term) vs new (short-term) baseline immunity ===\n")
  print(verify_s0_shift(), digits = 4, row.names = FALSE)

  cat("\n=== 4) Data hygiene (saved flat-mode draw objects) ===\n")
  print(verify_draw_hygiene(), row.names = FALSE)

  cat("\n=== 5) Burden ranking sanity (all_draws_daly_true, flat mode) ===\n")
  print(verify_burden_ranking(), row.names = FALSE)

  cat("\nSection 3 (65+ mechanism headline) must be run in the live session --",
      "see verify_65plus_mechanism_headline().\n")
}
