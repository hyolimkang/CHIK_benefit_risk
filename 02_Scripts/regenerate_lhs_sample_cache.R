# =============================================================================
# regenerate_lhs_sample_cache.R
#
# Regenerates 01_Data/lhs_sample.RData from the CURRENT 01_setup.R +
# 02_setup_age_props.R, so 03_brazil_all_draws_ori_v3.R's cache-load (line
# ~177) reflects the current LHS formulas instead of a stale snapshot.
#
# Background (chat record 2026-08-31): the cache on disk was last saved
# 2026-06-15, before 01_setup.R was edited 2026-08-27 (LHS dimensionality
# changed from k=62 -- including 11 now-unused per-state ar_* columns -- down
# to k=51). Since lhs::randomLHS()'s design depends on total k, that change
# silently shifted EVERY column's values, not just the removed ones --
# verified empirically: draw 142's p_sae_vacc_u65/dw_hosp/dw_subac/etc all
# differed between the stale cache and a fresh run of current 01_setup.R,
# while finite_weeksweep_brr.R (which regenerates lhs_sample inline, never
# caches) matched the fresh run exactly.
#
# Run this from the project root; safe to re-run any time 01_setup.R or
# 02_setup_age_props.R changes.
# =============================================================================

setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")

source("02_Scripts/01_setup.R")
source("02_Scripts/02_setup_age_props.R")

stopifnot(is.data.frame(lhs_sample))
stopifnot(all(c("subac_o40", "chr6m_o40", "chr12m_o40", "chr30m_o40") %in% colnames(lhs_sample)))
stale_state_cols <- c("ar_ce", "ar_bh", "ar_pa", "ar_pn", "ar_rg", "ar_pi", "ar_tc", "ar_ag", "ar_mg", "ar_se", "ar_go")
stopifnot(!any(stale_state_cols %in% colnames(lhs_sample)))  # confirm the stale per-state ar_* columns are gone
stopifnot(all(lhs_sample$p_death_vacc_u65 == 0))       # confirm the 18-64 death=0 fix is present

save(lhs_sample, file = "01_Data/lhs_sample.RData")
message("Saved fresh 01_Data/lhs_sample.RData -- ", nrow(lhs_sample), " rows x ", ncol(lhs_sample), " cols")
message("Columns: ", paste(colnames(lhs_sample), collapse = ", "))
