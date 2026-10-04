### -----------------------------------------------------------------------
### BRR (benefit-risk ratio) step for the week x coverage sweep -- the layer
### between combined_nnv_df_region_coverage_model_finite_weeksweep.RData
### (finite_weeksweep_nnv.R output) and a final (week x coverage x region x
### scenario) BRR table for the heatmap.
###
### Adapts 03_brazil_all_draws_ori_v3.R's Section 02-08B logic (in
### C:/Users/user/OneDrive/CHIK_benefit_risk/02_Scripts/03_brazil_all_draws_ori_v3.R),
### generalized in two ways that script hardcodes:
###   1) make_averted_draws_true() there defaults to coverage_keep = "cov50"
###      and tot_vacc_map_true explicitly filters VC == "cov50" -- both
###      dropped here so all 4 coverage tags survive.
###   2) brr_draw_summary_true's group_by() there has no Coverage column (it
###      never needed one, since only cov50 ever reached it) -- Coverage and
###      week are added here.
###
### NOT replicated (out of scope for a BRR heatmap): Section 08C (Word/
### flextable export), 08D (CEAC probability-curve plots), Section 09
### (age-bin diagnostic). Those are reporting/diagnostic layers on top of
### draw_level_xy_serostatus_all, which this script does produce -- rebuild
### them from that object later if needed.
###
### External function/data dependencies copied in verbatim (self-contained,
### no side effects, unlike sourcing 01_setup.R which starts with
### rm(list = ls(all = TRUE))):
###   - lhs_sample construction + compute_daly_one()   <- 01_setup.R
###   - age-specific extension + compute_daly_one_age_specific()
###                                                     <- 02_setup_age_props.R
### If those upstream files change, re-sync the copies below.
###
### Run this AFTER finite_weeksweep_fulldraws.R (postsim_vc_ixchiq_model_finite_week*.RData)
### AND finite_weeksweep_nnv.R (combined_nnv_df_region_coverage_model_finite_weeksweep.RData)
### have both completed.
### -----------------------------------------------------------------------

library(purrr)
library(dplyr)
library(tidyr)
library(lhs)
library(truncnorm)

WEEKS <- c(1, 8, 16, 24, 32, 42, 52)  # must match the other two sweep scripts
FILE_SUFFIX <- "_ve0"  # must match the VE_TAG run of finite_weeksweep_fulldraws.R/_nnv.R being consumed ("" for VE98.9, "_ve0" for VE0)

load("00_Data/0_2_Processed/lhs_combined_finite.RData")   # posterior_list_finite
load("00_Data/0_2_Processed/rho_df_finite.RData")          # rho_df (m3 refit)
load(sprintf("00_Data/0_2_Processed/combined_nnv_df_region_coverage_model_finite_weeksweep%s.RData", FILE_SUFFIX))  # all_weeks_nnv (region/VE/VC/week/scenario/AgeGroup grain)

region_names <- c("Bahia", "Ceará", "Minas Gerais", "Pernambuco", "Paraíba",
                  "Rio Grande do Norte", "Piauí", "Alagoas", "Tocantins",
                  "Sergipe", "Goiás")

lhs_sample_young <- readRDS("00_Data/0_2_Processed/lhs_sample_young.RDS")
lhs_old          <- readRDS("00_Data/0_2_Processed/lhs_old.RDS")

# preui_flat_all.RData = pre-vaccination baseline, per-draw (age x week)
# symptomatic-case matrices, one list per region. Built by
# 01_Script/1_4_Analysis/preui_flat_run.R via simulate_pre_ui_age() with the
# same finite-history fits/LHS this whole sweep uses. Week/coverage
# independent (pre-vaccination baseline doesn't depend on the vaccination
# campaign), so it's reused as-is across every week x coverage cell -- no
# need to regenerate it per week.
load("00_Data/0_2_Processed/preui_flat_all.RData")
preui_all <- list(
  "Ceará" = preui_ce, "Bahia" = preui_bh, "Paraíba" = preui_pa,
  "Pernambuco" = preui_pn, "Rio Grande do Norte" = preui_rg,
  "Piauí" = preui_pi, "Alagoas" = preui_ag, "Tocantins" = preui_tc,
  "Minas Gerais" = preui_mg, "Sergipe" = preui_se, "Goiás" = preui_go
)

setting_key <- c(
  "Ceará"               = "High",
  "Bahia"               = "Low",
  "Paraíba"             = "High",
  "Pernambuco"          = "Moderate",
  "Rio Grande do Norte" = "Low",
  "Piauí"               = "High",
  "Tocantins"           = "Moderate",
  "Alagoas"             = "High",
  "Minas Gerais"        = "Low",
  "Sergipe"             = "Low",
  "Goiás"               = "Low"
)

age_map <- data.frame(
  age_index = 1:20,
  age_gr = c(
    "<1", "1-4", "5-9", "10-11", "12-17", "18-19", "20-24", "25-29",
    "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
    "60-64", "65-69", "70-74", "75-79", "80-84", "85+"
  ),
  AgeCat = c(
    "1-11", "1-11", "1-11", "1-11",
    "12-17",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "18-64", "18-64", "18-64", "18-64", "18-64",
    "65+", "65+", "65+", "65+", "65+"
  )
)

rr_vals <- c(1.0, 0.5, 0.1, 0.0)

map_scenario_agecat_int <- function(scen_int) {
  dplyr::case_when(
    scen_int == 1L ~ "1-11",
    scen_int == 2L ~ "12-17",
    scen_int == 3L ~ "18-64",
    scen_int == 4L ~ "65+",
    TRUE ~ NA_character_
  )
}


## ---- lhs_sample (RISK/DALY parameters) + compute_daly_one() ----------
## Verbatim copy of C:/Users/user/OneDrive/CHIK_benefit_risk/02_Scripts/01_setup.R
## L112-244 and L248-362. NOTE this "lhs_sample" is unrelated to
## lhs_combined_finite / lhs_sample_young -- it's the benefit-risk repo's
## separate vaccine-risk/DALY-parameter LHS draw table (p_sae_vacc_*,
## le_lost_*, dw_*, dur_*, acute/subac/chr*m, etc.), 1000 draws.

set.seed(123)
runs <- 1000
A <- lhs::randomLHS(n = runs, k = 51)
lhs_sample <- matrix(nrow = nrow(A), ncol = ncol(A))

lhs_sample[,1]  <- qbeta(A[,1],  shape1 = 6+0.5,  shape2 = 32949-6+0.5)
# p_sae_vacc_65: 19 non-fatal SAEs / 18,445 65+ vaccinees, matching
# 01_Data/all_risk_four_age.csv exactly (was 20 -- an unexplained +1 not
# matching the CSV source, and double-counting the 1 fatal case since
# sae_10k_base = p_sae_vacc_base + p_death_vacc_base already adds the death
# back on top; see 01_setup.R's identical fix, chat record 2026-09-02).
lhs_sample[,2]  <- qbeta(A[,2],  shape1 = 19+0.5, shape2 = 18445-19+0.5)
# p_death_vacc_u65: fixed at exactly 0 -- see 01_setup.R's identical fix
# (chat record 2026-08-27): 18-64 vaccine-attributable death assumed
# impossible, not merely rare (was qbeta(0.5, 32949.5), small but nonzero).
lhs_sample[,3]  <- 0
lhs_sample[,4]  <- qbeta(A[,4],  shape1 = 1+0.5,  shape2 = 18445-1+0.5)

lhs_sample[,5]  <- qbeta(A[,5],  shape1 = 3703+1,  shape2 = 92763-3703+1)
lhs_sample[,6]  <- qbeta(A[,6],  shape1 = 1568+1,  shape2 = 50598-1568+1)
lhs_sample[,7]  <- qbeta(A[,7],  shape1 = 11344+1, shape2 = 396351-11344+1)
lhs_sample[,8]  <- qbeta(A[,8],  shape1 = 3690+1,  shape2 = 303588-3690+1)

lhs_sample[,9]  <- qbeta(A[,9],  shape1 = 32+1,  shape2 = 92763-32+1)
lhs_sample[,10] <- qbeta(A[,10], shape1 = 22+1,  shape2 = 50598-22+1)
lhs_sample[,11] <- qbeta(A[,11], shape1 = 312+1, shape2 = 396351-312+1)
lhs_sample[,12] <- qbeta(A[,12], shape1 = 534+1, shape2 = 303588-534+1)

lhs_sample[,13] <- qtruncnorm(A[,13], a = 0.967, b = 0.998, mean = 0.989, sd = (0.998 - 0.967) / (2 * qnorm(0.975)))
lhs_sample[,14] <- qunif(A[,14], min = 0.1 * 0.9, max = 0.1 * 1.1)
lhs_sample[,15] <- qunif(A[,15], min = 0.2 * 0.9, max = 0.2 * 1.1)
lhs_sample[,16] <- qunif(A[,16], min = 0.3 * 0.9, max = 0.3 * 1.1)
lhs_sample[,17] <- qunif(A[,17], min = 9 * 0.9,   max = 9 * 1.1)
lhs_sample[,18] <- qunif(A[,18], min = 7 * 0.9,   max = 7 * 1.1)
lhs_sample[,19] <- qunif(A[,19], min = 14 * 0.9,  max = 14 * 1.1)
lhs_sample[,20] <- qunif(A[,20], min = 30 * 0.9,  max = 30 * 1.1)
lhs_sample[,21] <- qunif(A[,21], min = 90 * 0.9,  max = 90 * 1.1)

lhs_sample[,22] <- qbeta(A[,22], shape1 = 49.14034, shape2 = 34.14837)
lhs_sample[,23] <- qbeta(A[,23], shape1 = 34.21298, shape2 = 31.963)
lhs_sample[,24] <- qbeta(A[,24], shape1 = 36.77819, shape2 = 32.7458)
lhs_sample[,25] <- qbeta(A[,25], shape1 = 35.84287, shape2 = 32.55955)
lhs_sample[,26] <- qbeta(A[,26], shape1 = 539.2823, shape2 = 14152.25)
lhs_sample[,27] <- qbeta(A[,27], shape1 = 58.96698, shape2 = 1415.207)
lhs_sample[,28] <- qbeta(A[,28], shape1 = 115.4225, shape2 = 111.383)

lhs_sample[,29] <- qlnorm(A[,29], meanlog = 4.26127,   sdlog = 0.0511915)
lhs_sample[,30] <- qlnorm(A[,30], meanlog = 4.119037,  sdlog = 0.0511915)
lhs_sample[,31] <- qlnorm(A[,31], meanlog = 3.583519,  sdlog = 0.0511915)
lhs_sample[,32] <- qlnorm(A[,32], meanlog = 1.808289,  sdlog = 0.0511915)

lhs_sample[,33] <- qbeta(A[,33],  shape1 = 4.581639, shape2 = 9.87148)
lhs_sample[,34] <- qlnorm(A[,34], meanlog = -0.6301724, sdlog = 0.0852)
lhs_sample[,35] <- qbeta(A[,35],  shape1 = 22.51835, shape2 = 146.7926)
lhs_sample[,36] <- qlnorm(A[,36], meanlog = -3.734278,  sdlog = 0.3877361)
lhs_sample[,37] <- qbeta(A[,37],  shape1 = 21.45106, shape2 = 399.158)
lhs_sample[,38] <- qlnorm(A[,38], meanlog = -4.108138,  sdlog = 0.5998406)
lhs_sample[,39] <- qbeta(A[,39],  shape1 = 393.9252, shape2 = 547.0268)
lhs_sample[,40] <- qbeta(A[,40],  shape1 = 17875.92, shape2 = 42754.63)
lhs_sample[,41] <- qbeta(A[,41],  shape1 = 799.2911, shape2 = 3306.976)
lhs_sample[,42] <- qbeta(A[,42],  shape1 = 77.56885, shape2 = 836.687)
lhs_sample[,43] <- qbeta(A[,43],  shape1 = 7.605944, shape2 = 1074.944)
lhs_sample[,44] <- qlnorm(A[,44], meanlog = -2.145581,  sdlog = 0.1815621)
lhs_sample[,45] <- qlnorm(A[,45], meanlog = -0.5430045, sdlog = 0.154684)
lhs_sample[,46] <- qlnorm(A[,46], meanlog = -2.262311,  sdlog = 0.4746817)
lhs_sample[,47] <- qlnorm(A[,47], meanlog = -1.148854,  sdlog = 0.1815042)
lhs_sample[,48] <- qlnorm(A[,48], meanlog = -0.6931472, sdlog = 0.0511915)
lhs_sample[,49] <- qlnorm(A[,49], meanlog = 0,          sdlog = 0.0511915)
lhs_sample[,50] <- qlnorm(A[,50], meanlog = 0.6931472,  sdlog = 0.08380206)
lhs_sample[,51] <- qbeta(A[,51],  shape1 = 164.6044, shape2 = 618335.4)

colnames(lhs_sample) <- c(
  "p_sae_vacc_u65", "p_sae_vacc_65", "p_death_vacc_u65", "p_death_vacc_65",
  "p_sae_nat_11", "p_sae_nat_17", "p_sae_nat_64", "p_sae_nat_65",
  "p_death_nat_11", "p_death_nat_17", "p_death_nat_64", "p_death_nat_65",
  "ve", "ar_small", "ar_med", "ar_large", "epi_months",
  "trav_7d", "trav_14d", "trav_30d", "trav_90d",
  "symp_asia", "symp_africa", "symp_america", "symp_overall", "fatal_hosp", "hosp", "lt",
  "le_lost_1_11", "le_lost_12_17", "le_lost_18_64", "le_lost_65",
  "dw_chronic", "dur_chronic", "dw_hosp", "dur_acute", "dw_nonhosp", "dur_nonhosp",
  "acute", "subac", "chr6m", "chr12m", "chr30m",
  "dw_chronic_mild", "dw_chronic_severe", "dur_subac", "dw_subac",
  "dur_6m", "dur_12m", "dur_30m", "fatal_nonhosp"
)
lhs_sample <- as.data.frame(lhs_sample)

compute_daly_one <- function(age_group, deaths_10k = 0, hosp_10k = 0,
                             nonhosp_symp_10k = 0, symp_10k = NULL,
                             sae_10k = 0, deaths_sae_10k = 0, draw_pars) {

  life_expectancy <- dplyr::case_when(
    age_group == "1-11"  ~ draw_pars$le_lost_1_11,
    age_group == "12-17" ~ draw_pars$le_lost_12_17,
    age_group == "18-64" ~ draw_pars$le_lost_18_64,
    age_group == "65+"   ~ draw_pars$le_lost_65,
    TRUE ~ NA_real_
  )

  p_acute  <- draw_pars$acute
  p_subac  <- draw_pars$subac
  p_chr6m  <- draw_pars$chr6m
  p_chr12m <- draw_pars$chr12m
  p_chr30m <- draw_pars$chr30m
  total_p <- p_acute + p_subac + p_chr6m + p_chr12m + p_chr30m
  p_acute  <- p_acute  / total_p
  p_subac  <- p_subac  / total_p
  p_chr6m  <- p_chr6m  / total_p
  p_chr12m <- p_chr12m / total_p
  p_chr30m <- p_chr30m / total_p

  dw_hosp    <- draw_pars$dw_hosp
  dw_nonhosp <- draw_pars$dw_nonhosp
  dw_subac   <- draw_pars$dw_subac
  dw_chronic <- draw_pars$dw_chronic
  dur_acute   <- draw_pars$dur_acute
  dur_nonhosp <- draw_pars$dur_nonhosp
  dur_subac   <- draw_pars$dur_subac
  dur_6m      <- draw_pars$dur_6m
  dur_12m     <- draw_pars$dur_12m
  dur_30m     <- draw_pars$dur_30m

  yll_dz <- deaths_10k * life_expectancy
  yld_dz_hosp_acute    <- hosp_10k * dw_hosp    * dur_acute
  yld_dz_nonhosp_acute <- nonhosp_symp_10k * dw_nonhosp * dur_nonhosp
  if (is.null(symp_10k)) symp_10k <- hosp_10k + nonhosp_symp_10k
  dz_subac_10k <- symp_10k * p_subac
  dz_6m_10k    <- symp_10k * p_chr6m
  dz_12m_10k   <- symp_10k * p_chr12m
  dz_30m_10k   <- symp_10k * p_chr30m
  yld_dz_subac <- dz_subac_10k * dw_subac   * dur_subac
  yld_dz_6m    <- dz_6m_10k    * dw_chronic * dur_6m
  yld_dz_12m   <- dz_12m_10k   * dw_chronic * dur_12m
  yld_dz_30m   <- dz_30m_10k   * dw_chronic * dur_30m
  yld_dz <- yld_dz_hosp_acute + yld_dz_nonhosp_acute + yld_dz_subac + yld_dz_6m + yld_dz_12m + yld_dz_30m
  daly_dz <- yll_dz + yld_dz

  yll_sae <- deaths_sae_10k * life_expectancy
  sae_acute_10k  <- sae_10k * p_acute
  sae_subac_10k  <- sae_10k * p_subac
  sae_6m_10k     <- sae_10k * p_chr6m
  sae_12m_10k    <- sae_10k * p_chr12m
  sae_30m_10k    <- sae_10k * p_chr30m
  yld_sae_acute <- sae_acute_10k * dw_hosp    * dur_acute
  yld_sae_subac <- sae_subac_10k * dw_subac   * dur_subac
  yld_sae_6m    <- sae_6m_10k    * dw_chronic * dur_6m
  yld_sae_12m   <- sae_12m_10k   * dw_chronic * dur_12m
  yld_sae_30m   <- sae_30m_10k   * dw_chronic * dur_30m
  yld_sae <- yld_sae_acute + yld_sae_subac + yld_sae_6m + yld_sae_12m + yld_sae_30m
  daly_sae <- yll_sae + yld_sae

  list(daly_dz = daly_dz, yll_dz = yll_dz, yld_dz = yld_dz,
       daly_sae = daly_sae, yll_sae = yll_sae, yld_sae = yld_sae)
}


## ---- Age-specific extension + compute_daly_one_age_specific() --------
## Verbatim copy of 02_Scripts/02_setup_age_props.R (self-contained given
## lhs_sample + compute_daly_one() above).

w_u40_ibge2022 <- 92.6 / (92.6 + 70.5)
w_u40 <- w_u40_ibge2022
w_o40 <- 1 - w_u40

o40_beta_params <- list(
  acute  = c(shape1 = 388.2874, shape2 = 742.4988),
  subac  = c(shape1 = 2487.085, shape2 = 6450.815),
  chr6m  = c(shape1 = 5998.872, shape2 = 21654.86),
  chr12m = c(shape1 = 184.6674, shape2 = 1216.058),
  chr30m = c(shape1 = 15.74654, shape2 = 516.3419)
)

lhs_sample$acute_u40  <- lhs_sample$acute
lhs_sample$subac_u40  <- lhs_sample$subac
lhs_sample$chr6m_u40  <- lhs_sample$chr6m
lhs_sample$chr12m_u40 <- lhs_sample$chr12m
lhs_sample$chr30m_u40 <- lhs_sample$chr30m

set.seed(20260420)
A_o40 <- lhs::randomLHS(n = nrow(lhs_sample), k = 5)
lhs_sample$acute_o40  <- qbeta(A_o40[,1], shape1 = o40_beta_params$acute["shape1"],  shape2 = o40_beta_params$acute["shape2"])
lhs_sample$subac_o40  <- qbeta(A_o40[,2], shape1 = o40_beta_params$subac["shape1"],  shape2 = o40_beta_params$subac["shape2"])
lhs_sample$chr6m_o40  <- qbeta(A_o40[,3], shape1 = o40_beta_params$chr6m["shape1"],  shape2 = o40_beta_params$chr6m["shape2"])
lhs_sample$chr12m_o40 <- qbeta(A_o40[,4], shape1 = o40_beta_params$chr12m["shape1"], shape2 = o40_beta_params$chr12m["shape2"])
lhs_sample$chr30m_o40 <- qbeta(A_o40[,5], shape1 = o40_beta_params$chr30m["shape1"], shape2 = o40_beta_params$chr30m["shape2"])

compute_daly_one_age_specific <- function(age_group, deaths_10k = 0, hosp_10k = 0,
                                          nonhosp_symp_10k = 0, symp_10k = NULL,
                                          sae_10k = 0, deaths_sae_10k = 0,
                                          draw_pars, draw_id = NULL,
                                          w_u40_local = w_u40, w_o40_local = w_o40) {
  needed <- c("acute_u40","acute_o40","subac_u40","subac_o40","chr6m_u40","chr6m_o40",
              "chr12m_u40","chr12m_o40","chr30m_u40","chr30m_o40")
  missing_cols <- setdiff(needed, names(draw_pars))
  if (length(missing_cols) > 0) {
    if (is.null(draw_id)) stop("compute_daly_one_age_specific(): missing age-specific columns and no draw_id.")
    for (col in missing_cols) draw_pars[[col]] <- lhs_sample[[col]][draw_id]
  }

  blend <- function(u, o) u * w_u40_local + o * w_o40_local
  eff <- switch(as.character(age_group),
    "1-11"  = list(acute = draw_pars$acute_u40, subac = draw_pars$subac_u40, chr6m = draw_pars$chr6m_u40, chr12m = draw_pars$chr12m_u40, chr30m = draw_pars$chr30m_u40),
    "12-17" = list(acute = draw_pars$acute_u40, subac = draw_pars$subac_u40, chr6m = draw_pars$chr6m_u40, chr12m = draw_pars$chr12m_u40, chr30m = draw_pars$chr30m_u40),
    "18-64" = list(acute = blend(draw_pars$acute_u40, draw_pars$acute_o40), subac = blend(draw_pars$subac_u40, draw_pars$subac_o40),
                  chr6m = blend(draw_pars$chr6m_u40, draw_pars$chr6m_o40), chr12m = blend(draw_pars$chr12m_u40, draw_pars$chr12m_o40),
                  chr30m = blend(draw_pars$chr30m_u40, draw_pars$chr30m_o40)),
    "65+"   = list(acute = draw_pars$acute_o40, subac = draw_pars$subac_o40, chr6m = draw_pars$chr6m_o40, chr12m = draw_pars$chr12m_o40, chr30m = draw_pars$chr30m_o40),
    stop("Unknown age_group: ", age_group)
  )
  draw_pars$acute  <- eff$acute
  draw_pars$subac  <- eff$subac
  draw_pars$chr6m  <- eff$chr6m
  draw_pars$chr12m <- eff$chr12m
  draw_pars$chr30m <- eff$chr30m

  compute_daly_one(age_group = age_group, deaths_10k = deaths_10k, hosp_10k = hosp_10k,
                   nonhosp_symp_10k = nonhosp_symp_10k, symp_10k = symp_10k,
                   sae_10k = sae_10k, deaths_sae_10k = deaths_sae_10k, draw_pars = draw_pars)
}

daly_pars_true <- lhs_sample %>%
  dplyr::mutate(draw_id = dplyr::row_number()) %>%
  dplyr::select(
    draw_id,
    le_lost_1_11, le_lost_12_17, le_lost_18_64, le_lost_65,
    dw_hosp, dw_nonhosp, dw_subac, dw_chronic,
    dur_acute, dur_nonhosp, dur_subac, dur_6m, dur_12m, dur_30m,
    acute, subac, chr6m, chr12m, chr30m
  )

risk_draw_df <- bind_rows(
  tibble(draw_id = 1:nrow(lhs_sample), AgeCat = "1-11",  p_sae_vacc_base = lhs_sample[, "p_sae_vacc_u65"], p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]),
  tibble(draw_id = 1:nrow(lhs_sample), AgeCat = "12-17", p_sae_vacc_base = lhs_sample[, "p_sae_vacc_u65"], p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]),
  tibble(draw_id = 1:nrow(lhs_sample), AgeCat = "18-64", p_sae_vacc_base = lhs_sample[, "p_sae_vacc_u65"], p_death_vacc_base = lhs_sample[, "p_death_vacc_u65"]),
  tibble(draw_id = 1:nrow(lhs_sample), AgeCat = "65+",   p_sae_vacc_base = lhs_sample[, "p_sae_vacc_65"],  p_death_vacc_base = lhs_sample[, "p_death_vacc_65"])
)


## ---- Benefit-side draw helpers -----------------------------------------
## Verbatim copies of 03_brazil_all_draws_ori_v3.R Sections 03-06 helpers.

# NOTE (2026-08): now identity passthroughs -- pre_list/post_arr are already
# TRUE burden (p_sym-only, no rho) out of the m3-refit simulator, so dividing
# by rho here would re-inflate by ~1/rho (~7-10x). Same fix as applied to the
# identically-named helpers in 03_brazil_all_draws_ori_v3.R -- see chat
# record, 2026-08.
scale_symp_by_rho_pre <- function(pre_list, rho_vec) {
  pre_list
}
scale_symp_by_rho_post <- function(post_arr, rho_vec) {
  post_arr
}

make_hosp_draws <- function(symp_list, hosp_rate) {
  lapply(symp_list, function(m) sweep(m, 1, hosp_rate, `*`))
}
calc_total_hosp_draws <- function(pre_list, post_arr, hosp_rate) {
  pre_arr <- simplify2array(make_hosp_draws(pre_list, hosp_rate))
  if (length(dim(pre_arr)) == 2) pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  post_list <- lapply(seq_len(dim(post_arr)[3]), function(i) post_arr[,,i])
  post_arr2 <- simplify2array(make_hosp_draws(post_list, hosp_rate))
  if (length(dim(post_arr2)) == 2) post_arr2 <- array(post_arr2, dim = c(dim(post_arr2), 1))
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  tibble(draw_id = seq_len(n_draws),
         total_pre  = apply(pre_arr[,,1:n_draws, drop=FALSE],  3, sum, na.rm = TRUE),
         total_post = apply(post_arr2[,,1:n_draws, drop=FALSE], 3, sum, na.rm = TRUE))
}
calc_total_hosp_draws_rho <- function(pre_list, post_arr, hosp_rate, rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]; post_arr <- post_arr[,,1:n_draws, drop = FALSE]; rho_vec <- rho_vec[1:n_draws]
  calc_total_hosp_draws(scale_symp_by_rho_pre(pre_list, rho_vec), scale_symp_by_rho_post(post_arr, rho_vec), hosp_rate)
}

make_fatal_hosp_draws <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate) {
  lapply(symp_list, function(m) {
    hospitalised <- sweep(m, 1, hosp_rate, `*`)
    non_hospitalised <- m - hospitalised
    fatal <- sweep(hospitalised, 1, fatal_rate, `*`) + sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    list(hosp = hospitalised, fatal = fatal)
  })
}
calc_total_fatal_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate) {
  pre_conv <- make_fatal_hosp_draws(pre_list, hosp_rate, fatal_rate, nh_fatal_rate)
  pre_arr <- simplify2array(lapply(pre_conv, `[[`, "fatal"))
  if (length(dim(pre_arr)) == 2) pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[,,i])
  post_conv <- make_fatal_hosp_draws(post_list, hosp_rate, fatal_rate, nh_fatal_rate)
  post_arr2 <- simplify2array(lapply(post_conv, `[[`, "fatal"))
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  tibble(draw_id = 1:n_draws,
         total_pre  = apply(pre_arr[,,1:n_draws, drop=FALSE],  3, sum, na.rm = TRUE),
         total_post = apply(post_arr2[,,1:n_draws, drop=FALSE], 3, sum, na.rm = TRUE))
}
calc_total_fatal_draws_rho <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate, rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]; post_arr <- post_arr[,,1:n_draws, drop=FALSE]; rho_vec <- rho_vec[1:n_draws]
  calc_total_fatal_draws(scale_symp_by_rho_pre(pre_list, rho_vec), scale_symp_by_rho_post(post_arr, rho_vec), hosp_rate, fatal_rate, nh_fatal_rate)
}

make_daly_draws_age <- function(symp_list, hosp_rate, fatal_rate, nh_fatal_rate,
                                dw_hosp, dw_nonhosp, dw_chronic, dur_acute, dur_subacute,
                                dur_chronic, dw_subacute, subac_prop_vec, chr_prop_vec, le_left_vec) {
  lapply(symp_list, function(m) {
    hospitalised <- sweep(m, 1, hosp_rate, `*`)
    non_hospitalised <- m - hospitalised
    fatal <- sweep(hospitalised, 1, fatal_rate, `*`) + sweep(non_hospitalised, 1, nh_fatal_rate, `*`)
    yld_acute <- (hospitalised * dw_hosp * dur_acute) + (non_hospitalised * dw_nonhosp * dur_acute)
    yld_subacute <- sweep(hospitalised, 1, subac_prop_vec, `*`) * dw_subacute * dur_subacute +
      sweep(non_hospitalised, 1, chr_prop_vec, `*`) * dw_subacute * dur_subacute
    yld_chronic <- sweep(hospitalised, 1, chr_prop_vec, `*`) * dw_chronic * dur_chronic +
      sweep(non_hospitalised, 1, chr_prop_vec, `*`) * dw_chronic * dur_chronic
    yld_total <- yld_acute + yld_subacute + yld_chronic
    yll <- sweep(fatal, 1, le_left_vec, `*`)
    list(daly_tot = yld_total + yll)
  })
}
calc_total_daly_draws <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate,
                                  dw_hosp, dw_nonhosp, dw_chronic, dur_acute, dur_subacute, dur_chronic,
                                  dw_subacute, subac_prop, chr_prop, le_left_vec) {
  pre_conv <- make_daly_draws_age(pre_list, hosp_rate, fatal_rate, nh_fatal_rate, dw_hosp, dw_nonhosp, dw_chronic,
                                  dur_acute, dur_subacute, dur_chronic, dw_subacute, subac_prop, chr_prop, le_left_vec)
  pre_arr <- simplify2array(lapply(pre_conv, `[[`, "daly_tot"))
  if (length(dim(pre_arr)) == 2) pre_arr <- array(pre_arr, dim = c(dim(pre_arr), 1))
  post_list <- lapply(seq(dim(post_arr)[3]), function(i) post_arr[,,i])
  post_conv <- make_daly_draws_age(post_list, hosp_rate, fatal_rate, nh_fatal_rate, dw_hosp, dw_nonhosp, dw_chronic,
                                   dur_acute, dur_subacute, dur_chronic, dw_subacute, subac_prop, chr_prop, le_left_vec)
  post_arr2 <- simplify2array(lapply(post_conv, `[[`, "daly_tot"))
  n_draws <- min(dim(pre_arr)[3], dim(post_arr2)[3])
  tibble(draw_id = 1:n_draws,
         total_pre  = apply(pre_arr[,,1:n_draws, drop=FALSE],  3, sum, na.rm = TRUE),
         total_post = apply(post_arr2[,,1:n_draws, drop=FALSE], 3, sum, na.rm = TRUE))
}
calc_total_daly_draws_rho <- function(pre_list, post_arr, hosp_rate, fatal_rate, nh_fatal_rate,
                                      dw_hosp, dw_nonhosp, dw_chronic, dur_acute, dur_subacute, dur_chronic,
                                      dw_subacute, subac_prop, chr_prop, le_left_vec, rho_vec) {
  n_draws <- min(length(pre_list), dim(post_arr)[3], length(rho_vec))
  pre_list <- pre_list[1:n_draws]; post_arr <- post_arr[,,1:n_draws, drop=FALSE]; rho_vec <- rho_vec[1:n_draws]
  calc_total_daly_draws(scale_symp_by_rho_pre(pre_list, rho_vec), scale_symp_by_rho_post(post_arr, rho_vec),
                        hosp_rate, fatal_rate, nh_fatal_rate, dw_hosp, dw_nonhosp, dw_chronic,
                        dur_acute, dur_subacute, dur_chronic, dw_subacute, subac_prop, chr_prop, le_left_vec)
}

age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17),
                mean(18:19), mean(20:24), mean(25:29), mean(30:34), mean(35:39),
                mean(40:44), mean(45:49), mean(50:54), mean(55:59), mean(60:64),
                mean(65:69), mean(70:74), mean(75:79), mean(80:84), mean(85:89))
le_by_age <- function(age_numeric) {
  if (age_numeric <= 1) return(quantile(le_sample$le_1, 0.5))
  if (age_numeric < 20) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 30) return(quantile(le_sample$le_2, 0.5))
  if (age_numeric < 40) return(quantile(le_sample$le_3, 0.5))
  if (age_numeric < 50) return(quantile(le_sample$le_4, 0.5))
  if (age_numeric < 60) return(quantile(le_sample$le_5, 0.5))
  if (age_numeric < 70) return(quantile(le_sample$le_6, 0.5))
  if (age_numeric < 80) return(quantile(le_sample$le_7, 0.5))
  return(quantile(le_sample$le_8, 0.5))
}
le_sample <- readRDS("00_Data/0_2_Processed/le_sample.RDS")
le_left_vec <- sapply(age_groups, le_by_age)

subac_u40 <- stats::median(lhs_sample_young$subac, na.rm = TRUE)
chr_u40   <- stats::median(lhs_sample_young$chr6m,  na.rm = TRUE) +
  stats::median(lhs_sample_young$chr12m, na.rm = TRUE) + stats::median(lhs_sample_young$chr30m, na.rm = TRUE)
subac_o40 <- stats::median(lhs_sample$subac_o40, na.rm = TRUE)
chr_o40   <- stats::median(lhs_sample$chr6m_o40,  na.rm = TRUE) +
  stats::median(lhs_sample$chr12m_o40, na.rm = TRUE) + stats::median(lhs_sample$chr30m_o40, na.rm = TRUE)
u40_cutoff_bin <- 10L
idx_o40 <- seq.int(u40_cutoff_bin + 1L, 20L)
subac_prop_vec <- rep(subac_u40, 20L); subac_prop_vec[idx_o40] <- subac_o40
chr_prop_vec   <- rep(chr_u40,   20L); chr_prop_vec[idx_o40]   <- chr_o40

calc_q_seromix_for_scenarios_agecat <- function(sim_region_ve_cov, age_map) {
  n_scenarios <- length(sim_region_ve_cov)
  out <- vector("list", n_scenarios)
  for (sc in seq_len(n_scenarios)) {
    raw_alloc_array <- sim_region_ve_cov[[sc]]$sim_result$raw_allocation_array
    vacc_to_S_array <- sim_region_ve_cov[[sc]]$sim_result$vacc_to_S_array
    n_draws <- dim(raw_alloc_array)[3]
    sc_df <- vector("list", n_draws)
    for (d in seq_len(n_draws)) {
      df <- data.frame(
        age_index = seq_len(dim(raw_alloc_array)[1]),
        total_vacc_age = rowSums(raw_alloc_array[,,d, drop = FALSE], na.rm = TRUE),
        seroneg_vacc_age = rowSums(vacc_to_S_array[,,d, drop = FALSE], na.rm = TRUE)
      ) %>%
        dplyr::mutate(seropos_vacc_age = pmax(0, total_vacc_age - seroneg_vacc_age)) %>%
        dplyr::left_join(age_map, by = "age_index") %>%
        dplyr::group_by(AgeCat) %>%
        dplyr::summarise(total_vacc_age = sum(total_vacc_age, na.rm=TRUE),
                         seroneg_vacc_age = sum(seroneg_vacc_age, na.rm=TRUE),
                         seropos_vacc_age = sum(seropos_vacc_age, na.rm=TRUE), .groups = "drop") %>%
        dplyr::mutate(Scenario = sc, draw_id = d,
                     q_seroneg_vacc = ifelse(total_vacc_age > 0, seroneg_vacc_age / total_vacc_age, NA_real_),
                     q_seropos_vacc = ifelse(total_vacc_age > 0, seropos_vacc_age / total_vacc_age, NA_real_))
      sc_df[[d]] <- df
    }
    out[[sc]] <- dplyr::bind_rows(sc_df)
  }
  dplyr::bind_rows(out)
}

# Generalized (no coverage_keep filter -- Coverage/week already carried by df_true).
make_averted_draws_true <- function(df_true, outcome_name) {
  df_true %>%
    dplyr::mutate(
      Scenario = as.integer(Scenario),
      outcome  = outcome_name,
      baseline = total_pre,
      post     = total_post,
      averted  = total_pre - total_post
    ) %>%
    dplyr::select(draw_id, Region, VE, Coverage, week, Scenario, outcome, baseline, post, averted)
}


## ---- Per-week driver: Sections 02, 04-08B, generalized over Coverage/week ----

all_weeks_brr <- purrr::map_dfr(WEEKS, function(week_val) {

  message(sprintf("[finite_weeksweep_brr] week = %d -- loading postsim", week_val))
  load(sprintf("00_Data/0_2_Processed/postsim_vc_ixchiq_model_finite_week%d%s.RData", week_val, FILE_SUFFIX))
  # -> postsim_vc_ixchiq_model (region -> VE -> coverage -> postsim_all_ui output)

  hosp_fatal_env <- new.env()
  load("00_Data/0_2_Processed/chikv_fatal_hosp_rate.RData", envir = hosp_fatal_env)
  hosp <- hosp_fatal_env$hosp; fatal <- hosp_fatal_env$fatal; nh_fatal <- hosp_fatal_env$nh_fatal

  message(sprintf("[finite_weeksweep_brr] week = %d -- all_draws_hosp/fatal/daly_true", week_val))
  all_draws_hosp_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
    pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
    rho_pool_region <- as.numeric(posterior_list_finite[[region_name]]$rho)
    rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
    purrr::imap_dfr(region_list, function(ve_list, ve_name) {
      purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
        purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
          post_arr <- scen$sim_result$age_array_raw_symp
          n_draws <- min(length(pre_list), dim(post_arr)[3])
          rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
          calc_total_hosp_draws_rho(pre_list, post_arr, hosp, rho_vec) %>%
            dplyr::mutate(Region = region_name, VE = ve_name, Coverage = cov_name, week = week_val, Scenario = scen_id)
        })
      })
    })
  })

  all_draws_fatal_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
    pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
    rho_pool_region <- as.numeric(posterior_list_finite[[region_name]]$rho)
    rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
    purrr::imap_dfr(region_list, function(ve_list, ve_name) {
      purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
        purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
          post_arr <- scen$sim_result$age_array_raw_symp
          n_draws <- min(length(pre_list), dim(post_arr)[3])
          rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
          calc_total_fatal_draws_rho(pre_list, post_arr, hosp, fatal, nh_fatal, rho_vec) %>%
            dplyr::mutate(Region = region_name, VE = ve_name, Coverage = cov_name, week = week_val, Scenario = scen_id)
        })
      })
    })
  })

  all_draws_daly_true <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
    pre_list <- preui_all[[region_name]]$sim_results_list_rawsymp
    rho_pool_region <- as.numeric(posterior_list_finite[[region_name]]$rho)
    rho_pool_region <- rho_pool_region[is.finite(rho_pool_region) & rho_pool_region > 0]
    purrr::imap_dfr(region_list, function(ve_list, ve_name) {
      purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
        purrr::imap_dfr(cov_list$scenario_result, function(scen, scen_id) {
          post_arr <- scen$sim_result$age_array_raw_symp
          n_draws <- min(length(pre_list), dim(post_arr)[3])
          rho_vec <- sample(rho_pool_region, size = n_draws, replace = TRUE)
          calc_total_daly_draws_rho(
            pre_list, post_arr, hosp, fatal, nh_fatal,
            dw_hosp = quantile(lhs_sample_young$dw_hosp, 0.5), dw_nonhosp = quantile(lhs_sample_young$dw_nonhosp, 0.5),
            dw_chronic = quantile(lhs_sample_young$dw_chronic, 0.5), dur_acute = quantile(lhs_sample_young$dur_acute, 0.5),
            dur_subacute = quantile(lhs_sample_young$dur_subac, 0.5), dur_chronic = quantile(lhs_sample_young$dur_chronic, 0.5),
            dw_subacute = quantile(lhs_sample_young$dw_subac, 0.5),
            subac_prop = subac_prop_vec, chr_prop = chr_prop_vec, le_left_vec = le_left_vec, rho_vec = rho_vec
          ) %>%
            dplyr::mutate(Region = region_name, VE = ve_name, Coverage = cov_name, week = week_val, Scenario = scen_id)
        })
      })
    })
  })

  # SAE outcome = hosp + fatal combined (matches 03_brazil's Step 7.1).
  all_draws_sae_true <- all_draws_hosp_true %>%
    dplyr::rename(total_pre_hosp = total_pre, total_post_hosp = total_post) %>%
    dplyr::inner_join(
      all_draws_fatal_true %>% dplyr::rename(total_pre_fatal = total_pre, total_post_fatal = total_post),
      by = c("draw_id", "Region", "VE", "Coverage", "week", "Scenario")
    ) %>%
    dplyr::mutate(total_pre = total_pre_hosp + total_pre_fatal, total_post = total_post_hosp + total_post_fatal) %>%
    dplyr::select(draw_id, total_pre, total_post, Region, VE, Coverage, week, Scenario)

  message(sprintf("[finite_weeksweep_brr] week = %d -- q_all_regions (serostatus decomposition)", week_val))
  q_all_regions <- purrr::imap_dfr(postsim_vc_ixchiq_model, function(region_list, region_name) {
    purrr::imap_dfr(region_list, function(ve_list, ve_name) {
      purrr::imap_dfr(ve_list, function(cov_list, cov_name) {
        calc_q_seromix_for_scenarios_agecat(cov_list$scenario_result, age_map) %>%
          dplyr::mutate(Region = region_name, VE = ve_name, Coverage = cov_name, week = week_val) %>%
          dplyr::select(Region, VE, Coverage, week, Scenario, draw_id, AgeCat,
                       total_vacc_age, seroneg_vacc_age, seropos_vacc_age, q_seroneg_vacc, q_seropos_vacc)
      })
    })
  })

  message(sprintf("[finite_weeksweep_brr] week = %d -- benefit + risk join", week_val))
  benefit_base_true <- bind_rows(
    make_averted_draws_true(all_draws_daly_true,  "DALY"),
    make_averted_draws_true(all_draws_sae_true,   "SAE"),
    make_averted_draws_true(all_draws_fatal_true, "Death")
  ) %>%
    dplyr::mutate(Scenario = as.integer(Scenario), AgeCat = map_scenario_agecat_int(Scenario)) %>%
    dplyr::left_join(daly_pars_true, by = "draw_id")

  tot_vacc_map_true <- all_weeks_nnv %>%
    dplyr::filter(week == week_val) %>%
    dplyr::mutate(
      AgeCat = case_when(
        AgeGroup %in% 2:4   ~ "1-11",
        AgeGroup == 5       ~ "12-17",
        AgeGroup %in% 6:15  ~ "18-64",
        AgeGroup %in% 16:20 ~ "65+",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(AgeCat)) %>%
    dplyr::group_by(scenario, region, AgeCat, VE, VC, week) %>%
    dplyr::summarise(tot_vacc_grp = sum(tot_vacc, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      target = case_when(
        scenario == "Scenario_1" & AgeCat == "1-11"  ~ 1L,
        scenario == "Scenario_2" & AgeCat == "12-17" ~ 1L,
        scenario == "Scenario_3" & AgeCat == "18-64" ~ 1L,
        scenario == "Scenario_4" & AgeCat == "65+"   ~ 1L,
        TRUE ~ 0L
      )
    ) %>%
    dplyr::filter(target == 1L) %>%
    dplyr::transmute(
      scenario = as.integer(gsub("Scenario_", "", scenario)),
      AgeCat, VE, Coverage = VC, week, tot_vacc_grp, Region = region
    )

  rr_vals_local <- rr_vals

  risk_components_all <- q_all_regions %>%
    tidyr::crossing(RR_seropos = rr_vals_local) %>%
    dplyr::left_join(risk_draw_df, by = c("draw_id", "AgeCat")) %>%
    dplyr::mutate(
      p_sae_vacc_seroneg = p_sae_vacc_base, p_death_vacc_seroneg = p_death_vacc_base,
      p_sae_vacc_seropos = p_sae_vacc_base * RR_seropos, p_death_vacc_seropos = p_death_vacc_base * RR_seropos,
      p_sae_vacc_seroneg_contrib = dplyr::if_else(total_vacc_age > 0, q_seroneg_vacc * p_sae_vacc_seroneg, 0),
      p_death_vacc_seroneg_contrib = dplyr::if_else(total_vacc_age > 0, q_seroneg_vacc * p_death_vacc_seroneg, 0),
      p_sae_vacc_seropos_contrib = dplyr::if_else(total_vacc_age > 0, q_seropos_vacc * p_sae_vacc_seropos, 0),
      p_death_vacc_seropos_contrib = dplyr::if_else(total_vacc_age > 0, q_seropos_vacc * p_death_vacc_seropos, 0),
      p_sae_vacc_adj = p_sae_vacc_seroneg_contrib + p_sae_vacc_seropos_contrib,
      p_death_vacc_adj = p_death_vacc_seroneg_contrib + p_death_vacc_seropos_contrib
    ) %>%
    dplyr::mutate(
      sae_10k_base   = 1e4 * p_sae_vacc_base + 1e4 * p_death_vacc_base,
      death_10k_base = 1e4 * p_death_vacc_base,
      sae_10k_seroneg   = 1e4 * p_sae_vacc_seroneg_contrib + 1e4 * p_death_vacc_seroneg_contrib,
      death_10k_seroneg = 1e4 * p_death_vacc_seroneg_contrib,
      sae_10k_seropos   = 1e4 * p_sae_vacc_seropos_contrib + 1e4 * p_death_vacc_seropos_contrib,
      death_10k_seropos = 1e4 * p_death_vacc_seropos_contrib,
      sae_10k_adj   = sae_10k_seroneg + sae_10k_seropos,
      death_10k_adj = death_10k_seroneg + death_10k_seropos
    )

  risk_join_all <- risk_components_all %>%
    dplyr::select(Region, VE, Coverage, week, Scenario, draw_id, AgeCat, RR_seropos,
                 total_vacc_age, sae_10k_base, death_10k_base, sae_10k_adj, death_10k_adj,
                 sae_10k_seroneg, death_10k_seroneg, sae_10k_seropos, death_10k_seropos)

  benefit_draw_df <- benefit_base_true %>%
    dplyr::left_join(tot_vacc_map_true, by = c("Region", "Scenario" = "scenario", "VE", "Coverage", "week", "AgeCat")) %>%
    dplyr::mutate(
      Scenario = as.integer(Scenario),
      averted_10k  = (averted / tot_vacc_grp) * 1e4
    )
  benefit_draw_df <- tidyr::crossing(benefit_draw_df, RR_seropos = rr_vals_local)

  draw_level_xy_serostatus <- benefit_draw_df %>%
    dplyr::left_join(risk_join_all, by = c("Region", "VE", "Coverage", "week", "Scenario", "draw_id", "AgeCat", "RR_seropos"))

  # ---- Vectorized (not rowwise -- "very slow" per the original's own
  # comment) DALY-of-vaccine-SAE computation, split by AgeCat since
  # compute_daly_one_age_specific()'s switch() only takes a scalar
  # age_group, but every OTHER input is already a plain vector -- so each
  # AgeCat's rows can be computed in one vectorized call instead of row by row.
  message(sprintf("[finite_weeksweep_brr] week = %d -- vaccine-SAE DALY (vectorized by AgeCat)", week_val))
  draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
    dplyr::filter(outcome == "DALY", !is.na(AgeCat)) %>%
    dplyr::group_split(AgeCat) %>%
    purrr::map_dfr(function(df) {
      ag <- df$AgeCat[1]
      dp <- list(
        le_lost_1_11 = df$le_lost_1_11, le_lost_12_17 = df$le_lost_12_17,
        le_lost_18_64 = df$le_lost_18_64, le_lost_65 = df$le_lost_65,
        dw_hosp = df$dw_hosp, dw_nonhosp = df$dw_nonhosp,
        dw_subac = df$dw_subac, dw_chronic = df$dw_chronic,
        dur_acute = df$dur_acute, dur_nonhosp = df$dur_nonhosp, dur_subac = df$dur_subac,
        dur_6m = df$dur_6m, dur_12m = df$dur_12m, dur_30m = df$dur_30m,
        acute = df$acute, subac = df$subac, chr6m = df$chr6m, chr12m = df$chr12m, chr30m = df$chr30m
      )
      base_res    <- compute_daly_one_age_specific(ag, sae_10k = df$sae_10k_base,    deaths_sae_10k = df$death_10k_base,    draw_id = df$draw_id, draw_pars = dp)
      seroneg_res <- compute_daly_one_age_specific(ag, sae_10k = df$sae_10k_seroneg, deaths_sae_10k = df$death_10k_seroneg, draw_id = df$draw_id, draw_pars = dp)
      seropos_res <- compute_daly_one_age_specific(ag, sae_10k = df$sae_10k_seropos, deaths_sae_10k = df$death_10k_seropos, draw_id = df$draw_id, draw_pars = dp)
      df$daly_10k_base    <- base_res$daly_sae
      df$daly_10k_seroneg <- seroneg_res$daly_sae
      df$daly_10k_seropos <- seropos_res$daly_sae
      df$daly_10k_adj     <- df$daly_10k_seroneg + df$daly_10k_seropos
      df
    }) %>%
    dplyr::bind_rows(
      draw_level_xy_serostatus %>%
        dplyr::filter(outcome != "DALY" | is.na(AgeCat)) %>%
        dplyr::mutate(daly_10k_base = NA_real_, daly_10k_seroneg = NA_real_, daly_10k_seropos = NA_real_, daly_10k_adj = NA_real_)
    )

  draw_level_xy_serostatus <- draw_level_xy_serostatus %>%
    dplyr::mutate(
      VE_label = factor(VE, levels = c("VE0", "VE98.9"), labels = c("Disease blocking only", "Disease and infection blocking")),
      setting = unname(setting_key[Region])
    ) %>%
    dplyr::mutate(
      x_10k_base = dplyr::case_when(outcome == "SAE" ~ sae_10k_base, outcome == "Death" ~ death_10k_base, outcome == "DALY" ~ daly_10k_base, TRUE ~ NA_real_),
      x_10k_adj  = dplyr::case_when(outcome == "SAE" ~ sae_10k_adj,  outcome == "Death" ~ death_10k_adj,  outcome == "DALY" ~ daly_10k_adj,  TRUE ~ NA_real_),
      brr_base = dplyr::if_else(is.na(x_10k_base) | x_10k_base == 0, NA_real_, averted_10k / x_10k_base),
      brr_adj  = dplyr::if_else(is.na(x_10k_adj)  | x_10k_adj  == 0, NA_real_, averted_10k / x_10k_adj)
    )

  draw_level_xy_serostatus
})

save(all_weeks_brr, file = sprintf("00_Data/0_2_Processed/draw_level_xy_serostatus_finite_weeksweep%s.RData", FILE_SUFFIX))


## ---- Summarize across draws: median/95%UI BRR per (week, coverage, region/setting, scenario, AgeCat, VE, RR_seropos) ----

brr_weeksweep_summary <- all_weeks_brr %>%
  dplyr::mutate(
    brr_base = ifelse(is.infinite(brr_base), NA, brr_base),
    brr_adj  = ifelse(is.infinite(brr_adj),  NA, brr_adj),
    setting  = factor(setting, levels = c("Low", "Moderate", "High"))
  ) %>%
  dplyr::group_by(outcome, Scenario, AgeCat, VE_label, RR_seropos, setting, Coverage, week) %>%
  dplyr::summarise(
    brr_base_med = quantile(brr_base, 0.50,  na.rm = TRUE),
    brr_base_lo  = quantile(brr_base, 0.025, na.rm = TRUE),
    brr_base_hi  = quantile(brr_base, 0.975, na.rm = TRUE),
    brr_adj_med  = quantile(brr_adj,  0.50,  na.rm = TRUE),
    brr_adj_lo   = quantile(brr_adj,  0.025, na.rm = TRUE),
    brr_adj_hi   = quantile(brr_adj,  0.975, na.rm = TRUE),
    # NOTE (2026-08): added -- finite_weeksweep_heatmap_probgt1.R expects
    # these two columns (Pr(BRR>1) across draws per cell) and errored with
    # "Column `brr_base_prob_gt1` not found" without them. See chat record.
    brr_base_prob_gt1 = mean(brr_base > 1, na.rm = TRUE),
    brr_adj_prob_gt1  = mean(brr_adj  > 1, na.rm = TRUE),
    .groups = "drop"
  )

save(brr_weeksweep_summary, file = sprintf("00_Data/0_2_Processed/brr_weeksweep_summary_finite%s.RData", FILE_SUFFIX))

dir.create("02_Outputs/2_2_Tables", showWarnings = FALSE, recursive = TRUE)
writexl::write_xlsx(brr_weeksweep_summary, path = sprintf("02_Outputs/2_2_Tables/brr_weeksweep_summary_finite%s.xlsx", FILE_SUFFIX))

message("[finite_weeksweep_brr] done. brr_weeksweep_summary rows: ", nrow(brr_weeksweep_summary))
message("Columns: outcome, Scenario, AgeCat, VE_label, RR_seropos, setting, Coverage, week, brr_base_med/lo/hi, brr_adj_med/lo/hi, brr_base_prob_gt1, brr_adj_prob_gt1")
message("For the heatmap: filter RR_seropos == 0, pick outcome/AgeCat/VE_label, then x=week, y=Coverage, fill=brr_adj_med (or brr_base_med).")
