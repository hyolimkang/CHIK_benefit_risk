# Unit tests of the infection-risk and benefit-risk functions defined in
# 02_Scripts/01_setup.R. The functions are taken from the script itself
# (load_functions_from), not re-implemented. Expected values follow directly
# from the definitions used in the manuscript methods:
#   cumulative hazard H = -log(1 - AR);  P(infection) = 1 - exp(-H)
#   benefit per 10,000 = 1e4 * AR * p_outcome * VE
#   risk    per 10,000 = 1e4 * p_vaccine_outcome
#   BRR = benefit / risk (undefined when risk = 0)

fn <- load_functions_from("02_Scripts/01_setup.R", c(
  "compute_ar", "compute_outcome", "make_weights_foi", "weekly_to_daily_foi",
  "compute_daily_foi", "compute_ar_travel"))

test_that("compute_ar: annual hazard to attack rate", {
  expect_equal(fn$compute_ar(lambda = 0.1, days = 365), 1 - exp(-0.1))
  expect_equal(fn$compute_ar(lambda = 0.1, days = 0), 0)
  expect_equal(fn$compute_ar(lambda = 0.1, s0 = 0.5, days = 365), 0.5 * (1 - exp(-0.1)))
  ar <- fn$compute_ar(lambda = 0.3, days = c(7, 14, 30, 90))
  expect_true(all(diff(ar) > 0))
  expect_true(all(ar >= 0 & ar <= 1))
})

test_that("make_weights_foi: non-negative, normalised, symmetric", {
  for (L in c(52, 364, 365)) {
    w <- fn$make_weights_foi(L)
    expect_length(w, L)
    expect_equal(sum(w), 1)
    expect_true(all(w >= 0))
    expect_equal(w, rev(w))
  }
})

test_that("weekly_to_daily_foi: 7 days per week, total hazard preserved", {
  wk <- c(0.01, 0.05, 0.02)
  d <- fn$weekly_to_daily_foi(wk)
  expect_length(d, 21)
  expect_equal(sum(d), sum(wk))
  expect_equal(d[1:7], rep(0.01 / 7, 7))
})

test_that("compute_daily_foi: integrates back to the total attack rate", {
  for (ar in c(0.001, 0.05, 0.3)) {
    f <- fn$compute_daily_foi(AR_total = ar, L = 364)
    expect_length(f, 364)
    expect_equal(1 - exp(-sum(f)), ar)
  }
})

test_that("compute_ar_travel: window sum of the daily hazard", {
  f <- fn$compute_daily_foi(AR_total = 0.2, L = 364)
  # whole outbreak covered -> total attack rate
  expect_equal(fn$compute_ar_travel(f, entry_day = 1, D = 364), 0.2)
  expect_equal(fn$compute_ar_travel(f, entry_day = 50, D = 1000), 0.2)
  # explicit window
  expect_equal(fn$compute_ar_travel(f, entry_day = 100, D = 30), 1 - exp(-sum(f[100:129])))
  # entry day is clamped so the trip stays inside the outbreak
  expect_equal(fn$compute_ar_travel(f, entry_day = 360, D = 30), 1 - exp(-sum(f[335:364])))
  # longer trips (same entry) never lower the risk
  r <- vapply(c(7, 14, 30, 90), function(D) fn$compute_ar_travel(f, 120, D), 0)
  expect_true(all(diff(r) >= 0))
  # trips over the epidemic peak carry more risk than trips at the tail
  expect_gt(fn$compute_ar_travel(f, 160, 30), fn$compute_ar_travel(f, 5, 30))
})

test_that("compute_outcome: benefits, risks and BRR per 10,000", {
  o <- fn$compute_outcome(AR = 0.01, p_hosp = 0.02, p_death = 0.001,
                          p_sae_vacc = 1e-4, p_death_vacc = 2e-5, VE_hosp = 0.9, VE_death = 0.9)
  expect_equal(o$risk_nv_hosp, 1e4 * 0.01 * 0.02)
  expect_equal(o$averted_10k_sae, 1e4 * 0.01 * (0.02 + 0.001) * 0.9)
  expect_equal(o$averted_10k_death, 1e4 * 0.01 * 0.001 * 0.9)
  expect_equal(o$excess_10k_sae, 1e4 * 1e-4)
  expect_equal(o$excess_10k_death, 1e4 * 2e-5)
  expect_equal(o$brr_sae, o$averted_10k_sae / o$excess_10k_sae)
  expect_equal(o$brr_death, o$averted_10k_death / o$excess_10k_death)
})

test_that("compute_outcome: zero vaccine risk gives an undefined (NA) BRR", {
  o <- fn$compute_outcome(AR = 0.01, p_hosp = 0.02, p_death = 0.001,
                          p_sae_vacc = 1e-4, p_death_vacc = 0, VE_hosp = 0.9, VE_death = 0.9)
  expect_true(is.na(o$brr_death))      # 18-64: vaccine-attributable death risk fixed at 0
  expect_false(is.na(o$brr_sae))
})

test_that("compute_outcome: no infection, no benefit", {
  o <- fn$compute_outcome(AR = 0, p_hosp = 0.02, p_death = 0.001,
                          p_sae_vacc = 1e-4, p_death_vacc = 2e-5, VE_hosp = 0.9, VE_death = 0.9)
  expect_equal(o$averted_10k_sae, 0)
  expect_equal(o$brr_sae, 0)
})

test_that("compute_outcome with vector AR: BRR length behaviour is documented", {
  # Characterisation test. compute_outcome() uses ifelse() with a scalar test
  # (excess_10k_* does not depend on AR), so brr_* has length 1 when AR is a
  # vector. 06_brazil_travel_final_finite.R recomputes brr_sae/brr_death
  # elementwise for exactly this reason; this test guards that assumption.
  o <- fn$compute_outcome(AR = c(0.01, 0.02), p_hosp = 0.02, p_death = 0.001,
                          p_sae_vacc = 1e-4, p_death_vacc = 2e-5, VE_hosp = 0.9, VE_death = 0.9)
  expect_length(o$averted_10k_sae, 2)
  expect_length(o$brr_sae, 1)
  expect_equal(o$averted_10k_sae / o$excess_10k_sae, c(0.01, 0.02) * 1e4 * 0.021 * 0.9 / 1)
})
