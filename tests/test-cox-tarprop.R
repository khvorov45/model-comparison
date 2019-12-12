# Testing cox-tarprop
# Arseniy Khvorov
# Created 2019/11/27
# Last edit 2019/11/27

library(testthat)

test_that("cox simulation works", {
  coxdat <- sim_cox(
    nsam = 100, beta_0 = -5, beta_logTitre = 1.5, 
    risk_prop_expected = 0.7, long_prop_expected = 0.1, max_follow = 100
  )
  expect_equal(
    with(coxdat, time_recorded[!is_long & status]), 
    with(coxdat, (time_at_risk_required / risk_proportion)[!is_long & status])
  )
  expect_equal(
    with(coxdat, time_recorded[!coxdat$status & !is_long]), 
    rep(100, sum(!coxdat$status & !coxdat$is_long))
  )
})

test_that("fitting function works", {
  coxdat <- sim_cox(
    nsam = 1e5, beta_0 = -5, beta_logTitre = 1.5, 
    risk_prop_expected = 1, long_prop_expected = 0, max_follow = 100
  )
  fit <- fit_cox(coxdat, Surv(time_recorded, status == 1) ~ logTitre)
  expect_equal(fit$estimate, fit$true_val, tol = 0.02)
})

test_that("dictionaries work", {
  coxdat <- sim_one_dd("std", data_dict)
  expect_equal(nrow(coxdat), data_dict$std$sim_args$nsam)
  fit <- fit_one_dd("std", model_dict, coxdat)
  fit <- sim_fit_one("std", data_dict, "std", model_dict)
})

test_that("many simulations work", {
  sims <- sim_fit_many(10, "std", data_dict, "std", model_dict, 100)
  expect_equal(sims$seed, seq(100, 109))
})

test_that("parameter variation works", {
  sims <- vary_par(
    c(0, 0.5, 1), "long_prop_expected",
    10, "std", data_dict, "std", model_dict, 100
  )
  expect_equal(sims$seed, seq(100, 129))
})

test_that("multiple (one at a time) parameter variation works", {
  var_list <- list(
    long_prop_expected = c(0, 0.5, 1),
    risk_prop_expected = c(0.1, 0.4, 0.7)
  )
  sims <- vary_par_1aat(
    var_list,
    10, "std", data_dict, "std", model_dict, 100
  )
  expect_equal(sims$seed, seq(100, 159))
})


