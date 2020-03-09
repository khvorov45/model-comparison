# Summary of simualtions in logistic
# Arseniy Khvorov
# Created 2019/12/03
# Last edit 2019/12/03

library(tidyverse)

# Directories to be used later
logistic_dir <- "logistic"
logistic_summ_dir <- "logistic-summary"

# Functions ===================================================================

negate_lr_est <- function(dat) {
  mutate(dat, estimate = if_else(model == "logistic", -estimate, estimate))
}

sclr_curve <- function(x, l, b0, b1) l / (1 + exp(b0 + b1 * x))

calc_pred <- function(x, dat) {
  mutate(
    dat,
    logTitre = x,
    fit = sclr_curve(logTitre, lambda, beta_0, beta_logTitre)
  )
}

save_csv <- function(dat, name) {
  write_csv(dat, file.path(logistic_summ_dir, glue::glue("{name}.csv")))
}

# Script ======================================================================

res <- read_csv(file.path(logistic_dir, "result-10000sims.csv"))

summ <- res %>%
  negate_lr_est() %>%
  group_by(
    term, model, true_value, nsam, theta, beta_0, beta_logTitre, par_varied
  ) %>%
  summarise(
    est_mean = mean(estimate),
    est_sd = sd(estimate),
    se_mean = mean(std_error),
    coverage = sum(estimate <= conf_high & estimate >= conf_low) / n(),
    converged = length(unique(seed)) / 1e4
  )

save_csv(summ, "summary-10000sims")

# Predictions
preds <- res %>%
  negate_lr_est() %>%
  select(
    term, estimate, model, seed, nsam,
    theta_true = theta, beta0_true = beta_0,
    beta_logTitre_true = beta_logTitre, par_varied
  ) %>%
  pivot_wider(
    names_from = "term", values_from = "estimate"
  ) %>%
  mutate(
    lambda = if_else(model == "logistic", 1, 1 - 1 / (1 + exp(theta)))
  ) %>%
  map_dfr(
    seq(0, 7, length.out = 101), calc_pred, .
  ) %>%
  group_by(
    model, nsam, logTitre, theta_true, beta0_true, beta_logTitre_true,
    par_varied
  ) %>%
  summarise(fit = list(quantile(fit, c(0.025, 0.5, 0.975)))) %>%
  ungroup() %>%
  unnest_wider(fit)

save_csv(preds, "preds-10000sims")
