# Summary of simualtions in logistic
# Arseniy Khvorov
# Created 2019/12/03
# Last edit 2019/12/03

library(tidyverse)

# Directories to be used later
logistic_dir <- "logistic"
logistic_summ_dir <- "logistic-summary"

# Script ======================================================================

res <- read_csv(file.path(logistic_dir, "result-10000sims.csv"))

summ <- res %>%
  mutate(
    estimate = if_else(model == "logistic", -estimate, estimate)
  ) %>%
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

write_csv(summ, file.path(logistic_summ_dir, "summary-10000sims.csv"))

