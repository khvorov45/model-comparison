# Summary of Cox simulations
# Arseniy Khvorov
# Created 2019/12/02
# Last edit 2019/12/02

library(tidyverse)

# Directories to be used later
cox_tarprop_dir <- "cox-tarprop"
cox_tarprop_summ_dir <- "cox-tarprop-summary"

# Functions ===================================================================

dont_group_by <- function(.tbl, var_names) {
  group_by(
    .tbl, 
    !!!syms(names(.tbl)[!names(.tbl) %in% var_names])
  )
}

# Script ======================================================================

res <- read_csv(file.path(cox_tarprop_dir, "result-10000sims.csv"))

res_sum <- res %>%
  mutate(
    conf_low = estimate - qnorm(0.975) * std_error,
    conf_high = estimate + qnorm(0.975) * std_error
  ) %>%
  dont_group_by(c("estimate", "std_error", "seed", "conf_low", "conf_high")) %>%
  summarise(
    est_mean = mean(estimate), 
    est_sd = sd(estimate), 
    se_mean = mean(std_error),
    coverage = sum((true_val >= conf_low) & (true_val <= conf_high)) / n(),
    nsim = length(unique(seed))
  ) %>%
  ungroup()

write_csv(
  res_sum,
  file.path(
    cox_tarprop_summ_dir, 
    paste0("summary-10000sims.csv")
  )
)

