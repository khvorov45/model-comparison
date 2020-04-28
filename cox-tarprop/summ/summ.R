# Summary of Cox simulations

library(tidyverse)

# Directories used
sim_dir <- here::here("sim")
summ_dir <- here::here("summ")

# Functions ===================================================================

dont_group_by <- function(.tbl, var_names) {
  group_by(
    .tbl,
    !!!syms(names(.tbl)[!names(.tbl) %in% var_names])
  )
}

# Script ======================================================================

res <- read_csv(file.path(sim_dir, "sim.csv"), col_types = cols())

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

write_csv(res_sum, file.path(summ_dir, "summ.csv"))
