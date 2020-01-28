# Summarising MCMC output
# Arseniy Khvorov
# Created 2019/09/24
# Last edit 2019/11/26

library(tidyverse)

# Directories to be used later
fit_sclr_ba_dir <- "fit-sclr-bayesian"
fit_sclr_ba_summ_dir <- "fit-sclr-bayesian-summary"

# Functions ===================================================================

# Read one set of results
read_one <- function(filepath) {
  read_csv(filepath, col_types = cols())
}

calc_probs <- function(logHI, out) {
  out %>%
    mutate(
      logHI = logHI,
      prot = 1 - 1 / (1 + exp(beta_0 + beta_HI * logHI)),
      inf = lambda * (1 - prot),
    ) %>%
    pivot_longer(c(prot, inf), names_to = "prob_type", values_to = "prob")
}

calc_all_probs <- function(out, logHIs) {
  map_dfr(logHIs, calc_probs, out)
}

sum_outprob <- function(outprob) {
  outprob %>%
    group_by(logHI, prob_type) %>%
    summarise(
      prob_lb = quantile(prob, 0.025),
      prob_med = quantile(prob, 0.5),
      prob_ub = quantile(prob, 0.975)
    ) %>%
    ungroup()
}

# Save the summary file
save_summ <- function(summ, name) {
  write_csv(
    summ, 
    file.path(fit_sclr_ba_summ_dir, paste0(name, ".csv"))
  )
}

# Script ======================================================================

# Model output
out_files <- tools::list_files_with_exts(fit_sclr_ba_dir, "csv")
out <- map(out_files, read_one) %>% 
  map(filter, niter > 6e4) %>% 
  map(sample_n, 5e4)
names(out) <- str_replace(basename(out_files), ".csv", "")

# Summary
logHIs <- seq(0, 7.5, length.out = 101)
outprobs <- map(out, calc_all_probs, logHIs)
outprobs_sum <- map(outprobs, sum_outprob)

# Prior distributions
han_priors_r <- tibble(
  lambda = runif(5e4, 0, 1),
  beta_0 = rnorm(5e4, -15, 10),
  beta_HI = rnorm(5e4, 5, 5)
) %>% calc_all_probs(logHIs) %>% sum_outprob()
names(han_priors_r)[3:5] <- paste0(names(han_priors_r)[3:5], "_prior")

# Add prior distributions
outprobs_sum_prior <- map(
  outprobs_sum, ~ inner_join(.x, han_priors_r, by = c("logHI", "prob_type"))
)

iwalk(outprobs_sum_prior, save_summ)
