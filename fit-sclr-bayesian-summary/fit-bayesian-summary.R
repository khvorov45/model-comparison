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

calc_probs <- function(loghi, out) {
  out %>%
    mutate(
      loghi = loghi,
      prot = 1 - 1 / (1 + exp(beta_0 + beta_HI * loghi)),
      inf = lambda * (1 - prot),
    ) %>%
    pivot_longer(c(prot, inf), names_to = "prob_type", values_to = "prob")
}

calc_all_probs <- function(out, loghis) {
  map_dfr(loghis, calc_probs, out)
}

sum_outprob <- function(outprob) {
  outprob %>%
    group_by(loghi, prob_type) %>%
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
loghis <- seq(0, 8.7, length.out = 101)
outprobs <- map(out, calc_all_probs, loghis)
outprobs_sum <- map(outprobs, sum_outprob)

# Prior distributions
han_priors_r <- tibble(
  lambda = runif(5e4, 0, 1),
  beta_0 = rnorm(5e4, -15, 10),
  beta_HI = rnorm(5e4, 5, 5)
) %>% calc_all_probs(loghis) %>% sum_outprob()
names(han_priors_r)[3:5] <- paste0(names(han_priors_r)[3:5], "_prior")

# Add prior distributions
outprobs_sum_prior <- map(
  outprobs_sum, ~ inner_join(.x, han_priors_r, by = c("loghi", "prob_type"))
)

iwalk(outprobs_sum_prior, save_summ)
