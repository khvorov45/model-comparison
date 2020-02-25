# Summary of bootstraps of sclr fit
# Arseniy Khvorov
# Created 2020-02-25
# Last edit 2020-02-25

library(tidyverse)

# Directories to be used later
fit_sclr_boot_dir <- "fit-sclr-boot"
fit_sclr_bs_dir <- "fit-sclr-boot-summary"

# Functions ===================================================================

# Reads a set of results from path
read_res <- function(path) {
  res <- read_csv(path, col_types = cols())
  if (!"population" %in% names(res)) res$population <- "General"
  res %>%
    mutate(
      population = factor(
        population,
        levels = c("General", "Exposed") 
      ),
      filename = basename(path) 
    )
}

read_res_all <- function(dirpath) {
  map_dfr(tools::list_files_with_exts(dirpath, "csv"), read_res)
}

widen_bootstraps <- function(bootstr) {
  bootstr %>%
    pivot_wider(ind, names_from = "term", values_from = "estimate")
}

calc_probs <- function(out) {
  out %>%
    mutate(
      loghi = loghi,
      inf = 1 - 1 / (1 + exp(beta_0 + beta_loghimid * loghi)),
      prot = 1 - inf,
      prot_rel = 1 - inf / (1 - 1 / (1 + exp(beta_0 + beta_loghimid * log(5))))
    ) %>%
    pivot_longer(
      c(prot, inf, prot_rel), names_to = "prob_type", values_to = "prob"
    )
}

summ_probs <- function(dat) {
  dat %>%
    group_by(virus, population, filename, prob_type, loghi) %>%
    summarise(
      prob_lb = quantile(prob, 0.025),
      prob_med = median(prob),
      prob_ub = quantile(prob, 0.975)
    ) %>%
    ungroup()
}

split_summ <- function(dat) {
  dat <- dat %>%
    group_by(filename)
  dat_keys <- group_keys(dat) %>% pull(filename)
  dat <- group_split(dat, keep = FALSE)
  names(dat) <- dat_keys
  dat
}

save_summ <- function(dat, name) {
  write_csv(dat, file.path(fit_sclr_bs_dir, name))
}

# Script ======================================================================

# Bootstrap samples
sclr_boot <- read_res_all(fit_sclr_boot_dir)

# Log HI's for which to calculate infection/protection probabilities
loghis <- seq(0, 8.7, length.out = 101)

# Calculate predicated probabilities
preds <- sclr_boot %>%
  group_by(virus, population, filename) %>%
  group_modify(~ widen_bootstraps(.x)) %>%
  slice(rep(1:n(), each = length(loghis))) %>%
  mutate(loghi = rep_len(loghis, length.out = n())) %>%
  group_modify(~ calc_probs(.x)) %>%
  summ_probs() %>%
  split_summ()

iwalk(preds, save_summ)
