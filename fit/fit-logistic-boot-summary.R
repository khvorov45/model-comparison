# Summary of bootstraps of logistic fit

library(tidyverse)

# Directories to be used later
fit_dir <- here::here("fit")
fit_lr_boot_out_dir <- file.path(fit_dir, "out-logistic-boot")
fit_sclr_boot_dir <- file.path(fit_dir, "out-sclr-boot")

# Functions ===================================================================

# Reads a set of results from path
read_res <- function(path, model = "lr") {
  res <- read_csv(path, col_types = cols())
  if (!"population" %in% names(res)) res$population <- "General"
  res %>%
    mutate(
      population = factor(
        population,
        levels = c("General", "Exposed")
      ),
      filename = paste0(
        str_replace(basename(path), ".csv", ""),
        paste0("-preds-", model, "-boot")
      )
    )
}

read_res_all <- function(dirpath, model = "lr") {
  map_dfr(tools::list_files_with_exts(dirpath, "csv"), read_res, model)
}

widen_bootstraps <- function(bootstr) {
  bootstr %>%
    pivot_wider(ind, names_from = "term", values_from = "estimate")
}

calc_probs <- function(out) {
  out %>%
    mutate(
      loghi = loghi,
      inf = 1 - 1 / (1 + exp(b0 + b1 * loghi)),
      prot = 1 - inf,
      prot_rel = 1 - inf / (1 - 1 / (1 + exp(b0 + b1 * log(5))))
    ) %>%
    pivot_longer(
      c(prot, inf, prot_rel),
      names_to = "prob_type", values_to = "prob"
    )
}

calc_probs_sclr <- function(out) {
  invlogit <- function(x) 1 - 1 / (1 + exp(x))
  out %>%
    mutate(
      loghi = loghi,
      prot = invlogit(beta_0 + beta_loghimid * loghi),
      inf = invlogit(theta) * (1 - prot),
    ) %>%
    pivot_longer(
      c(prot, inf),
      names_to = "prob_type", values_to = "prob"
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

full_summ <- function(boot_data, prob_fun = calc_probs) {
  boot_data %>%
    group_by(virus, population, filename) %>%
    group_modify(~ widen_bootstraps(.x)) %>%
    slice(rep(1:n(), each = length(loghis))) %>%
    mutate(loghi = rep_len(loghis, length.out = n())) %>%
    group_modify(~ prob_fun(.x)) %>%
    summ_probs() %>%
    split_summ()
}

save_summ <- function(dat, name) {
  write_csv(dat, file.path(fit_dir, paste0(name, ".csv")))
}

# Script ======================================================================

# Bootstrap samples
lr_boot <- read_res_all(fit_lr_boot_out_dir)
sclr_boot <- read_res_all(fit_sclr_boot_dir, "sclr")

# Log HI's for which to calculate infection/protection probabilities
loghis <- seq(0, 8.7, length.out = 101)

# Calculate predicated probabilities
preds <- full_summ(lr_boot)
preds_sclr <- full_summ(sclr_boot, calc_probs_sclr)

iwalk(preds, save_summ)
iwalk(preds_sclr, save_summ)
