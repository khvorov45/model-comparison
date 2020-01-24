# Fit the logistic model to data
# Arseniy Khvorov
# Created 2020/01/16
# Last edit 2020/01/24

library(tidyverse)
library(furrr)
library(rsample)
library(broom)

plan(multiprocess)

# Directories to be used later
fit_logistic_boot_dir <- "fit-logistic-boot"
data_dir <- "data"

# Functions ===================================================================

read_one <- function(name) {
  read_csv(file.path(data_dir, paste0(name, ".csv")))
}

read_kiddyvax <- function(name) {
  read_csv(
    file.path(data_dir, paste0(name, ".csv")),
    col_types = cols(
      id = col_integer(),
      virus = col_character(),
      status = col_integer(),
      infection_date = col_date("%Y-%m-%d"),
      start_date = col_date("%Y-%m-%d"),
      end_date = col_date("%Y-%m-%d"),
      hi = col_integer()
    )
  )
}

boot_fit <- function(dat, formula, n_res = 1000) {
  resamples <- bootstraps(dat, times = n_res)
  fit_one <- function(split, ind) {
    dat <- analysis(split)
    fit <- tryCatch(
      tidy(glm(formula, dat, family = binomial(link = "logit"))),
      error = function(e) return(NULL),
      warning = function(w) return(NULL)
    )
    if (!is.null(fit)) fit$ind <- ind
    fit
  }
  future_imap_dfr(
    resamples$splits, fit_one, .options = future_options(packages = "broom")
  )
}

save_lr_boot <- function(samples, name) {
  write_csv(
    samples,
    file.path(fit_logistic_boot_dir, paste0(name, ".csv"))
  )
}

# Script ======================================================================

# Hanam data
han <- read_one("hanam-HI-exp") %>%
  bind_rows(read_one("hanam-HI-gen")) %>%
  filter(virus != "H1N1seas")

# Kiddyvax data
kv <- read_kiddyvax("kiddyvax-main") %>%
  filter(virus %in% c("h1pdm", "bvic"))

# Fit hanam to viruses and populations separately
fit_bootstraps_han <- han %>%
  group_by(virus, population) %>%
  group_modify(~ boot_fit(.x, status_bin ~ logHImid, 1e3))
save_lr_boot(fit_bootstraps_han, "samples")

# Fit kiddyvax to viruses separately
fit_lr_boot_kv <- kv %>%
  mutate(logHImid = if_else(hi == 5, log(5), log(hi) + log(2) / 2)) %>%
  group_by(virus) %>%
  group_modify(~ boot_fit(.x, status ~ logHImid, 1e3))
save_lr_boot(fit_lr_boot_kv, "kiddyvax-samples")