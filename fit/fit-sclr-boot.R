# Fit the scaled logit model, bootstraps

library(tidyverse)
library(sclr)
library(furrr)
library(rsample)

plan(multiprocess)

# Directories to be used later
fit_sclr_boot_dir <- here::here("fit", "out-sclr-boot")
data_dir <- here::here("data")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

boot_fit <- function(dat, formula, n_res = 1000) {
  resamples <- bootstraps(dat, times = n_res)
  fit_one <- function(split, ind) {
    dat <- analysis(split)
    fit <- tryCatch(
      tidy(sclr(formula, dat, algorithm = "newton-raphson")),
      error = function(e) {
        return(NULL)
      },
      warning = function(w) {
        return(NULL)
      }
    )
    if (!is.null(fit)) fit$ind <- ind
    fit
  }
  future_imap_dfr(resamples$splits, fit_one)
}

save_lr_boot <- function(samples, name) {
  if (!dir.exists(fit_sclr_boot_dir)) dir.create(fit_sclr_boot_dir)
  write_csv(
    samples,
    file.path(fit_sclr_boot_dir, paste0(name, ".csv"))
  )
}

# Script ======================================================================

# Hanam data
han <- map_dfr(c("hanam-hi-exp", "hanam-hi-gen"), read_data) %>%
  filter(virus != "H1N1seas")

# Kiddyvax data
kv <- read_data("kiddyvaxmain") %>%
  filter(virus %in% c("bvic", "h1pdm"))

# Fit hanam to viruses and populations separately
fit_bootstraps_han <- han %>%
  group_by(virus, population) %>%
  group_modify(~ boot_fit(.x, status_bin ~ loghimid, 1e4))
save_lr_boot(fit_bootstraps_han, "hanam-hi")

# Fit kiddyvax to viruses separately
fit_lr_boot_kv <- kv %>%
  group_by(virus) %>%
  group_modify(~ boot_fit(.x, status ~ loghimid, 1e4))
save_lr_boot(fit_lr_boot_kv, "kiddyvaxmain")
