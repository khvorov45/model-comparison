# Fit the scaled logit model to hanam data
# Arseniy Khvorov
# Created 2020/01/16
# Last edit 2020/01/16

library(tidyverse)
library(sclr)
library(furrr)
library(rsample)
library(broom)

plan(multiprocess)

# Directories to be used later
fit_logistic_boot_dir <- "fit-logistic-boot"
data_dir <- "data"

# Functions ===================================================================

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

# Script ======================================================================

# Hanam data
han <- read_csv(file.path(data_dir, "hanam-HI-exp.csv")) %>%
  bind_rows(read_csv(file.path(data_dir, "hanam-HI-gen.csv"))) %>%
  filter(virus != "H1N1seas")

# Fit to viruses and populations separately
fit_bootstraps <- han %>%
  group_by(virus, population) %>%
  group_modify(~ boot_fit(.x, status_bin ~ logHImid, 1e3))
write_csv(
  fit_bootstraps,
  file.path(fit_logistic_boot_dir, "samples.csv")
)
