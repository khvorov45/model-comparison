# Fit the scaled logit model to hanam data
# Arseniy Khvorov
# Created 2019/12/04
# Last edit 2019/12/04

library(tidyverse)
library(sclr)
library(furrr)
library(rsample)

plan(multiprocess)

# Directories to be used later
fit_sclr_boot_dir <- "fit-sclr-boot"
data_dir <- "data"

# Functions ===================================================================

boot_fit <- function(dat, formula, n_res = 1000) {
  resamples <- bootstraps(dat, times = n_res)
  fit_one <- function(split, ind) {
    dat <- analysis(split)
    fit <- tryCatch(
      tidy(sclr(formula, dat, algorithm = "newton-raphson")),
      error = function(e) return(NULL),
      warning = function(w) return(NULL)
    )
    if (!is.null(fit)) fit$ind <- ind
    fit
  }
  future_imap_dfr(resamples$splits, fit_one)
}

# Script ======================================================================

# Hanam data
han <- read_csv(file.path(data_dir, "hanam.csv"))

# Prepare
han_prep <- han %>%
  filter(!is.na(status) | !is.na(preHI), virus != "H1N1seas") %>%
  mutate(
    logHI = log(preHI),
    logHImid = case_when(
      preHI == 5 ~ log(5),
      preHI == 1280 ~ log(2560),
      TRUE ~ logHI + log(2) / 2
    ),
    status_bin = if_else(status == "Not infected", 0 ,1)
  ) %>%
  select(status, status_bin, preHI, logHI, logHImid, virus)

# Fit to viruses separately
fit_bootstraps <- han_prep %>% # To get se's for the infection curve
  group_by(virus) %>%
  group_modify(~ boot_fit(.x, status_bin ~ logHImid, 1e3))
write_csv(
  fit_bootstraps,
  file.path(fit_sclr_boot_dir, "samples.csv")
)
