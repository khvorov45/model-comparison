# Fit the scaled logit model

library(tidyverse)
library(sclr)

# Directories to be used later
fit_dir <- here::here("fit")
data_dir <- here::here("data")

# Functions ===================================================================

read_hanam <- function(nme) {
  read_csv(
    file.path(data_dir, paste0(nme, ".csv")),
    col_types = cols_only(
      status_bin = col_integer(),
      loghimid = col_double(),
      virus = col_character(),
      population = col_character()
    )
  ) %>%
    rename(status = status_bin)
}

read_kiddyvax <- function(nme) {
  read_csv(
    file.path(data_dir, paste0(nme, ".csv")),
    col_types = cols_only(
      virus = col_character(),
      status = col_integer(),
      loghimid = col_double()
    )
  )
}

predict_sclr_one <- function(fit_one, loghis) {
  predict(fit_one, data.frame(loghimid = loghis)) %>%
    select(prot_point, prot_l, prot_u) %>%
    mutate(
      loghimid = loghis,
      population = attr(fit_one, "population"),
      virus = attr(fit_one, "virus")
    )
}

fit_sclr_one <- function(data, formula) {
  fit <- sclr(formula, data)
  attr(fit, "virus") <- unique(data$virus)
  if ("population" %in% names(data)) {
    attr(fit, "population") <- unique(data$population)
  }
  fit
}

save_preds <- function(preds, nme) {
  write_csv(
    preds,
    file.path(fit_dir, paste0(nme, "-preds-sclr.csv"))
  )
}

# Script ======================================================================

# Log HI's for which to calculate infection/protection probabilities
loghis <- seq(0, 8.7, length.out = 101)

# Hanam data
han <- bind_rows(read_hanam("hanam-hi-exp"), read_hanam("hanam-hi-gen")) %>%
  filter(virus != "H1N1seas")

# Kiddyvax data
kv <- read_kiddyvax("kiddyvaxmain") %>%
  filter(virus %in% c("bvic", "h1pdm"))

# Regular fit
fits_han <- han %>%
  group_split(virus, population) %>%
  map(fit_sclr_one, status ~ loghimid)

fits_kv <- kv %>%
  group_split(virus) %>%
  map(fit_sclr_one, status ~ loghimid)

# Predictions from regular fits
all_preds <- map(
  list("hanam-hi" = fits_han, "kiddyvaxmain" = fits_kv),
  ~ map_dfr(.x, predict_sclr_one, loghis)
)

iwalk(all_preds, save_preds)
