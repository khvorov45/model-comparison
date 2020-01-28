# Fit the orinary logistic model
# Arseniy Khvorov
# Created 2020-01-16
# Last edit 2020-01-28

library(tidyverse)

# Directories to be used later
fit_logistic_dir <- "fit-logistic"
data_dir <- "data"

# Functions ===================================================================

read_hanam <- function(nme) {
  read_csv(
    file.path(data_dir, paste0(nme, ".csv")), 
    col_types = cols_only(
      status_bin = col_integer(),
      loghilb = col_double(),
      loghiub = col_double(),
      virus = col_character(),
      population = col_character()
    )
  ) %>%
    rename(status = status_bin)
}

predict_lr_one <- function(fit_one, logHIs) {
  predict(fit_one, data.frame(logHImid = logHIs), se.fit = TRUE) %>% 
    as_tibble() %>%
    mutate(
      fit_high = invlogit(fit + qnorm(0.975) * se.fit),
      fit_low = invlogit(fit - qnorm(0.975) * se.fit),
      fit = invlogit(fit),
      logHImid = logHIs,
      population = attr(fit_one, "population"),
      virus = attr(fit_one, "virus")
    )
}

fit_lr_one <- function(data, formula) {
  fit <- glm(formula, data, family = binomial(link = "logit"))
  attr(fit, "virus") <- unique(data$virus)
  attr(fit, "population") <- unique(data$population)
  fit
}

# Script ======================================================================

# Hanam data
han <- read_res(file.path(data_dir, "hanam-HI-exp.csv")) %>%
  bind_rows(read_res(file.path(data_dir, "hanam-HI-gen.csv")))

# Log HI's for which to calculate infection/protection probabilities
logHIs <- seq(0, 7.5, length.out = 101)

# Regular fit
fits <- han %>%
  group_split(virus, population) %>%
  map(fit_lr_one, status_bin ~ logHImid)

# Predictions from regular fits
preds <- fits %>% map_dfr(predict_lr_one, logHIs)
