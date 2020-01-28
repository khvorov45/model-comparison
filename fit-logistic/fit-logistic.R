# Fit the orinary logistic model
# Arseniy Khvorov
# Created 2020-01-16
# Last edit 2020-01-29

library(tidyverse)

# Directories to be used later
fit_logistic_dir <- "fit-logistic"
data_dir <- "data"

# Functions ===================================================================

invlogit <- function(x) 1 - 1 / (1 + exp(x))

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

predict_lr_one <- function(fit_one, loghis) {
  predict(fit_one, data.frame(loghimid = loghis), se.fit = TRUE) %>% 
    as_tibble() %>%
    select(fit, se = se.fit) %>%
    mutate(
      fit_high = invlogit(fit + qnorm(0.975) * se),
      fit_low = invlogit(fit - qnorm(0.975) * se),
      fit = invlogit(fit),
      loghimid = loghis,
      population = attr(fit_one, "population"),
      virus = attr(fit_one, "virus")
    )
}

fit_lr_one <- function(data, formula) {
  fit <- glm(formula, data, family = binomial(link = "logit"))
  attr(fit, "virus") <- unique(data$virus)
  if ("population" %in% names(data)) 
    attr(fit, "population") <- unique(data$population)
  fit
}

save_preds <- function(preds, nme) {
  write_csv(
    preds,
    file.path(fit_logistic_dir, paste0(nme, ".csv"))
  )
}

# Script ======================================================================

# Log HI's for which to calculate infection/protection probabilities
loghis <- seq(0, 8.7, length.out = 101)

# Hanam data
han <- map_dfr(c("hanam-HI-exp", "hanam-HI-gen"), read_hanam) %>%
  filter(virus != "H1N1seas")

# Kiddyvax data
kv <- read_kiddyvax("kiddyvaxmain") %>%
  filter(virus %in% c("bvic", "h1pdm"))

# Regular fit
fits_han <- han %>%
  group_split(virus, population) %>%
  map(fit_lr_one, status ~ loghimid)

fits_kv <- kv %>%
  group_split(virus) %>%
  map(fit_lr_one, status ~ loghimid)

# Predictions from regular fits
preds_han <- fits_han %>% map_dfr(predict_lr_one, loghis)
save_preds(preds_han, "hanam-hi")

preds_kv <- fits_kv %>% map_dfr(predict_lr_one, loghis)
save_preds(preds_kv, "kiddyvaxmain")
