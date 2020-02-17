# Fitting the cox model to the kiddyvax data
# Arseniy Khvorov
# Created 2020-02-17
# Last edit 2020-02-17

library(tidyverse)
library(survival)

# Directories used
fit_cox_dir <- "fit-cox"
data_dir <- "data"

# Functions ===================================================================

read_kiddyvax <- function(nme) {
  read_csv(
    file.path(data_dir, paste0(nme, ".csv")),
    col_types = cols_only(
      virus = col_character(),
      status = col_integer(),
      loghimid = col_double(),
      infection_date = col_date(),
      start_date = col_date(),
      end_date = col_date()
    )
  )
}

fit_cox_one <- function(dat, formula) {
  coxph(formula, dat, model = TRUE)
}

predict_cox_one <- function(fit, loghis, loghirel) {
  tibble(
    loghimid = loghis,
    loghr_rel = coef(fit) * (loghimid - loghirel),
    loghr_rel_se = sqrt((loghimid - loghirel)^2 * vcov(fit)[1, 1]),
    loghr_rel_low = loghr_rel - qnorm(0.975) * loghr_rel_se,
    loghr_rel_high = loghr_rel + qnorm(0.975) * loghr_rel_se,
    hr_rel = exp(loghr_rel),
    hr_rel_low = exp(loghr_rel_low),
    hr_rel_high = exp(loghr_rel_high),
    prot_rel = 1 - hr_rel,
    prot_rel_low = 1 - hr_rel_high,
    prot_rel_high = 1 - hr_rel_low
  )
}

save_cox_pred <- function(pred, name) {
  write_csv(pred, file.path(fit_cox_dir, glue::glue("{name}.csv")))
}

# Script ======================================================================

# Kiddyvax data
kv <- read_kiddyvax("kiddyvaxmain") %>%
  filter(virus %in% c("bvic", "h1pdm"))

# Add follow-up
kv_fup <- kv %>%
  filter(!is.na(status)) %>%
  mutate(
    end_of_risk_date = if_else(status == 1, infection_date, end_date),
    fup_days = end_of_risk_date - start_date
  )

# Fit the models
fits_cox <- kv_fup %>%
  group_split(virus) %>%
  map(
    ~ fit_cox_one(
      .x, Surv(fup_days, event = status == 1, type = "right") ~ loghimid
    )
  )
names(fits_cox) <- group_keys(kv_fup, virus)$virus

# Predict
loghis <- seq(0, 8.7, length.out = 101)
loghirel <- log(5)
preds <- fits_cox %>%
  map_dfr(predict_cox_one, loghis, loghirel, .id = "virus")
save_cox_pred(preds, "kiddyvaxmain")



temp <- predict(fits_cox[[1]], data.frame(loghimid = log(5)), type = "lp")
attr(temp, "constant")
