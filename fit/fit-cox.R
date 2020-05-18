# Fitting the cox model to the kiddyvax data

library(tidyverse)
library(survival)

# Directories used
fit_dir <- here::here("fit")
data_dir <- here::here("data")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

fit_cox_one <- function(dat, formula) {
  coxph(formula, dat, model = TRUE)
}

fit_cox_manymod <- function(dat, formulas) {
  map(formulas, ~ fit_cox_one(dat, .x))
}

predict_cox_one <- function(fit, newdata) {
  terms_noy <- delete.response(fit$terms)
  model_mat <- model.frame(terms_noy, newdata)
  ests_mat <- matrix(coef(fit), ncol = 1)
  loghr <- apply(model_mat, 1, function(x) x %*% ests_mat)
  sds <- map_dbl(1:nrow(model_mat), function(i) {
    x <- as.matrix(model_mat[i, ])
    variance <- x %*% vcov(fit) %*% t(x)
    sqrt(variance)
  })
  newdata %>%
    mutate(
      loghr = loghr,
      se = sds,
      loghr_low = loghr - qnorm(0.975) * se,
      loghr_high = loghr + qnorm(0.975) * se,
      prot = 1 - exp(loghr),
      prot_low = 1 - exp(loghr_high),
      prot_high = 1 - exp(loghr_low),
      se_wrong = sqrt(sum(vcov(fit))),
      loghr_low_wrong = loghr - qnorm(0.975) * se_wrong,
      loghr_high_wrong = loghr + qnorm(0.975) * se_wrong,
      prot_low_wrong = 1 - exp(loghr_high_wrong),
      prot_high_wrong = 1 - exp(loghr_low_wrong)
    )
}

save_cox_pred <- function(pred, name) {
  write_csv(pred, file.path(fit_dir, glue::glue("{name}-preds-cox.csv")))
}

# Script ======================================================================

# Kiddyvax data
kv <- read_data("kiddyvaxmain") %>%
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
loghis <- seq(0, 8.7, length.out = 1001)
loghis_rel <- loghis - log(5)
preds <- fits_cox %>%
  map_dfr(
    predict_cox_one,
    tibble(loghi = loghis, loghimid = loghis_rel),
    .id = "virus"
  )
save_cox_pred(preds, "kiddyvaxmain")

# Sophia's data
sophia <- read_data("sophia")

models <- list(
  sophia = formula(Surv(t, event == 1) ~ postvax + proxy),
  me = formula(Surv(start, end, event == 1) ~ postvax + proxy + cluster(hhid))
)

sophia_fits <- sophia %>%
  group_split(model) %>%
  map(~ fit_cox_manymod(.x, models))
names(sophia_fits) <- group_keys(sophia, model)$model

sophia_preds <- sophia_fits %>%
  map_dfr(function(fits, virus_name) {
    map_dfr(fits, function(fit, model_name) {
      predict_cox_one(
        fit, tibble(loghi = loghis, postvax = loghis_rel, proxy = 0)
      )
    }, .id = "model")
  }, .id = "virus")
save_cox_pred(sophia_preds, "sophia")
