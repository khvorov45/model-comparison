# Simulations to see if Cox produces reliable results with time at risk
# being proportional to follow-up time
#
# Arseniy Khvorov
# Created 2019/09/20
# Last edit 2019/12/02

library(tidyverse)
library(broom)
library(survival)
library(furrr)
library(extraDistr)

plan(multiprocess)

# Directories to be used later
cox_tarprop_dir <- "cox-tarprop"

# Settings ====================================================================

nsim <- 5

vary_list <- list(
  risk_prop_expected = c(0.01, seq(0.1, 0.9, 0.1), 1),
  long_prop_expected = seq(0, 1, 0.1)
)

# Functions ===================================================================

# Survival data with exponential survival time (hazard constant with time)
sim_cox <- function(nsam, beta_0, beta_logTitre, risk_prop_expected,
                    long_prop_expected, max_follow, kappa,
                    seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)

  pop <- tibble(
    logTitre = rnorm(nsam, 2, 2),
    loghazard = beta_0 + beta_logTitre * logTitre,
    time_at_risk_required = rexp(nsam, rate = exp(loghazard)),
    risk_proportion = case_when(
      risk_prop_expected == 1 ~ rep(1, nsam),
      risk_prop_expected == 0 ~ rep(0, nsam),
      TRUE ~ rbeta(
        nsam, risk_prop_expected * kappa, (1 - risk_prop_expected) * kappa
      )
    ),
    time_at_risk_max = max_follow * risk_proportion,
    status = as.integer(time_at_risk_required <= time_at_risk_max),
    is_long = rbern(nsam, long_prop_expected),
    time_recorded = if_else(
      status == 1,
      time_at_risk_required / risk_proportion,
      max_follow
    )
  ) %>% mutate(
    time_recorded = if_else(
      is_long == 1,
      time_recorded + runif(nsam, 0, 200),
      time_recorded
    )
  )
  
  attr(pop, "true_values") <- tibble(
    nsam, beta_0, beta_logTitre, risk_prop_expected, 
    long_prop_expected, max_follow, kappa
  ) 
  attr(pop, "seed") <- seed
  
  pop
}

# Fit the cox model
fit_cox <- function(data, formula) {
  coxph(formula, data) %>% 
    tidy() %>%
    select(term, estimate, std_error = std.error) %>%
    mutate(
      term = paste0("beta_", term),
      seed = attr(data, "seed")
    ) %>%
    inner_join(
      attr(data, "true_values") %>%
        pivot_longer(everything(), names_to = "term", values_to = "true_val"), 
      by = "term"
    ) %>%
    bind_cols(attr(data, "true_values") %>% slice(rep(1, n())))
}

# Simulate one dataset (from a dictionary)
sim_one_dd <- function(data_name, data_dict,
                       seed = sample.int(.Machine$integer.max, 1)) {
  sim_fun <- data_dict[[data_name]]$sim_fun
  sim_args <- data_dict[[data_name]]$sim_args
  data <- do.call(sim_fun, c(sim_args, seed = seed))
  attr(data, "name") <- data_name
  data
}

# Fit a model (from a dictionary)
fit_one_dd <- function(model_name, model_dict, data) {
  fit_fun <- model_dict[[model_name]]$fit_fun
  fit_args <- model_dict[[model_name]]$fit_args
  do.call(fit_fun, c(fit_args, list(data = data))) %>%
    mutate(model_name = model_name, data_name = attr(data, "name"))
}

# Simulate one dataset and fit one model from dictionaries
sim_fit_one <- function(data_name, data_dict, model_name, model_dict,
                        seed = sample.int(.Machine$integer.max, 1)) {
  fit_one_dd(model_name, model_dict, sim_one_dd(data_name, data_dict, seed))
}

# Simulate many datasets and fit models to them
sim_fit_many <- function(nsim, data_name, data_dict, model_name, model_dict,
                         init_seed = sample.int(.Machine$integer.max, 1)) {
  future_map_dfr(
    1:nsim, 
    function(i) sim_fit_one(
      data_name, data_dict, model_name, model_dict, init_seed + (i - 1)
    ),
    .options = future_options(
      packages = c("extraDistr", "broom", "survival", "tidyr")
    )
  )
}

vary_par <- function(par_vals, par_name, nsim, data_name, data_dict, 
                     model_name, model_dict,
                     init_seed = sample.int(.Machine$integer.max, 1)) {
  imap_dfr(
    par_vals,
    function(val, i) {
      data_dict[[data_name]]$sim_args[[par_name]] <- val
      sim_fit_many(
        nsim, data_name, data_dict, model_name, model_dict,
        init_seed + (i - 1) * nsim
      )
    }
  ) %>%
    mutate(par_varied = par_name)
}

vary_par_1aat <- function(vary_list, nsim, data_name, data_dict, 
                          model_name, model_dict,
                          init_seed = sample.int(.Machine$integer.max, 1)) {
  imap_dfr(
    vary_list,
    function(par_vals, par_name) {
      i <- which(names(vary_list) == par_name)
      vary_par(
        par_vals, par_name, 
        nsim, data_name, data_dict, 
        model_name, model_dict,
        init_seed + (i - 1) * nsim * length(vary_list[[par_name]])
      )
      
    }
    
  )
}

std_data <- function(kappa) {
  list(
    sim_fun = sim_cox,
    sim_args = list(
      nsam = 1e4, beta_0 = -3, beta_logTitre = -1.5, risk_prop_expected = 1, 
      long_prop_expected = 0, max_follow = 100, kappa = kappa
    )
  )
}

sim_all_data <- function(vary_list, nsim, data_dict, 
                         model_name, model_dict) {
  total_sims <- sum(map_dbl(vary_list, length)) * nsim
  imap_dfr(
    names(data_dict), 
    ~ vary_par_1aat(
      vary_list,
      nsim, .x, data_dict, model_name, model_dict, (.y - 1) * total_sims
    )
  )
}

# Script ======================================================================

data_dict <- list(
  std_1000 = std_data(1000),
  std_100 = std_data(100),
  std_10 = std_data(10),
  std_1 = std_data(1),
  std_05 = std_data(0.5)
)

model_dict <- list(
  std = list(
    fit_fun = fit_cox,
    fit_args = list(formula = Surv(time_recorded, status == 1) ~ logTitre)
  )
)

fits_summs <- sim_all_data(vary_list, nsim, data_dict, "std", model_dict)

write_csv(
  fits_summs,
  file.path(
    cox_tarprop_dir, 
    paste0("result-", nsim, "sims.csv")
  )
)
