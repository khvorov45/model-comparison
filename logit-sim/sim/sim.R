# Simulations to see how logistic regression behaves when fit to data with
# baseline below 1.

library(tidyverse)
library(sclr)
library(extraDistr)
library(broom)
library(future)
library(furrr)

plan(multiprocess)

# Directories to be used later
sim_dir <- here::here("sim")

# Settings ====================================================================

nsim <- 1e4 # Number of simulations

# Simulation parameters

simpars <- list(
  std = list(nsam = 1e3, lambda = 0.5, beta_0 = -5, beta_logTitre = 1.5),
  small = list(nsam = 100, lambda = 0.5, beta_0 = -5, beta_logTitre = 1.5)
)

vary_list <- list(
  nsam = c(50, seq(100, 1e3, 100), 1e4),
  lambda = seq(0.1, 1, 0.1)
)

# Functions ===================================================================

sim_pop <- function(nsam, lambda, beta_0, beta_logTitre,
                    seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  pop <- tibble(
    logTitre = rnorm(nsam, 2, 2),
    prob = lambda / (1 + exp(beta_0 + beta_logTitre * logTitre)),
    status = rbern(nsam, prob)
  )
  attr(pop, "seed") <- seed
  theta <- log(lambda / (1 - lambda))
  attr(pop, "true_vals") <- tibble(nsam, theta, beta_0, beta_logTitre)
  pop
}

attach_true <- function(.tbl, data) {
  true_vals <- pivot_longer(
    attr(data, "true_vals"), everything(),
    names_to = "term", values_to = "true_value"
  )
  inner_join(.tbl, true_vals, by = "term") %>%
    bind_cols(attr(data, "true_vals") %>% slice(rep(1, nrow(.tbl))))
}

fit_lr <- function(data) {
  lr_fit <- tryCatch(
    withCallingHandlers(
      {
        glm(status ~ logTitre, binomial(link = "logit"), data) %>%
          tidy(conf.int = TRUE)
      },
      warning = function(w) {
        if (conditionMessage(w) ==
          "glm.fit: fitted probabilities numerically 0 or 1 occurred") {
          invokeRestart("muffleWarning")
        }
      }
    ),
    warning = function(w) {
      return(NULL)
    },
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(lr_fit)) {
    return(NULL)
  }
  lr_fit %>%
    mutate(
      model = "logistic", seed = attr(data, "seed"),
      term = recode(
        term,
        "(Intercept)" = "beta_0", "logTitre" = "beta_logTitre"
      )
    ) %>%
    select(
      term, estimate,
      std_error = std.error, conf_low = conf.low,
      conf_high = conf.high, model, seed
    ) %>%
    attach_true(data)
}

fit_sclr <- function(data) {
  sclr_fit <- tryCatch(
    sclr(status ~ logTitre, data, algorithm = "newton-raphson"),
    warning = function(w) {
      return(NULL)
    },
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(sclr_fit)) {
    return(NULL)
  }
  sclr_fit %>%
    tidy() %>%
    mutate(model = "scaled_logit", seed = attr(data, "seed")) %>%
    attach_true(data)
}

fit_both <- function(data) {
  bind_rows(fit_lr(data), fit_sclr(data))
}

simfit <- function(data_name, data_dict,
                   seed = sample.int(.Machine$integer.max, 1)) {
  do.call(sim_pop, c(data_dict[[data_name]], seed = seed)) %>% fit_both()
}

simfit_many <- function(nsim, data_name, data_dict,
                        init_seed = sample.int(.Machine$integer.max, 1)) {
  future_map_dfr(
    1:nsim,
    ~ simfit(data_name, data_dict, seed = init_seed + (.x - 1))
  )
}

vary_par <- function(par_vals, par_name, nsim, data_name, data_dict,
                     init_seed = sample.int(.Machine$integer.max, 1)) {
  imap_dfr(
    par_vals,
    function(val, ind) {
      data_dict[[data_name]][[par_name]] <- val
      simfit_many(nsim, data_name, data_dict, init_seed + (ind - 1) * nsim)
    }
  ) %>% mutate(par_varied = par_name)
}

vary_pars_1aat <- function(vary_list, nsim, data_name, data_dict,
                           init_seed = sample.int(.Machine$integer.max, 1)) {
  imap_dfr(
    vary_list,
    function(vals, name) {
      ind <- which(names(vary_list) == name)
      vary_par(
        vals, name, nsim, data_name, data_dict,
        init_seed = init_seed + (ind - 1) * nsim * length(vals)
      )
    }
  )
}

# Script ======================================================================

sims <- vary_par(vary_list$nsam, "nsam", nsim, "std", simpars, 20191203) %>%
  bind_rows(vary_par(
    vary_list$lambda, "lambda", nsim, "small", simpars,
    20191203 + length(vary_list$nsam) * nsim
  ))

write_csv(sims, file.path(sim_dir, "sim.csv"))
