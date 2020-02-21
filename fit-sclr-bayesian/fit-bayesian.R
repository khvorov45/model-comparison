# Fitting Hanam data in a Baysian way
# Arseniy Khvorov
# Created 2019-09-18
# Last edit 2019-01-28

library(tidyverse)
library(dclone)
library(rjags)

# Directories to be used later
fit_sclr_ba_dir <- "fit-sclr-bayesian"
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

read_kiddyvax <- function(nme) {
  read_csv(
    file.path(data_dir, paste0(nme, ".csv")), 
    col_types = cols_only(
      status = col_integer(),
      loghilb = col_double(),
      loghiub = col_double(),
      virus = col_character()
    )
  )
}

# Convert one chain's output to table (JAGS only)
tidy_mcmc_chain <- function(mcmcchain, nchain) {
  as_tibble(as.data.frame(mcmcchain)) %>%
    mutate(nchain = nchain, niter = row_number())
}

# Convert multiple chains output to table (JAGS only)
tidy_mcmc <- function(mcmcout) {
  imap_dfr(mcmcout, tidy_mcmc_chain)
}

# Fit JAGS model
fit_bayesian <- function(model, n_adapt, n_iter) {
  nchains <- length(model$inits)
  cl <- parallel::makeCluster(min(nchains, parallel::detectCores()))
  on.exit(parallel::stopCluster(cl))
  dclone::jags.parfit(
    cl = cl,
    data = c(model$data, n = nrow(model$data)),
    params = model$pars,
    model = model$filepath,
    inits = model$inits,
    n.chains = nchains,
    n.adapt = n_adapt,
    n.update = 0, # Will cut off manually
    n.iter = n_iter
  )
}

# Fit many models
fit_bayesian_many <- function(models, n_adapt, n_iter) {
  cl <- parallel::makeCluster(min(length(models), parallel::detectCores()))
  on.exit(parallel::stopCluster(cl))
  parLapply(
    cl, models, fit_bayesian, n_adapt = n_adapt, n_iter = n_iter
  ) %>% map(tidy_mcmc)
}

# Saves table output
save_model_output <- function(tidyout, model_name) {
  write_csv(tidyout, file.path(fit_sclr_ba_dir, paste0(model_name, ".csv")))
}

# Create a list with model data
create_jags_model <- function(dat, ...) {
  dat <- filter(dat, ...)
  list(
    filepath = file.path(fit_sclr_ba_dir, "sclr.jags"),
    data = dat,
    pars = c("lambda", "beta_0", "beta_hi", "loghimean", "loghisd"),
    inits = list(
      list(lambda = 0.1, beta_0 = 0, beta_HI = 0, logHImean = 0, logHIsd = 1),
      list(lambda = 0.5, beta_0 = 5, beta_HI = 2.5, logHImean = 2, logHIsd = 2),
      list(lambda = 0.9, beta_0 = -5, beta_HI = 5, logHImean = 4, logHIsd = 4)
    )
  )
}

# Script ======================================================================

# HI subset of Hanam data
han_hi <- map_dfr(c("hanam-hi-gen", "hanam-hi-exp"), read_hanam)

# Kiddyvax main study
kiddyvax_main <- read_kiddyvax("kiddyvaxmain")

# Model definition
models_hanam <- list(
  `hanam-hi-gen-h1pdm` = create_jags_model(
    han_hi, virus == "H1N1pdm", population == "General"
  ),
  `hanam-hi-gen-h3` = create_jags_model(
    han_hi, virus == "H3N2", population == "General"
  ),
  `hanam-hi-exp-h1pdm` = create_jags_model(
    han_hi, virus == "H1N1pdm", population == "Exposed"
  ),
  `hanam-hi-exp-h3` = create_jags_model(
    han_hi, virus == "H3N2", population == "Exposed"
  ),
  `kiddyvaxmain-h1pdm` = create_jags_model(kiddyvax_main, virus == "h1pdm"),
  `kiddyvaxmain-bvic` = create_jags_model(kiddyvax_main, virus == "bvic")
)

# Fit models
out <- fit_bayesian_many(models_hanam, n_adapt = 1e2, n_iter = 1e2)

# Save output
iwalk(out, save_model_output)
