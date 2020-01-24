# Fitting Hanam data in a Baysian way
# Arseniy Khvorov
# Created 2019/09/18
# Last edit 2019/12/11

library(tidyverse)
library(dclone)
library(rjags)

# Directories to be used later
fit_bayesian_dir <- "fit-bayesian"
data_dir <- "data"

# Functions ===================================================================

# Reads one csv file
read_one <- function(nme) {
  read_csv(file.path(data_dir, paste0(nme, ".csv")), col_types = cols())
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
  write_csv(tidyout, file.path(fit_bayesian_dir, paste0(model_name, ".csv")))
}

# Create a list with model data
han_model <- function(dat) {
  dat <- dat %>%
    select(status = status_bin, logHIlb, logHIub)
  list(
    filepath = file.path(fit_bayesian_dir, "hanam.jags"),
    data = dat,
    pars = c("lambda", "beta_0", "beta_HI", "logHImean", "logHIsd"),
    inits = list(
      list(lambda = 0.1, beta_0 = 0, beta_HI = 0, logHImean = 0, logHIsd = 1),
      list(lambda = 0.5, beta_0 = 5, beta_HI = 2.5, logHImean = 2, logHIsd = 2),
      list(lambda = 0.9, beta_0 = -5, beta_HI = 5, logHImean = 4, logHIsd = 4)
    )
  )
}
kv_model <- function(dat) {
  dat <- dat %>%
    select(status, logHIlb = loghilb, logHIub = loghiub)
  list(
    filepath = file.path(fit_bayesian_dir, "hanam.jags"),
    data = dat,
    pars = c("lambda", "beta_0", "beta_HI", "logHImean", "logHIsd"),
    inits = list(
      list(lambda = 0.1, beta_0 = 0, beta_HI = 0, logHImean = 0, logHIsd = 1),
      list(lambda = 0.5, beta_0 = 5, beta_HI = 2.5, logHImean = 2, logHIsd = 2),
      list(lambda = 0.9, beta_0 = -5, beta_HI = 5, logHImean = 4, logHIsd = 4)
    )
  )
}

# Script ======================================================================

# HI subset
han_hi_gen <- read_one("hanam-HI-gen")

# Exposed HI
han_hi_exp <- read_one("hanam-HI-exp")

# Kiddyvax main study
kiddyvax <- read_one("kiddyvax-main")

# Model definition
models_hanam <- list(
  #H1N1pdm = han_model(filter(han_hi_gen, virus == "H1N1pdm")),
  #H3N2 = han_model(filter(han_hi_gen, virus == "H3N2")),
  #H1N1pdmExp = han_model(filter(han_hi_exp, virus == "H1N1pdm")),
  #H3N2Exp = han_model(filter(han_hi_exp, virus == "H3N2")),
  `kiddyvax-H1N1pdm` = kv_model(filter(kiddyvax, virus == "h1pdm")),
  `kiddyvax-BVic` = kv_model(filter(kiddyvax, virus == "bvic"))
)

# Fit models
out <- fit_bayesian_many(models_hanam, n_adapt = 1e2, n_iter = 1e2)

# Save output
iwalk(out, save_model_output)
