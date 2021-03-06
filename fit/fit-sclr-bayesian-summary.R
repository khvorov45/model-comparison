# Summarising MCMC output

library(tidyverse)

# Directories to be used later
fit_dir <- here::here("fit")
sclr_ba_out_dir <- file.path(fit_dir, "out-sclr-bayesian")

# Functions ===================================================================

# Read one set of results
read_one <- function(filepath) {
  read_csv(filepath, col_types = cols())
}

calc_probs <- function(loghi, out) {
  out %>%
    mutate(
      loghi = loghi,
      prot = 1 - 1 / (1 + exp(beta_0 + beta_hi * loghi)),
      inf = lambda * (1 - prot),
    ) %>%
    pivot_longer(c(prot, inf), names_to = "prob_type", values_to = "prob")
}

calc_all_probs <- function(out, loghis) {
  map_dfr(loghis, calc_probs, out)
}

sum_outprob <- function(outprob) {
  outprob %>%
    group_by(loghi, prob_type) %>%
    summarise(
      prob_lb = quantile(prob, 0.025),
      prob_med = quantile(prob, 0.5),
      prob_ub = quantile(prob, 0.975)
    ) %>%
    ungroup()
}

find_prot_titre_levels <- function(out, prob_name, level = 0.5) {
  low <- -100
  high <- 100

  get_prot_prob <- function(var_name, titre, out) {
    calc_probs(titre, out) %>%
      sum_outprob() %>%
      filter(prob_type == "prot") %>%
      pull(var_name)
  }

  cnt <- 0
  med <- 0
  while (TRUE) {
    med_prev <- med
    med <- median(c(low, high))
    print(med)
    if (med == med_prev) cnt <- cnt + 1
    if (cnt == 5) {
      return(med)
    }
    pred_med <- get_prot_prob(prob_name, med, out)
    if (near(pred_med, level)) {
      return(med)
    }
    if (pred_med < level) {
      low <- med
    } else {
      high <- med
    }
  }
}

find_prot_titre_levels_all <- function(out, level = 0.5) {
  tibble(
    level = level,
    loghi_lb = find_prot_titre_levels(out, "prob_ub", level),
    loghi_med = find_prot_titre_levels(out, "prob_med", level),
    loghi_ub = find_prot_titre_levels(out, "prob_lb", level),
    hi_lb = exp(loghi_lb),
    hi_med = exp(loghi_med),
    hi_ub = exp(loghi_ub)
  )
}

addinfo <- function(summ, name) {
  summ %>%
    mutate(
      dataset = if_else(str_detect(name, "hanam"), "hanam-hi", "kiddyvaxmain"),
      virus = case_when(
        str_detect(name, "bvic") ~ "bvic",
        str_detect(name, "h1pdm") ~ "h1pdm",
        str_detect(name, "h3") ~ "h3",
      ),
      population = if_else(str_detect(name, "exp"), "Exposed", "General")
    )
}

# Save the summary file
save_summ <- function(summ, name, type = "preds") {
  write_csv(
    summ,
    file.path(fit_dir, glue::glue("{name}-{type}-sclr-bayesian.csv"))
  )
}

# Script ======================================================================

# Model output
out_files <- tools::list_files_with_exts(sclr_ba_out_dir, "csv")
out <- map(out_files, read_one) %>%
  map(filter, niter > 6e4) %>%
  map(sample_n, 5e4)
names(out) <- str_replace(basename(out_files), ".csv", "")

# Summary
loghis <- seq(0, 8.7, length.out = 101)
outprobs <- map(out, calc_all_probs, loghis)
outprobs_sum <- map(outprobs, sum_outprob)

# Protective titre levels
prot_his <- out %>%
  map(find_prot_titre_levels_all) %>%
  bind_rows(.id = "name")
save_summ(prot_his, "hanam-hi", "prothis")

# Prior distributions
han_priors_r <- tibble(
  lambda = runif(5e4, 0, 1),
  beta_0 = rnorm(5e4, -15, 10),
  beta_hi = rnorm(5e4, 5, 5)
) %>%
  calc_all_probs(loghis) %>%
  sum_outprob()
names(han_priors_r)[3:5] <- paste0(names(han_priors_r)[3:5], "_prior")

# Add prior distributions
outprobs_sum_prior <- map(
  outprobs_sum, ~ inner_join(.x, han_priors_r, by = c("loghi", "prob_type"))
)

# Split by dataset
outprobs_sum_prior_dt <- outprobs_sum_prior %>%
  imap_dfr(addinfo) %>%
  group_split(dataset)
names(outprobs_sum_prior_dt) <- map_chr(
  outprobs_sum_prior_dt, ~ unique(.x$dataset)
)
outprobs_sum_prior_dt_final <- map(
  outprobs_sum_prior_dt, ~ select(.x, -dataset)
)

iwalk(outprobs_sum_prior_dt_final, save_summ)
