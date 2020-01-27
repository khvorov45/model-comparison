# Graphing MCMC output
# Arseniy Khvorov
# Created 2019/09/24
# Last edit 2019/12/12

library(tidyverse)
library(ggpubr)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories to be used later
fit_bayesian_dir <- "fit-bayesian"
fit_diag_dir <- "fit-bayesian-diag"

# Functions ===================================================================

# Read one set of results
read_one <- function(filepath) {
  read_csv(filepath, col_types = cols())
}

# Converts output to long format
out_to_long <- function(tidyout) {
  tidyout %>%
    pivot_longer(
      c(-niter, -nchain), names_to = "parameter", values_to = "value"
    )
}

# Common plot theme
common_theme <- function() {
  dark_theme_bw(verbose = FALSE) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "null"),
    panel.spacing = unit(0, "null"),
    strip.background = element_rect(fill = NA),
    axis.text.y = element_text(angle = 90)
  )
}

trace_plot <- function(tidyout) {
  out_to_long(tidyout) %>%
    ggplot(aes(value, niter, col = as.factor(nchain))) +
    common_theme() +
    theme(
      axis.text.x = element_text(angle = 90),
    ) +
    scale_y_continuous(
      "Iteration", breaks = c(min(tidyout$niter), max(tidyout$niter)),
      labels = scales::comma_format()
    ) +
    xlab("Value") +
    facet_wrap(~ parameter, nrow = 1, scale = "free_x") +
    geom_path(alpha = 0.5)
}

add_prior <- function(out_par, prior_dists) {
  prior_dist <- prior_dists[[unique(out_par$parameter)]]
  mutate(out_par, prior = prior_dist(value))
}

density_plot <- function(tidyout, prior_dists) {
  tidyout_long <- out_to_long(tidyout) %>%
    group_split(parameter) %>%
    map_dfr(add_prior, prior_dists = prior_dists) %>%
    ggplot(aes(value, col = as.factor(nchain))) +
    common_theme() +
    theme(
      axis.text.y = element_text(angle = 90, hjust = 0),
      axis.ticks.length.x = unit(0, "null")
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    xlab("Value") +
    ylab("Density") +
    facet_wrap(~ parameter, nrow = 1, scales = "free") +
    geom_freqpoly(alpha = 0.5, stat = "density") +
    geom_line(aes(x = value, y = prior), lty = "3131", col = "green", lwd = 1.1)
}

arrange_itdens_one <- function(dens, iter) {
  ggarrange(
    dens + rremove("x.axis") + rremove("xlab") +
      rremove("x.text") + rremove("x.ticks"),
    iter + theme(
      strip.background = element_blank(), strip.text = element_blank()
    ),
    ncol = 1, align = "v"
  )
}

save_trace_dens <- function(pl, model_name) {
  ggsave_dark(
    pl, 
    filename = file.path(fit_diag_dir, paste0(model_name, ".pdf")),
    dark = FALSE,
    width = 20, height = 8, units = "cm",
    device = "pdf"
  )
}

# Script ======================================================================

# Prior distributions
han_priors <- list(
  logHImean = function(x) dnorm(x, 2, 2),
  logHIsd = function(x) dexp(x, rate = 0.1),
  lambda = function(x) dunif(x, 0, 1),
  beta_0 = function(x) dnorm(x, -15, 10),
  beta_HI = function(x) dnorm(x, 5, 5)
)

# Model output
out_files <- tools::list_files_with_exts(fit_bayesian_dir, "csv")
out <- map(out_files, read_one)
names(out) <- str_replace(basename(out_files), ".csv", "")

# Trace
trace_plots <- map(out, trace_plot)

# Density
density_plots <- map(out, density_plot, han_priors)

trace_density_ar <- map2(density_plots, trace_plots, arrange_itdens_one)

# Save plots
iwalk(trace_density_ar, save_trace_dens)
