# Plausible parameter ranges

library(tidyverse)

# Directories used
plausible_dir <- here::here("plausible")
sclr_dir <- file.path(plausible_dir, "sclr")

# Functions ===================================================================

# Grid of expected curves at different parameter values
expected_lines <- function(dat_pars) {
  breaks <- 5 * 2^(0:8)
  dat_pars %>%
    ggplot(aes(exp(x), exp_prop)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      axis.text.x = element_text(angle = 90)
    ) +
    scale_x_log10("HI titre", breaks = breaks, labels = breaks) +
    scale_y_continuous(
      "Expected infected proportion",
      breaks = c(0.05, 0.1, 0.15)
    ) +
    coord_cartesian(xlim = c(5, 1280), ylim = c(0, 0.2)) +
    facet_grid(
      rows = vars(b1), cols = vars(b0),
      labeller = label_bquote(rows = b1 == .(b1), cols = b0 == .(b0))
    ) +
    geom_line() +
    ggtitle(bquote(lambda == .(unique(dat_pars$lmd))))
}

# Creates all the combinations
create_combos <- function(uniq) {
  combos <- tibble(Int = 1)
  for (nm in names(uniq)) {
    vals <- unique(uniq[[nm]])
    times_to_rep <- nrow(combos)
    combos <- slice(combos, rep(1:n(), each = length(vals)))
    combos[[nm]] <- rep(vals, times_to_rep)
  }
  select(combos, -Int)
}

# Expected curve from fitting the scaled logit model
expected_curve <- function(x, b0, bloghi, lmd) {
  lmd / (1 + exp(b0 + bloghi * x))
}

# Expected titre proportions
expected_int_proportion <- function(x, mn, sd) {
  this_dnorm <- function(x) dnorm(x, mean = mn, sd = sd)
  case_when(
    near(x, 5) ~ integrate(this_dnorm, -Inf, log(10))$value,
    near(x, 1280) ~ integrate(this_dnorm, log(1280), Inf)$value,
    TRUE ~ integrate(this_dnorm, log(x), log(2 * x))$value
  )
}
expected_int_proportion_vec <- function(x_vec, mn_vec, sd_vec) {
  pmap_dbl(list(x_vec, mn_vec, sd_vec), expected_int_proportion)
}

# Plot of titre distributions at different parameters
plot_titre_dist <- function(dat) {
  plbreaks <- 5 * 2^(0:8)
  titre_combs %>%
    ggplot(aes(x, exp_prop)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "lines")
    ) +
    coord_cartesian(xlim = c(5, 1280), ylim = c(0, 0.75)) +
    scale_x_log10(
      "HI titre",
      breaks = plbreaks, labels = plbreaks
    ) +
    ylab("Proportion infected") +
    facet_grid(
      cols = vars(mn), rows = vars(sd),
      labeller = label_bquote(rows = "sd" == .(sd), cols = "mean" == .(mn))
    ) +
    geom_step(col = "white")
}

# Saves the various parameter plots
save_par_plots <- function(pl, ind) {
  if (!dir.exists(sclr_dir)) dir.create(sclr_dir)
  plpath <- file.path(sclr_dir, paste0("pl_", ind, ".png"))
  ggdark::ggsave_dark(
    pl,
    dark = TRUE, filename = plpath, device = "png",
    width = 20, height = 20, units = "cm"
  )
}

# Script ======================================================================

# Expected lines for different parameter sets
pars <- list(
  lmd = seq(0.05, 1, 0.05),
  b0 = seq(-30, 2, 8),
  b1 = seq(0, 10, 2),
  x = seq(0, 7, length.out = 101)
)

pars_combs <- create_combos(pars) %>%
  mutate(exp_prop = expected_curve(x, b0, b1, lmd)) %>%
  group_split(lmd)

pars_pls <- map(pars_combs, expected_lines)
iwalk(pars_pls, save_par_plots)

# Expected titre distributions at different means and sds
titre_dist <- list(
  mn = seq(-1, 4, 1),
  sd = seq(0.5, 5, 0.5),
  x = 5 * 2^(0:8)
)

titre_combs <- create_combos(titre_dist) %>%
  mutate(exp_prop = expected_int_proportion_vec(x, mn, sd))

prop_inf <- plot_titre_dist(titre_combos)

ggdark::ggsave_dark(
  prop_inf,
  dark = FALSE,
  filename = file.path(plausible_dir, "plausible-titres.pdf"),
  width = 20, height = 20, units = "cm"
)
