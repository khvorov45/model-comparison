# Plausible values of sclr model parameters
# Arseniy Khvorov
# Created 2019/09/17
# Last edit 2019/12/06

library(tidyverse)
library(ggdark) # devtools:install_github("khvorov45/ggdark")

# Directories to be used later
plaus_sclr_dir <- "plausible-sclr"

# Functions ===================================================================

# Grid of expected curves at different parameter values
expected_lines <- function(dat_pars) {
  breaks = 5 * 2^(0:8)
  dat_pars %>%
    ggplot(aes(exp(x), exp_prop)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      axis.text.x = element_text(angle = 90)
    ) +
    scale_x_log10("HI titre", breaks = breaks, labels = breaks) +
    scale_y_continuous(
      "Expected infected proportion", breaks = c(0.05, 0.1, 0.15)
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
    combos <- slice(combos, rep(1:n(), each = length(vals)))
    combos[[nm]] <- uniq[[nm]]
  }
  select(combos, -Int)
}

# Expected curve from fitting the scaled logit model
expected_curve <- function(x, b0, bloghi, lmd) {
  lmd / (1 + exp(b0 + bloghi * x))
}

# Saves the various parameter plots
save_par_plots <- function(pl, ind) {
  plpath = file.path(plaus_sclr_dir, paste0("pl_", ind, ".png"))
  ggsave_dark(
    pl, dark = TRUE, filename = plpath, device = "png",
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
