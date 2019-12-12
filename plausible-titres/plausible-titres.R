# Plausible titre distribution
# Arseniy Khvorov
# Created 2019/09/17
# Last edit 2019/12/06

library(tidyverse)
library(ggdark) # devtools:install_github("khvorov45/ggdark")

# Directories to be used later
plaus_titres_dir <- "plausible-titres"

# Functions ===================================================================

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
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.spacing = unit(0, "lines"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "lines")
    ) +
    coord_cartesian(xlim = c(5, 1280), ylim = c(0, 0.75)) +
    scale_x_log10(
      "HI titre", breaks = plbreaks, labels = plbreaks
    ) +
    ylab("Proportion infected") +
    facet_grid(
      cols = vars(mn), rows = vars(sd),
      labeller = label_bquote(rows = "sd" == .(sd), cols = "mean" == .(mn))
    ) +
    geom_step(col = "white")
}

# Script ======================================================================

# Expected titre distributions at different means and sds
titre_dist <- list(
  mn = seq(-1, 4, 1),
  sd = seq(0.5, 5, 0.5),
  x = 5 * 2^(0:8)
)

titre_combs <- create_combos(titre_dist) %>%
  mutate(exp_prop = expected_int_proportion_vec(x, mn, sd))

prop_inf <- plot_titre_dist(titre_combos)

ggsave_dark(
  prop_inf, dark = TRUE, 
  filename = file.path(plaus_titres_dir, "plausible-titres.pdf"),
  width = 20, height = 20, units = "cm"
)
