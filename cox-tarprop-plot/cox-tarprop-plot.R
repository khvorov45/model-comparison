# Graphing of simulation results
# Arseniy Khvorov
# Created 2019-09-20
# Last edit 2020-03-03

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories to be used later
cox_tarprop_summ_dir <- "cox-tarprop-summary"
cox_tarprop_plot_dir <- "cox-tarprop-plot"

# Settings ===================================================================

y_var_labs <- c(
  "est_mean" = "Mean estimate",
  "rel_bias" = "Relative bias",
  "est_sd" = "Estimate SD"
)

# Functions ===================================================================

read_res <- function(name) {
  read_csv(
    file.path(cox_tarprop_summ_dir, glue::glue("{name}.csv")),
    col_types = cols_only(
      est_mean = col_double(),
      est_sd = col_double(),
      true_val = col_double(),
      par_varied = col_character(),
      risk_prop_expected = col_double(),
      long_prop_expected = col_double()
    )
  )
}

add_hline <- function(name, val, y_vars) {
  geom_hline(
    data = tibble(name = factor(name, levels = y_vars), val = val),
    mapping = aes(yintercept = val),
    lty = "1111"
  )
}

add_plot_els <- function(y_vars, true_value_extr) {
  add_els <- list()
  if ("est_mean" %in% y_vars) {
    add_els <- c(add_els, list(
      add_hline("est_mean", true_value_extr, y_vars)
    ))
  }
  if ("rel_bias" %in% y_vars) {
    add_els <- c(add_els, list(
      add_hline("rel_bias", 0, y_vars)
    ))
  }
  add_els
}

plot_fun <- function(res, x_name, x_lab, y_vars = names(y_var_labs)) {
  res %>%
    mutate(rel_bias = (est_mean - true_val) / abs(true_val)) %>%
    filter(par_varied == tidyselect::all_of(x_name)) %>%
    select(
      tidyselect::all_of(x_name), true_val, tidyselect::all_of(y_vars)
    ) %>%
    pivot_longer(
      tidyselect::all_of(y_vars), names_to = "name", values_to = "value"
    ) %>%
    mutate(name = factor(name, levels = y_vars)) %>%
    ggplot(aes(!!sym(x_name), value)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(0, "null"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    facet_wrap(
      ~name, scales = "free_y", nrow = 1, strip.position = "left",
      labeller = as_labeller(y_var_labs)
    ) +
    xlab(x_lab) +
    ylab(NULL) +
    scale_x_continuous(
      breaks = unique(res[[x_name]]), labels = scales::percent_format(1)
    ) +
    geom_line() +
    geom_point(shape = 18) +
    add_plot_els(y_vars, first(res$true_val))
}

save_plot <- function(pl, name) {
  ggsave_dark(
    pl, 
    filename = file.path(cox_tarprop_plot_dir, paste0(name, ".pdf")),
    dark = FALSE, width = 15, height = 7, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

summ <- read_res("summary-10000sims")

pl_risk <- plot_fun(
  summ, 
  "risk_prop_expected", "Expected proportion of time at risk",
  c("rel_bias", "est_sd")
)
save_plot(pl_risk, "risk")

pl_long <- plot_fun(
  summ, 
  "long_prop_expected", "Expected proportion with earlier follow-up start",
  c("rel_bias", "est_sd")
)
save_plot(pl_long, "long")
