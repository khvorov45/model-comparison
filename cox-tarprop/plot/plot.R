# Graphing of simulation results

library(tidyverse)
library(latex2exp)

# Directories used
summ_dir <- here::here("summ")
plot_dir <- here::here("plot")

# Settings ===================================================================

y_var_labs <- c(
  "est_mean" = "Mean estimate",
  "rel_bias" = "Relative bias",
  "est_sd" = "Estimate SD"
)

pop_labs <- c(
  "std_05" = "$\\kappa$ = 0.5",
  "std_1" = "$\\kappa$ = 1",
  "std_10" = "$\\kappa$ = 10",
  "std_100" = "$\\kappa$ = 100",
  "std_1000" = "$\\kappa$ = 1000"
)

pops <- c("std_1", "std_10", "std_100")

# Functions ===================================================================

label_facets <- function(lbls) {
  if ("data_name" == names(lbls)) {
    lbls$data_name <- recode(
      lbls$data_name, !!!pop_labs
    ) %>% TeX()
  } else {
    lbls$name <- recode(
      lbls$name, !!!y_var_labs
    ) %>% TeX()
  }
  label_parsed(lbls)
}

label_kappa <- function(brks) {
  brks <- str_replace(brks, "std_", "")
  recode(brks, "05" = "0.5")
}

add_hline <- function(name, val, y_vars) {
  geom_hline(
    data = tibble(name = factor(name, levels = y_vars), val = val),
    mapping = aes(yintercept = val),
    lty = "1111"
  )
}

plot_fun <- function(res, x_name, x_lab, y_vars = names(y_var_labs),
                     pops = names(pop_labs)) {
  if (length(unique(res$data_name)) > 1L) {
    my_ggplot <- function(dat) {
      ggplot(dat, aes(!!sym(x_name), value, color = data_name))
    }
  } else {
    my_ggplot <- function(dat) {
      ggplot(dat, aes(!!sym(x_name), value))
    }
  }
  res %>%
    filter(data_name %in% pops) %>%
    mutate(rel_bias = (est_mean - true_val) / abs(true_val)) %>%
    filter(par_varied == tidyselect::all_of(x_name)) %>%
    select(
      tidyselect::all_of(x_name), true_val, tidyselect::all_of(y_vars),
      data_name
    ) %>%
    pivot_longer(
      tidyselect::all_of(y_vars),
      names_to = "name", values_to = "value"
    ) %>%
    mutate(name = factor(name, levels = y_vars)) %>%
    my_ggplot() +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.x = unit(0, "null"),
      panel.spacing.y = unit(0.4, "lines"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    facet_grid(
      name ~ .,
      scales = "free_y",
      labeller = label_facets, switch = "y"
    ) +
    xlab(x_lab) +
    ylab(NULL) +
    scale_x_continuous(
      breaks = unique(res[[x_name]]), labels = scales::percent_format(1)
    ) +
    scale_color_discrete(TeX("$\\kappa$"), labels = label_kappa) +
    geom_line() +
    geom_point(shape = 18)
}

save_plot <- function(pl, name, width = 15, height = 7) {
  ggdark::ggsave_dark(
    pl,
    filename = file.path(plot_dir, paste0(name, ".pdf")),
    dark = FALSE, width = width, height = height, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

summ <- read_csv(file.path(summ_dir, "summ.csv"), col_types = cols())

pl_risk <- plot_fun(
  summ,
  "risk_prop_expected", "Expected proportion of time at risk",
  c("rel_bias", "est_sd")
)
save_plot(pl_risk, "risk")

pl_long <- plot_fun(
  filter(summ, data_name == "std_05"),
  "long_prop_expected", "Expected proportion with earlier follow-up start",
  c("rel_bias", "est_sd")
)
save_plot(pl_long, "long")
