# Graphing of simulation results
# Arseniy Khvorov
# Created 2019/09/20
# Last edit 2019/12/02

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(ggpubr)

# Directories to be used later
cox_tarprop_summ_dir <- "cox-tarprop-summary"
cox_tarprop_plot_dir <- "cox-tarprop-plot"

# Functions ===================================================================

plot_fun <- function(res, x_name, x_lab) {
  true_value_extr <- first(res$true_val)
  res %>%
    filter(par_varied == x_name) %>%
    select(x_name, est_mean, est_sd, true_val) %>%
    pivot_longer(
      c(est_mean, est_sd), names_to = "name", values_to = "value"
    ) %>%
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
      labeller = as_labeller(
        c("est_mean" = "Estimate mean", "est_sd" = "Estimate SD")
      )
    ) +
    xlab(x_lab) +
    ylab(NULL) +
    scale_x_continuous(
      breaks = unique(res[[x_name]]), labels = scales::percent_format(1)
    ) +
    geom_line() +
    geom_point(shape = 18) +
    geom_hline(
      data = tibble(name = "est_mean", val = true_value_extr),
      mapping = aes(yintercept = val),
      lty = "1111"
    )
}

save_plot <- function(pl, name) {
  ggsave_dark(
    pl, 
    filename = file.path(cox_tarprop_plot_dir, paste0(name, ".pdf")),
    dark = FALSE, width = 15, height = 7, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

summ <- read_csv(file.path(cox_tarprop_summ_dir, "summary-10000sims.csv"))

pl_risk <- plot_fun(
  summ, 
  "risk_prop_expected", "Expected proportion of time at risk"
)
save_plot(pl_risk, "risk")

pl_long <- plot_fun(
  summ, 
  "long_prop_expected", "Expected proportion with earlier follow-up start"
)
save_plot(pl_long, "long")
