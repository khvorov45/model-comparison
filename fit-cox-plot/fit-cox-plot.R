# Plotting cox predictions
# Arseniy Khvorov
# Created 2020-02-17
# Last edit 2020-02-17

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories used
fit_cox_dir <- "fit-cox"
fit_cox_plot_dir <- "fit-cox-plot"

# Functions ===================================================================

read_pred <- function(name) {
  read_csv(
    file.path(fit_cox_dir, glue::glue("{name}.csv")),
    col_types = cols_only(
      virus = col_character(),
      loghimid = col_double(),
      prot_rel = col_double(),
      prot_rel_low = col_double(),
      prot_rel_high = col_double()
    )
  )
}

recode_viruses <- function(dat) {
  dat %>%
    mutate(
      virus = recode(
        virus, "h1pdm" = "H1N1pdm", "bvic" = "B Vic"
      )
    )
}

plot_pred <- function(dat) {
  xbreaks <- c(5 * 2^(0:10))
  ggplot(dat, aes(loghimid, prot_rel)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(
      xlim = c(log(5), log(5120)),
      ylim = c(0, 1)
    ) +
    facet_wrap(~ virus, nrow = 1) +
    scale_x_continuous("HI titre", breaks = log(xbreaks), labels = xbreaks) +
    scale_y_continuous(
      "Relative protection", breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(1)
    ) +
    geom_ribbon(aes(ymin = prot_rel_low, ymax = prot_rel_high), alpha = 0.5) +
    geom_line()
}

save_plot <- function(pl, name, width, height) {
  ggsave_dark(
    file.path(fit_cox_plot_dir, glue::glue("{name}.pdf")), pl,
    device = "pdf", width = width, height = height, units = "cm"
  )
}

# Script ======================================================================

preds <- read_pred("kiddyvaxmain") %>% recode_viruses()

pl <- plot_pred(preds)
save_plot(pl, "kiddyvaxmain", 12, 7.5)
