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
    col_types = cols()
  )
}

recode_viruses <- function(dat) {
  dat %>%
    mutate(
      virus = factor(
        virus,
        levels = c("h1pdm", "bvic"), labels = c("H1N1pdm", "B Vic")
      )
    )
}

plot_pred <- function(dat, facet_by_virus = TRUE) {
  xbreaks <- c(5 * 2^(0:10))
  facets <- if (facet_by_virus) facet_wrap(~virus, nrow = 1) else NULL
  dat %>%
    ggplot(aes(loghi, prot)) +
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
    facets +
    scale_x_continuous("HI titre", breaks = log(xbreaks), labels = xbreaks) +
    scale_y_continuous(
      "Relative protection",
      breaks = seq(0, 1, 0.1),
      labels = scales::percent_format(1)
    ) +
    geom_ribbon(aes(ymin = prot_low, ymax = prot_high), alpha = 0.5) +
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

preds_sophia <- read_pred("sophia")

pl <- plot_pred(preds)
# save_plot(pl, "kiddyvaxmain", 12, 7.5)

pl_bvic <- plot_pred(filter(preds, virus == "B Vic"), facet_by_virus = FALSE)
save_plot(pl_bvic, "kiddyvaxmain-bvic", 7.5, 7.5)

pl_soph_og <- plot_pred(
  filter(preds_sophia, model == "sophia") %>%
    mutate(prot_low = prot_low_wrong, prot_high = prot_high_wrong)
)
# save_plot(pl_soph_og, "sophia-og", 12, 7.5)

pl_soph_ci <- plot_pred(
  filter(preds_sophia, model == "sophia")
)
# save_plot(pl_soph_ci, "sophia-ci", 12, 7.5)

pl_soph_ci_mod <- plot_pred(
  filter(preds_sophia, model == "me")
)
# save_plot(pl_soph_ci_mod, "sophia-ci-mod", 12, 7.5)
