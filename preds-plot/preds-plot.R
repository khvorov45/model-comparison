# Plotting protections curves

library(tidyverse)

# Directories used
fit_dir <- here::here("fit")
preds_plot_dir <- here::here("preds-plot")

# Functions ===================================================================

read_pred <- function(name) {
  read_csv(
    file.path(fit_dir, glue::glue("{name}.csv")),
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
    ggdark::dark_theme_bw(verbose = FALSE) +
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

save_plot <- function(pl, name, width = 12, height = 7.5) {
  ggdark::ggsave_dark(
    file.path(preds_plot_dir, glue::glue("{name}.pdf")), pl,
    device = "pdf", width = width, height = height, units = "cm"
  )
}

# Script ======================================================================

kv_cox_preds <- read_pred("kiddyvaxmain-preds-cox") %>% recode_viruses()
sophia_preds <- read_pred("sophia-preds-cox")

all_plots <- list(
  "kiddyvaxmain-cox" = plot_pred(kv_cox_preds),
  "kiddyvaxmain-cox-bvic" = plot_pred(
    filter(kv_cox_preds, virus == "B Vic"),
    facet_by_virus = FALSE
  ),
  "sophia-cox-og" = plot_pred(
    filter(sophia_preds, model == "sophia") %>%
      mutate(prot_low = prot_low_wrong, prot_high = prot_high_wrong)
  ),
  "sophia-cox-fixci" = plot_pred(
    filter(preds_sophia, model == "sophia")
  ),
  "sophia-cox-fixci-fixmod" = plot_pred(
    filter(preds_sophia, model == "me")
  )
)

iwalk(all_plots, save_plot)
