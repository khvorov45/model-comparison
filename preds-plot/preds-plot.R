# Plotting protections curves

library(tidyverse)

# Directories used
fit_dir <- here::here("fit")
preds_plot_dir <- here::here("preds-plot")
data_dir <- here::here("data")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

read_pred <- function(name) {
  read_csv(
    file.path(fit_dir, glue::glue("{name}.csv")),
    col_types = cols()
  )
}

recode_viruses <- function(dat) {
  dat %>%
    mutate(
      virus_lbl = factor(
        virus,
        levels = c("H1N1pdm", "h1pdm", "H3N2", "bvic"),
        labels = c("H1N1pdm", "H1N1pdm", "H3N2", "B Vic")
      )
    )
}

plot_pred <- function(dat, facets = "vir",
                      xmin = 5, xmax = 5120, ylab = "Relative protection") {
  xbreaks <- c(5 * 2^(0:10))
  facets <- rlang::arg_match(facets, c("vir", "virpop", "pop", "none"))
  if (facets == "virpop") {
    faceting <- facet_grid(population ~ virus_lbl)
  } else if (facets == "vir") {
    faceting <- facet_wrap(~virus_lbl, nrow = 1)
  } else if (facets == "pop") {
    faceting <- facet_wrap(~population, nrow = 1)
  } else {
    faceting <- NULL
  }
  dat %>%
    ggplot(aes(loghi, prot)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    faceting +
    coord_cartesian(xlim = c(log(xmin), log(xmax)), ylim = c(0, 1)) +
    scale_x_continuous("HI titre", breaks = log(xbreaks), labels = xbreaks) +
    scale_y_continuous(
      ylab,
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

kv_main_summ <- read_data("kiddyvaxmain-summ")
han_hi_summ <- read_data("hanam-hi-summ") %>%
  filter(virus != "H1N1seas") %>%
  recode_viruses()

kv_cox_preds <- read_pred("kiddyvaxmain-preds-cox")
sophia_preds <- read_pred("sophia-preds-cox")
han_hi_lr <- read_pred("hanam-hi-preds-lr") %>% recode_viruses()
kvm_lr <- read_pred("kiddyvaxmain-preds-lr") %>% recode_viruses()

all_plots <- list(
  "kiddyvaxmain-cox" = plot_pred(kv_cox_preds),
  "kiddyvaxmain-lr" = plot_pred(kvm_lr, ylab = "Protection"),
  "kiddyvaxmain-cox-bvic" = plot_pred(
    filter(kv_cox_preds, virus_lbl == "B Vic"),
    facets = "none",
  ),
  "sophia-cox-og" = plot_pred(
    filter(sophia_preds, model == "sophia") %>%
      mutate(prot_low = prot_low_wrong, prot_high = prot_high_wrong)
  ),
  "sophia-cox-fixci" = plot_pred(
    filter(sophia_preds, model == "sophia")
  ),
  "sophia-cox-fixci-fixmod" = plot_pred(
    filter(sophia_preds, model == "me")
  ),
  "hanam-hi-lr" = plot_pred(han_hi_lr, "virpop", ylab = "Protection")
)

iwalk(all_plots, save_plot)
