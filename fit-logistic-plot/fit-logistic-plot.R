# Plots of the ordinary logistic fit
# Arseniy Khvorov
# Created 2020-01-28
# Last edit 2020-02-03

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(facetscales) # devtools::install_github("khvorov45/facetscales")
library(ggrepel)

# Directories to be used later
fit_logistic_dir <- "fit-logistic"
data_dir <- "data"
fit_logistic_plot_dir <- "fit-logistic-plot"

# Functions ===================================================================

read_hanam_summ <- function(name) {
  read_csv(
    file.path(data_dir, paste0(name, ".csv")),
    col_types = cols(
      population = col_factor(levels = c("General", "Exposed")),
      virus = col_character(),
      prehi = col_integer(),
      loghimid = col_double(),
      ntot = col_integer(),
      inf_prop = col_double()
    )
  )
}

read_kv_summ <- function(name) {
  read_csv(
    file.path(data_dir, paste0(name, ".csv")),
    col_types = cols(
      virus = col_character(),
      hi = col_integer(),
      loghimid = col_double(),
      ntot = col_integer(),
      inf_prop = col_double()
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

read_hanam_lrres <- function(name) {
  read_csv(
    file.path(fit_logistic_dir, paste0(name, ".csv")),
    col_types = cols_only(
      fit = col_double(),
      fit_high = col_double(),
      fit_low = col_double(),
      loghimid = col_double(),
      population = col_factor(levels = c("General", "Exposed")),
      virus = col_character()
    )
  )
}

read_kv_lrres <- function(name) {
  read_csv(
    file.path(fit_logistic_dir, paste0(name, ".csv")),
    col_types = cols_only(
      fit = col_double(),
      fit_high = col_double(),
      fit_low = col_double(),
      loghimid = col_double(),
      virus = col_character()
    )
  )
}

# Elements common to both plots
common_plot_els <- function() {
  xbreaks <- c(5 * 2^(0:10))
  ybreaks <- seq(0, 1, 0.1)
  list(
    dark_theme_bw(verbose = FALSE),
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.y = element_blank()
    ),
    xlab("HI titre"),
    scale_x_continuous(
      breaks = log(xbreaks), labels = xbreaks,
      minor_breaks = log((xbreaks + xbreaks * 2) / 2)
    ),
    scale_y_continuous(breaks = ybreaks, labels = scales::percent_format(1)),
    geom_ribbon(alpha = 0.5),
    geom_line()
  )
}

prot_curve_fun <- function(outsum, facets = "virpop") {
  facets <- rlang::arg_match(facets, c("vir", "virpop"))
  if (facets == "virpop") {
    faceting <- facet_grid(population ~ virus)
  } else {
    faceting <- facet_wrap(~virus, nrow = 1)
  }
  outsum %>%
    ggplot(aes(loghi, prob_med, ymin = prob_lb, ymax = prob_ub)) +
    ylab("Protection") +
    coord_cartesian(xlim = c(log(5), log(1280))) +
    faceting +
    geom_line(aes(y = prob_ub_prior), lty = "3333", col = "gray50") +
    common_plot_els()
}

inf_curve_fun <- function(outsum, data, facets = "virpop", 
                          xmin = 5, xmax = 5120, ymin = 0, ymax = 0.5) {
  facets <- rlang::arg_match(facets, c("vir", "virpop"))
  if (facets == "virpop") {
    facets_spec <- list(
      facet_grid_sc(
        population ~ virus, 
        scales = list(
          y = list(
            General = scale_y_continuous(
              limits = c(0, 0.3), labels = scales::percent_format(1)
            ),
            Exposed = scale_y_continuous(
              limits = c(0, 0.5), labels = scales::percent_format(1)
            )
          )
        )
      ),
      coord_cartesian(xlim = c(log(xmin), log(xmax)))
    )
  } else {
    facets_spec <- list(
      facet_wrap(~virus, nrow = 1),
      coord_cartesian(xlim = c(log(xmin), log(xmax)), ylim = c(ymin, ymax))
    )
  }
  data <- filter(data, virus %in% unique(outsum$virus))
  outsum %>%
    ggplot(aes(loghimid, fit, ymin = fit_low, ymax = fit_high)) +
    ylab("Infection probability") +
    facets_spec +
    common_plot_els() +
    geom_point(
      data = data, mapping = aes(loghimid, inf_prop), inherit.aes = FALSE,
      shape = 18
    ) +
    geom_text_repel(
      data = data, mapping = aes(loghimid, inf_prop, label = ntot), 
      inherit.aes = FALSE
    )
}

save_plot <- function(pl, name, width = 10, height = 10) {
  plpath <- file.path(fit_logistic_plot_dir, paste0(name, ".pdf"))
  ggsave_dark(
    plot = pl, filename = plpath, dark = FALSE,
    width = width, height = height, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

# Data
han_hi_summ <- read_hanam_summ("hanam-hi-summ") %>% filter(virus != "H1N1seas")
kv_main_summ <- read_kv_summ("kiddyvaxmain-summ") %>% recode_viruses()

# Logistic fit predictions
han_hi_lr <- read_hanam_lrres("hanam-hi")
kvm_lr <- read_kv_lrres("kiddyvaxmain") %>% recode_viruses()

# Hanam plots
han_hi_lr_plot <- inf_curve_fun(han_hi_lr, han_hi_summ, "virpop", 5, 1280)
save_plot(han_hi_lr_plot, "hanam-hi")

# Kiddyvax plots
kvm_lr_plot <- inf_curve_fun(kvm_lr, kv_main_summ, "vir", 5, 5120, 0, 0.3)
save_plot(kvm_lr_plot, "kiddyvaxmain", 12.5, 8.5)
