# Plots of the ordinary logistic fit
# Arseniy Khvorov
# Created 2020-01-28
# Last edit 2020-02-04

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
      virus = factor(
        virus,
        levels = c("h1pdm", "bvic"), labels = c("H1N1pdm", "B Vic")
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

convert_to_prot <- function(dat) {
  dat %>%
    mutate(
      fit = 1 - fit,
      fit_low = 1 - fit_low,
      fit_high = 1 - fit_high
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

prot_curve_fun <- function(outsum, facets = "virpop", xmin = 5, xmax = 1280) {
  facets <- rlang::arg_match(facets, c("vir", "virpop", "pop", "none"))
  if (facets == "virpop") {
    faceting <- facet_grid(population ~ virus)
  } else if (facets == "vir") {
    faceting <- facet_wrap(~virus, nrow = 1)
  } else if (facets == "pop") {
    faceting <- facet_wrap(~population, nrow = 1)
  } else {
    faceting <- NULL
  }
  outsum %>%
    ggplot(aes(loghimid, fit, ymin = fit_low, ymax = fit_high)) +
    ylab("Protection") +
    coord_cartesian(xlim = c(log(xmin), log(xmax)), ylim = c(0, 1)) +
    faceting +
    common_plot_els()
}

inf_curve_fun <- function(outsum, data, facets = "virpop",
                          xmin = 5, xmax = 5120, ymin = 0, ymax = 0.5) {
  facets <- rlang::arg_match(facets, c("vir", "virpop", "pop", "none"))
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

han_hi_lr_prot <- convert_to_prot(han_hi_lr)
kvm_lr_prot <- convert_to_prot(kvm_lr)

# Hanam plots
han_hi_lr_plot <- inf_curve_fun(han_hi_lr, han_hi_summ, "virpop", 5, 1280)
# save_plot(han_hi_lr_plot, "hanam-hi-inf")

han_hi_lr_prot_plot <- prot_curve_fun(han_hi_lr_prot, "virpop", 5, 1280)
# save_plot(han_hi_lr_prot_plot, "hanam-hi-prot")

han_h3_lr_prot_plot <- prot_curve_fun(
  filter(han_hi_lr_prot, virus == "H3N2"), "pop", 5, 1280
)
save_plot(han_h3_lr_prot_plot, "hanam-h3-prot", 12, 7.5)

# Kiddyvax plots
kvm_lr_plot <- inf_curve_fun(kvm_lr, kv_main_summ, "vir", 5, 5120, 0, 0.27)
# save_plot(kvm_lr_plot, "kiddyvaxmain-inf", 10, 7.5)

kvm_lr_prot_plot <- prot_curve_fun(kvm_lr_prot, "vir", 5, 5120)
# save_plot(kvm_lr_prot_plot, "kiddyvaxmain-prot", 10, 7.5)

kvm_lr_prot_plot_bvic <- prot_curve_fun(
  filter(kvm_lr_prot, virus == "B Vic"), "none", 5, 5120
)
save_plot(kvm_lr_prot_plot_bvic, "kiddyvaxmain-bvic-prot", 7.5, 7.5)
