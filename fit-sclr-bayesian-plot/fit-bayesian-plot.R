# Graphing MCMC output
# Arseniy Khvorov
# Created 2019-09-24
# Last edit 2020-01-28

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(facetscales) # devtools::install_github("khvorov45/facetscales")
library(ggrepel)
library(ggpubr)

# Directories to be used later
fit_sclr_ba_summ_dir <- "fit-sclr-bayesian-summary"
fit_sclr_ba_plot_dir <- "fit-sclr-bayesian-plot"
data_dir <- "data"

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

# Read one set of results
read_one_summ <- function(filepath) {
  read_csv(filepath, col_types = cols()) %>% 
    mutate(
      name = str_replace(basename(filepath), ".csv", ""),
      virus = str_match(name, "-([[:alnum:]]*)$")[, 2],
      study = str_match(name, "^([[:alnum:]]*)-")[, 2],
      population = if_else(
        str_detect(name, "-exp-"), "Exposed", "General"
      ) %>% factor(levels = c("General", "Exposed"))
    )
}

recode_viruses <- function(dat) {
  dat %>%
    mutate(
      virus = recode(
        virus, "h1pdm" = "H1N1pdm", "h3" = "H3N2", "bvic" = "B Vic"
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
    geom_line(aes(y = prob_lb_prior), lty = "3333", col = "gray50"),
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
                          xmin = 5, xmax = 1280) {
  facets <- rlang::arg_match(facets, c("vir", "virpop"))
  if (facets == "virpop") {
    faceting <- facet_grid_sc(
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
    )
  } else {
    faceting <- facet_wrap(~virus, nrow = 1)
  }
  data <- filter(data, virus %in% unique(outsum$virus))
  outsum %>%
    ggplot(aes(loghi, prob_med, ymin = prob_lb, ymax = prob_ub)) +
    ylab("Infection probability") +
    coord_cartesian(xlim = c(log(xmin), log(xmax))) +
    faceting +
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

arrange_topbot <- function(top, bot) {
  ggdark::lighten_geoms()
  on.exit(ggdark::darken_geoms())
  top <- top + rremove("x.axis") + rremove("xlab") +
    rremove("x.text") + rremove("x.ticks") +
    theme(
      plot.margin = margin(5.5, 5.5, 0, 5.5),
      axis.ticks.length.x = unit(0, "null")
    )
  top$theme <- ggdark::lighten_theme(top$theme)
  print(top)
  bot <- bot + theme(
    strip.background = element_blank(), strip.text = element_blank(),
    plot.margin = margin(0, 5.5, 5.5, 5.5)
  )
  bot$theme <- ggdark::lighten_theme(bot$theme)
  ggarrange(top, bot, ncol = 1, align = "v")
}

save_plot <- function(pl, name, width = 10, height = 10) {
  plpath <- file.path(fit_sclr_ba_plot_dir, paste0(name, ".pdf"))
  ggsave_dark(
    plot = pl, filename = plpath, dark = FALSE,
    width = width, height = height, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

# Summary data
han_hi_summ <- read_hanam_summ("hanam-hi-summ") %>% filter(virus != "H1N1seas")
kv_main_summ <- read_kv_summ("kiddyvaxmain-summ") %>% recode_viruses()

# Model output summary
out_files <- tools::list_files_with_exts(fit_sclr_ba_summ_dir, "csv")
out_summ <- map_dfr(out_files, read_one_summ) %>%
  recode_viruses()

prot_data <- out_summ %>% filter(prob_type == "prot")
inf_data <- out_summ %>% filter(prob_type == "inf")

# Protection curves
prot_curve_han <- prot_curve_fun(filter(prot_data, study == "hanam"), "virpop")
save_plot(prot_curve_han, "hanam-hi-prot")

prot_curve_han_exp <- prot_curve_fun(
  filter(prot_data, study == "hanam", population == "Exposed"), "vir"
)
save_plot(prot_curve_han_exp, "hanam-hi-exp-prot", 10, 7.5)

prot_curve_kv <- prot_curve_fun(
  filter(prot_data, study == "kiddyvaxmain"), "vir"
)
save_plot(prot_curve_kv, "kiddyvaxmain-prot", 10, 7.5)

# Infection curves
inf_curve_han <- inf_curve_fun(
  filter(inf_data, study == "hanam"), han_hi_summ, "virpop"
)
save_plot(inf_curve_han, "hanam-hi-inf")

inf_curve_han_exp <- inf_curve_fun(
  filter(inf_data, study == "hanam", population == "Exposed"), 
  filter(han_hi_summ, population == "Exposed"), "vir"
)
save_plot(inf_curve_han_exp, "hanam-hi-exp-inf", 10, 7.5)

inf_curve_kv <- inf_curve_fun(
  filter(inf_data, study == "kiddyvaxmain"), kv_main_summ, "vir", xmax = 5120
)
save_plot(inf_curve_kv, "kiddyvaxmain-inf", 10, 7.5)

# Hanam exposed infection and protection
curve_han_exp <- arrange_topbot(inf_curve_han_exp, prot_curve_han_exp)
save_plot(curve_han_exp, "hanam-hi-exp")
