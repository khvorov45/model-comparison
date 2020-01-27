# Graphing MCMC output
# Arseniy Khvorov
# Created 2019/09/24
# Last edit 2019/01/28

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(facetscales) # devtools::install_github("khvorov45/facetscales")
library(ggrepel)

# Directories to be used later
fit_summ_dir <- "fit-bayesian-summary"
fit_plot_dir <- "fit-bayesian-plot"
data_dir <- "data"

# Functions ===================================================================

# Summarise (a subset of) hanam data
sum_han <- function(han) {
  han %>%
    group_by(preHI, virus) %>%
    summarise(
      prop_inf = sum(status != "Not infected", na.rm = TRUE) / n(), n_tot = n()
    ) %>%
    ungroup() %>%
    mutate(
      name = virus, 
      logHImid = case_when(
        near(preHI, 5) ~ log(5),
        near(preHI, 1280) ~ log(1280),
        TRUE ~ log(preHI) + log(2) / 2
      )
    )
}

# Read one set of results
read_one_summ <- function(filepath) {
  read_csv(filepath, col_types = cols()) %>% 
    mutate(
      name = str_replace(basename(filepath), ".csv", ""),
      virus = str_replace(name, "Exp", ""),
      population = if_else(str_detect(name, "Exp"), "Exposed", "General")
    )
}

# Elements common to both plots
common_plot_els <- function() {
  xbreaks <- c(5 * 2^(0:8))
  ybreaks <- seq(0, 1, 0.1)
  list(
    dark_theme_bw(verbose = FALSE),
    theme(
      strip.background = element_rect(fill = NA),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
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

prot_curve_fun <- function(outsum) {
  outsum %>%
    ggplot(aes(logHI, prob_med, ymin = prob_lb, ymax = prob_ub)) +
    ylab("Protection") +
    coord_cartesian(xlim = c(log(5), log(1280))) +
    facet_grid(population ~ virus) +
    geom_line(aes(y = prob_ub_prior), lty = "3333", col = "gray50") +
    common_plot_els()
}

inf_curve_fun <- function(outsum, data) {
  outsum %>%
    ggplot(aes(logHI, prob_med, ymin = prob_lb, ymax = prob_ub)) +
    ylab("Infection probability") +
    coord_cartesian(xlim = c(log(5), log(1280))) +
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
      ),
    ) +
    common_plot_els() +
    geom_point(
      data = data, mapping = aes(logHImid, inf_prop), inherit.aes = FALSE,
      shape = 18
    ) +
    geom_text_repel(
      data = data, mapping = aes(logHImid, inf_prop, label = ntot), 
      inherit.aes = FALSE
    )
}

save_plot <- function(pl, name, dark) {
  suff <- if_else(dark, "_dark", "_light")
  plpath <- file.path(fit_plot_dir, paste0(name, suff, ".pdf"))
  ggsave_dark(
    plot = pl, filename = plpath, dark = dark,
    width = 10, height = 10, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

# Hanam data
han <- read_csv(file.path(data_dir, "hanam-HI-summ.csv")) %>%
  filter(virus != "H1N1seas") %>%
  mutate(
    logHImid = case_when(
      near(preHI, 5) ~ log(5),
      near(preHI, 1280) ~ log(1280),
      TRUE ~ log(preHI) + log(2) / 2
    ),
    population = factor(
      population, levels = c("general", "exposed"), 
      labels = c("General", "Exposed")
    )
  )

# Model output summary
out_files <- tools::list_files_with_exts(fit_summ_dir, "csv")
out_summ <- map_dfr(out_files, read_one_summ) %>%
  mutate(population = factor(population, levels = c("General", "Exposed")))

# Curves

prot_data <- out_summ %>% filter(prob_type == "prot")
inf_data <- out_summ %>% filter(prob_type == "inf")

prot_curve <- prot_curve_fun(prot_data)
inf_curve <- inf_curve_fun(inf_data, han)

walk(
  c(TRUE, FALSE), 
  ~ walk2(
    list(prot_curve, inf_curve), c("protection", "infection"), save_plot, .x
  )
)
