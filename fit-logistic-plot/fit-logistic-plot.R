# Plots of the ordinary logistic fit
# Arseniy Khvorov
# Created 2020-01-28
# Last edit 2020-01-28

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(facetscales) # devtools::install_github("khvorov45/facetscales")
library(ggrepel)

# Directories to be used later
fit_logistic_dir <- "fit-logistic"
data_dir <- "data"
fit_logistic_boot_dir <- "fit-logistic-boot"

# Functions ===================================================================

invlogit <- function(x) 1 - 1 / (1 + exp(x))

# Reads a set of results from path
read_res <- function(path) {
  read_csv(path, col_types = cols()) %>%
    filter(virus != "H1N1seas") %>%
    mutate(
      population = factor(
        population, 
        levels = c("general", "exposed"), 
        labels = c("General", "Exposed"))
    )
}

fit_lr_one <- function(data, formula) {
  fit <- glm(formula, data, family = binomial(link = "logit"))
  attr(fit, "virus") <- unique(data$virus)
  attr(fit, "population") <- unique(data$population)
  fit
}

predict_lr_one <- function(fit_one, logHIs) {
  predict(fit_one, data.frame(logHImid = logHIs), se.fit = TRUE) %>% 
    as_tibble() %>%
    mutate(
      fit_high = invlogit(fit + qnorm(0.975) * se.fit),
      fit_low = invlogit(fit - qnorm(0.975) * se.fit),
      fit = invlogit(fit),
      logHImid = logHIs,
      population = attr(fit_one, "population"),
      virus = attr(fit_one, "virus")
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
    geom_ribbon(alpha = 0.5),
    geom_line()
  )
}

prot_curve_fun <- function(outsum, ylabl = "Protection") {
  outsum %>%
    ggplot(
      aes(logHImid, prot_med, ymin = prot_low, ymax = prot_high)
    ) +
    ylab(ylabl) +
    coord_cartesian(xlim = c(log(5), log(1280)), ylim = c(0, 1)) +
    facet_grid(population ~ virus) +
    common_plot_els()
}

inf_curve_fun <- function(outsum, data) {
  outsum %>%
    ggplot(aes(logHImid, fit, ymin = fit_low, ymax = fit_high)) +
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
      )
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

save_plot <- function(pl, name) {
  plpath <- file.path(fit_logistic_dir, paste0(name, ".pdf"))
  ggsave_dark(
    plot = pl, filename = plpath, dark = FALSE,
    width = 10, height = 10, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

# Hanam data
han <- read_res(file.path(data_dir, "hanam-HI-exp.csv")) %>%
  bind_rows(read_res(file.path(data_dir, "hanam-HI-gen.csv")))

# Hanam summary data
han_summ <- read_res(file.path(data_dir, "hanam-HI-summ.csv"))

# Bootstrap samples
han_lr_boot <- read_res(file.path(fit_logistic_boot_dir, "samples.csv"))

# Log HI's for which to calculate infection/protection probabilities
logHIs <- seq(0, 7.5, length.out = 101)

# Regular fit
fits <- han %>%
  group_split(virus, population) %>%
  map(fit_lr_one, status_bin ~ logHImid)

# Predictions from regular fits
preds <- fits %>% map_dfr(predict_lr_one, logHIs)

# Infection curves
inf_curve <- inf_curve_fun(preds, han_summ)
save_plot(inf_curve, "inf")

# Protection curves abolute
prot_abs <- preds %>%
  mutate(prot_med = 1 - fit, prot_low = 1 - fit_low, prot_high = 1 - fit_high)
prot_curve <- prot_curve_fun(prot_abs)
save_plot(prot_curve, "prot-abs")

# Protection curve relative
prot_rel <- han_lr_boot %>%
  mutate(term = recode(term, "(Intercept)" = "b0", "logHImid" = "bt")) %>%
  select(virus, population, term, estimate, ind) %>%
  pivot_wider(names_from = "term", values_from = estimate) %>%
  slice(rep(1:n(), each = length(logHIs))) %>%
  mutate(
    logHImid = rep_len(logHIs, length.out = n()),
    fit = invlogit(b0 + bt * logHImid),
    fit5 = invlogit(b0 + bt * log(5)),
    fit_rel = fit / fit5
  ) %>%
  group_by(virus, population, logHImid) %>%
  summarise(
    prot_med = 1 - median(fit_rel),
    prot_high = 1 - quantile(fit_rel, 0.975),
    prot_low = 1 - quantile(fit_rel, 0.025)
  ) %>%
  ungroup()

prot_rel_plot <- prot_curve_fun(prot_rel, "Protection relative to 5")
save_plot(prot_rel_plot, "prot-rel")
