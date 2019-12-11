# Fit the scaled logit model to hanam data
# Arseniy Khvorov
# Created 2019/12/04
# Last edit 2019/12/04

library(tidyverse)
library(sclr)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories to be used later
fit_sclr_dir <- "fit-sclr"
data_dir <- "data"
fit_sclr_boot_dir <- "fit-sclr-boot"

# Functions ===================================================================

inf_curve <- function(dat, count_data) {
  dat %>%
    ggplot(aes(logHImid, inf_point)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.spacing = unit(0, "null"),
      strip.background = element_rect(fill = NA),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    xlab("HI titre") +
    ylab("Infection probability") +
    coord_cartesian(xlim = c(log(5), log(1280)), ylim = c(0, 0.3)) +
    scale_x_continuous(
      labels = 5 * 2^(0:8), breaks = log(5 * 2^(0:8)),
      minor_breaks = log(7.5 * 2^(0:8))
    ) +
    scale_y_continuous(
      labels = scales::percent_format(1), breaks = seq(0, 1, 0.1)
    ) +
    facet_wrap(~virus, nrow = 1) +
    geom_ribbon(aes(ymin = inf_low, ymax = inf_high), alpha = 0.5) +
    geom_line() +
    geom_point(
      data = count_data, mapping = aes(logHImid, prop_inf), inherit.aes = FALSE,
      shape = 18, color = "wheat"
    ) +
    geom_text(
      data = count_data, mapping = aes(logHImid, prop_inf, label = n_tot), 
      inherit.aes = FALSE, 
      nudge_x = ifelse(near(count_data$logHImid, log(5)), 0.55, 0), 
      nudge_y = ifelse(near(count_data$logHImid, log(5)), 0, 0.016),
      color = "wheat"
    )
}

prot_curve <- function(dat) {
  dat %>%
    ggplot(aes(logHImid, prot_point)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.spacing = unit(0, "null"),
      strip.background = element_rect(fill = NA),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    xlab("HI titre") +
    ylab("Protection") +
    coord_cartesian(xlim = c(log(5), log(1280))) +
    scale_x_continuous(
      labels = 5 * 2^(0:8), breaks = log(5 * 2^(0:8)),
      minor_breaks = log(7.5 * 2^(0:8))
    ) +
    scale_y_continuous(
      labels = scales::percent_format(1), breaks = seq(0, 1, 0.1)
    ) +
    facet_wrap(~virus, nrow = 1) +
    geom_ribbon(aes(ymin = prot_l, ymax = prot_u), alpha = 0.5) +
    geom_line()
}

save_plot <- function(pl, name, dark) {
  suff <- if_else(dark, "_dark", "_light")
  plpath <- file.path(fit_sclr_dir, paste0(name, suff, ".pdf"))
  ggsave_dark(
    plot = pl, filename = plpath, dark = dark,
    width = 15, height = 7.5, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

# Hanam data
han <- read_csv(file.path(data_dir, "hanam.csv"))

# Prepare
han_prep <- han %>%
  filter(!is.na(status) | !is.na(preHI), virus != "H1N1seas") %>%
  mutate(
    logHI = log(preHI),
    logHImid = case_when(
      preHI == 5 ~ log(5),
      preHI == 1280 ~ log(2560),
      TRUE ~ logHI + log(2) / 2
    ),
    status_bin = if_else(status == "Not infected", 0 ,1)
  ) %>%
  select(status, status_bin, preHI, logHI, logHImid, virus)

# Fit to viruses separately
fit_normal <- han_prep %>%
  group_split(virus) %>%
  map(~ sclr(status_bin ~ logHImid, .x))
names(fit_normal) <- group_keys(han_prep, virus)$virus

logHIs <- tibble(logHImid = seq(1, 7.5, length.out = 101))

# Protection curves
pred_dat <- map_dfr(fit_normal, ~ predict(.x, logHIs), .id = "virus")
pred_dat_plot <- prot_curve(pred_dat)
walk(c(TRUE, FALSE), ~ save_plot(pred_dat_plot, "protection-one", .x))

# Infection curves (need bootstrap data for this)
boot_samples <- read_csv(file.path(fit_sclr_boot_dir, "samples.csv"))

boot_samples_wide <- boot_samples %>%
  select(virus, term, estimate, ind) %>%
  pivot_wider(names_from = "term", values_from = "estimate")

inf_pred <- boot_samples_wide %>%
  slice(rep(1:n(), nrow(logHIs))) %>%
  bind_cols(logHIs %>% slice(rep(1:n(), each = nrow(boot_samples_wide)))) %>%
  mutate(
    lambda = 1 - 1 / (1 + exp(theta)),
    inf_prob = lambda / (1 + exp(beta_0 + beta_logHImid * logHImid)),
    prot_prob = 1 - 1 / (1 + exp(beta_0 + beta_logHImid * logHImid))
  ) %>%
  group_by(logHImid, virus) %>%
  summarise(
    inf_point = median(inf_prob),
    inf_low = quantile(inf_prob, 0.025),
    inf_high = quantile(inf_prob, 0.975),
    prot_point = median(prot_prob),
    prot_l = quantile(prot_prob, 0.025),
    prot_u = quantile(prot_prob, 0.975)
  ) %>%
  ungroup()

han_binned <- han_prep %>%
  filter(!is.na(status), !is.na(logHImid)) %>%
  group_by(logHImid, virus) %>%
  summarise(prop_inf = sum(status != "Not infected") / n(), n_tot = n())

inf_boot <- inf_curve(inf_pred, han_binned)
prot_boot <- prot_curve(inf_pred)
walk(c(TRUE, FALSE), ~ save_plot(inf_boot, "infection-boot", .x))
walk(c(TRUE, FALSE), ~ save_plot(prot_boot, "protection-boot", .x))
