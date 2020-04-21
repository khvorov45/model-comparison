# Plots of kiddyvax data
# Arseniy Khvorov
# Created 2020-01-23
# Last edit 2020-01-24

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories to be used later
data_dir <- "data"
data_plot_dir <- "data-plot"

# Functions ===================================================================

read_kiddyvax <- function(name) {
  read_csv(
    file.path(data_dir, paste0(name, ".csv")),
    col_types = cols(
      id = col_integer(),
      virus = col_character(),
      status = col_integer(),
      hi = col_integer()
    )
  )
}

read_swab <- function(name) {
  read_csv(
    file.path(data_dir, paste0(name, ".csv")),
    col_types = cols(
      id = col_integer(),
      swab_date = col_date("%Y-%m-%d"),
      virus = col_character(),
      swab_result = col_integer(),
      start_date = col_date("%Y-%m-%d"),
      end_date = col_date("%Y-%m-%d")
    )
  )
}

lbl_virus <- function(dat) {
  dat %>%
    mutate(
      virus_lbl = factor(
        virus,
        levels = c("bvic", "byam", "h1pdm", "h1seas", "h3"),
        labels = c("B Vic", "B Yam", "A H1pdm", "A H1seas", "A H3")
      )
    )
}

lbl_status <- function(dat) {
  dat %>%
    mutate(
      status_lbl = factor(
        status,
        levels = c(0, 1), labels = c("Not infected", "Infected")
      )
    )
}

cmn_thm_els <- function() {
  theme(
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )
}

inf_plot <- function(kiddyvax) {
  counts <- kiddyvax %>%
    filter(!is.na(status)) %>%
    count(virus, status) %>%
    lbl_virus() %>%
    lbl_status()
  kiddyvax %>%
    filter(!is.na(status), !is.na(hi)) %>%
    lbl_virus() %>%
    lbl_status() %>%
    ggplot(aes(status_lbl, hi)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.grid.minor = element_blank()
    ) +
    cmn_thm_els() +
    ylab("HI titre") +
    xlab("Status") +
    scale_y_log10(breaks = 5 * 2^(0:12)) +
    facet_wrap(~virus_lbl, nrow = 1) +
    geom_jitter(alpha = 0.5, shape = 18) +
    geom_text(
      data = counts, mapping = aes(status_lbl, 10240, label = n),
      inherit.aes = FALSE
    )
}

swab_plot <- function(swab) {
  swab %>%
    lbl_virus() %>%
    group_by(start_date, id) %>%
    mutate(
      id2 = group_indices(),
      swab_result_lbl = factor(
        swab_result,
        levels = c(1, 0, NA),
        labels = c("Positive", "Negative", "Missing"),
        exclude = NULL
      )
    ) %>%
    ungroup() %>%
    ggplot(aes(start_date, id2)) +
    dark_theme_bw(verbose = FALSE) +
    cmn_thm_els() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    scale_color_discrete("Swab result") +
    scale_x_date(expand = c(0.01, 0.01), minor_breaks = "month") +
    scale_y_continuous(
      expand = c(0.01, 0.01), breaks = seq_len(length(unique(swab$id)))
    ) +
    facet_wrap(~virus_lbl, nrow = 1) +
    geom_segment(
      aes(x = start_date, xend = end_date, y = id2, yend = id2),
      lwd = 0.1, col = "gray50"
    ) +
    geom_point(aes(swab_date, id2, color = swab_result_lbl), shape = 17)
}

# Plots the antibody titre and infection status
plot_counts <- function(dat) {
  xbreaks <- c(5 * 2^(0:12))
  ybreaks <- seq(0, 1, 0.1)
  dat %>%
    ggplot(aes(log(hi), inf_prop)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor.y = element_blank()
    ) +
    xlab("HI titre") +
    ylab("Infected proportion") +
    scale_x_continuous(
      breaks = log(xbreaks), labels = xbreaks,
      minor_breaks = log((xbreaks + xbreaks * 2) / 2)
    ) +
    scale_y_continuous(labels = scales::percent_format(1)) +
    geom_point(
      data = dat, mapping = aes(log(hi), inf_prop), inherit.aes = FALSE,
      shape = 18
    ) +
    ggrepel::geom_text_repel(
      data = dat, mapping = aes(log(hi), inf_prop, label = ntot),
      inherit.aes = FALSE, color = "gray50"
    )
}

save_plot <- function(plot, name, width, height) {
  ggsave_dark(
    plot,
    filename = file.path(data_plot_dir, paste0(name, ".pdf")),
    dark = FALSE,
    width = width, height = height, units = "cm"
  )
}

# Script ======================================================================

kiddyvax <- read_kiddyvax("kiddyvaxmain")
swab <- read_swab("kiddyvaxmain-swab")
summ <- read_csv(file.path(data_dir, "kiddyvaxmain-summ.csv"))

summ_plot <- plot_counts(filter(summ, virus == "bvic"))
save_plot(summ_plot, "kiddyvax-main-summ", 7.5, 7.5)

kiddyvax_inf <- inf_plot(kiddyvax)
# save_plot(kiddyvax_inf, "kiddyvax-main-titre", 15, 7.5)

kiddyvax_swab <- swab_plot(swab)
# save_plot(last_plot(), "kiddyvax-swab", 50, 65)
