# Relevant graphs of data

library(tidyverse)

# Folders to be used later
data_dir <- here::here("data")
data_plot_dir <- here::here("data-plot")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

# Plots the antibody titre and infection status
plot_counts <- function(dat, facet_type = "vir") {
  xbreaks <- c(5 * 2^(0:8))
  ybreaks <- seq(0, 1, 0.1)
  facets <- list(
    "vir" = list(
      facet_wrap(vars(virus), nrow = 1),
      coord_cartesian(xlim = c(log(5), log(1280)), ylim = c(0, 0.3)),
      scale_y_continuous(breaks = ybreaks, labels = scales::percent_format(1))
    ),
    "virpop" = list(
      facet_grid(population ~ virus, scales = "free_y"),
      scale_y_continuous(labels = scales::percent_format(1)),
      theme(panel.spacing.y = unit(0.25, "lines"))
    ),
    "pop" = list(
      facet_wrap(vars(population), nrow = 1, scales = "free_y"),
      coord_cartesian(xlim = c(log(5), log(1280))),
      scale_y_continuous(labels = scales::percent_format(1)),
      theme(panel.spacing.x = unit(1, "lines"))
    ),
    "none" = list()
  )
  dat %>%
    ggplot(aes(log(prehi), inf_prop)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
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
    facets[[facet_type]] +
    geom_point(
      data = dat, mapping = aes(log(prehi), inf_prop), inherit.aes = FALSE,
      shape = 18
    ) +
    ggrepel::geom_text_repel(
      data = dat, mapping = aes(log(prehi), inf_prop, label = ntot),
      inherit.aes = FALSE, color = "gray50"
    )
}

plot_scatter <- function(dat, y_breaks = 5 * 2^(0:8), facet_type = "vir") {
  facets <- list(
    "vir" = facet_wrap(vars(virus), nrow = 1),
    "virpop" = facet_grid(population ~ virus)
  )
  facet_type <- rlang::arg_match(facet_type, names(facets))
  dat %>%
    ggplot(aes(status_bin_lbl, prehi)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines")
    ) +
    xlab("Status") +
    scale_y_log10("HI titre", breaks = y_breaks, labels = y_breaks) +
    coord_cartesian(ylim = c(4, 1280)) +
    facets[[facet_type]] +
    geom_point(alpha = 0.4, position = "jitter", shape = 16)
}

swab_plot <- function(swab) {
  swab %>%
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
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1)
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

# Saves the plot with the name into the folder
save_plot <- function(pl, name, height = 7.5) {
  plname <- file.path(data_plot_dir, paste0(name, ".pdf"))
  ggdark::ggsave_dark(
    plname, pl,
    width = 15, height = height, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

han_dat <- map_dfr(c("hanam-hi-gen", "hanam-hi-exp"), read_data) %>%
  filter(!is.na(status), !is.na(prehi), virus != "H1N1seas")
han_cnts <- read_data("hanam-hi-summ") %>%
  filter(virus != "H1N1seas")

kiddyvax <- read_data("kiddyvaxmain")
kv_swab <- read_data("kiddyvaxmain-swab")
kv_summ <- read_data("kiddyvaxmain-summ")

all_plots <- list(
  "hanam-hi-summ" = plot_counts(han_cnts, "virpop"),
  "hanam-hi-summ-h3n2" = plot_counts(filter(han_cnts, virus == "H3N2"), "pop"),
  "hanam-hi-gen-scatter" = plot_scatter(
    filter(han_dat, population == "General")
  ),
  "hanam-hi-scatter" = plot_scatter(han_dat, facet_type = "virpop"),
  "hanam-hi-summ-gen" = plot_counts(
    filter(han_cnts, population == "General"), "vir"
  ),
  "kiddyvax-main-summ" = plot_counts(
    filter(kv_summ, virus == "bvic") %>% mutate(prehi = hi),
    "none"
  ),
  "kiddyvax-main-swab" = swab_plot(kv_swab)
)

iwalk(all_plots, save_plot)
