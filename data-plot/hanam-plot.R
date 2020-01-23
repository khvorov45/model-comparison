# Relevant graphs of hanam data
# Arseniy Khvorov
# Created 2019/09/26
# Last edit 2019/11/22

library(tidyverse)
library(ggdark) # devtools::install("khvorov45/ggdark")

# Folders to be used later
data_dir <- "data"
graph_data_dir <- "data-plot"

# Functions ===================================================================

# Reads one csv file
read_one <- function(nme) {
  read_csv(file.path(data_dir, paste0(nme, ".csv")), col_types = cols())
}

# Plots the antibody titre and infection status
plot_counts <- function(dat, facet_type = "vir") {
  xbreaks <- c(5 * 2^(0:8))
  ybreaks <- seq(0, 1, 0.1)
  facets <- list(
    "vir" = list(
      facet_wrap(vars(virus), nrow = 1), 
      coord_cartesian(xlim = c(log(5), log(1280)), ylim = c(0, 0.3)),
      geom_text(
        data = dat, mapping = aes(log(preHI), inf_prop, label = ntot), 
        inherit.aes = FALSE, 
        nudge_x = ifelse(near(dat$preHI, 5), 0.55, 0), 
        nudge_y = ifelse(near(dat$preHI, 5), 0, 0.016)
      ),
      scale_y_continuous(breaks = ybreaks, labels = scales::percent_format(1))
    ),
    "virpop" = list(
      facet_grid(population ~ virus, scales = "free_y"),
      ggrepel::geom_text_repel(
        data = dat, mapping = aes(log(preHI), inf_prop, label = ntot), 
        inherit.aes = FALSE, color = "gray50"
      ),
      scale_y_continuous(labels = scales::percent_format(1)),
      theme(panel.spacing.y = unit(0.25, "lines"))
    )
  )
  dat %>%
    ggplot(aes(log(preHI), inf_prop)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_rect(fill = NA),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.grid.minor.y = element_blank()
    ) +
    xlab("HI titre") + ylab("Infected proportion") +
    scale_x_continuous(
      breaks = log(xbreaks), labels = xbreaks,
      minor_breaks = log((xbreaks + xbreaks * 2) / 2)
    ) +
    facets[[facet_type]] +
    geom_point(
      data = dat, mapping = aes(log(preHI), inf_prop), inherit.aes = FALSE,
      shape = 18
    )
}

plot_scatter <- function(dat, y_breaks = 5 * 2^(0:8), facet_type = "vir") {
  facets <- list(
    "vir" = facet_wrap(vars(virus), nrow = 1),
    "virpop" = facet_grid(population ~ virus)
  )
  facet_type <- rlang::arg_match(facet_type, names(facets))
  dat %>%
    ggplot(aes(status, preHI)) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      strip.background = element_rect(fill = NA),
      panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines")
    ) +
    xlab("Status") +
    scale_y_log10("HI titre", breaks = y_breaks, labels = y_breaks) +
    coord_cartesian(ylim = c(4, 1280)) +
    facets[[facet_type]] +
    geom_point(alpha = 0.4, position = "jitter", shape = 16)
}

# Saves the plot with the name into the folder
save_plot <- function(pl, dark, name, folder, height = 7.5) {
  postfix <- if_else(dark, "dark", "light")
  plname <- file.path(folder, paste0(name, "_", postfix, ".pdf"))
  ggsave_dark(
    plname, pl, dark,
    width = 15, height = height, units = "cm", device = "pdf"
  )
}

# Saves dark and light versions of the plot
save_plot_dl <- function(pl, name, folder, height = 7.5) {
  walk(c(TRUE, FALSE), ~ save_plot(pl, .x, name, folder, height))
}

# Script ======================================================================

han_dat <- map_dfr(c("hanam-HI-gen", "hanam-HI-exp"), read_one) %>%
  filter(!is.na(status), !is.na(preHI), virus != "H1N1seas") %>%
  mutate(
    status_bin = if_else(status_bin == 0, "Not infected", "Infected"),
    population = factor(
      population, levels = c("general", "exposed"), 
      labels = c("General", "Exposed")
    )
  )

han_cnts <- read_one("hanam-HI-summ") %>%
  filter(virus != "H1N1seas") %>%
  mutate(
    population = factor(
      population, levels = c("general", "exposed"), 
      labels = c("General", "Exposed")
    )
  )

pls_cnts <- plot_counts(han_cnts, "virpop")
save_plot_dl(pls_cnts, "hanam_counts", graph_data_dir, 10)

pls_scat_gen <- plot_scatter(filter(han_limited, population == "General"))
save_plot(pls_scat_gen, TRUE, "hanam_scatter_general", graph_data_dir)

pls_scat <- plot_scatter(han_limited, facet_type = "virpop")
save_plot_dl(pls_scat, "hanam_scatter", graph_data_dir)

pls_cnts_gen <- plot_counts(filter(han_cnts, population == "General"), "vir")
save_plot(pls_cnts_gen, TRUE, "hanam_counts_general", graph_data_dir)
