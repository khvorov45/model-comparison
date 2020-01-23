# Plots of kiddyvax data
# Arseniy Khvorov
# Created 2020-01-23
# Last edit 2020-01-23

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

# Script ======================================================================

kiddyvax <- read_kiddyvax("kiddyvax-main")

kiddyvax %>%
  filter(!is.na(status), !is.na(hi)) %>%
  mutate(
    status_lbl = factor(
      status, levels = c(0, 1), labels = c("Not infected", "Infected")
    ),
    virus_lbl = factor(
      virus, levels = c("bvic", "byam", "h1pdm", "h1seas", "h3"),
      labels = c("B Vic", "B Yam", "A H1pdm", "A H1seas", "A H3")
    )
  ) %>%
  ggplot(aes(status_lbl, hi)) +
  dark_theme_bw(verbose = FALSE) +
  theme(
    panel.spacing = unit(0, "lines"),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  ylab("HI titre") +
  xlab("Status") +
  scale_y_log10(breaks = 5 * 2^(0:12)) +
  facet_wrap(~virus_lbl, nrow = 1) +
  geom_jitter(alpha = 0.5, shape = 18)

ggsave_dark(
  filename = file.path(data_plot_dir, "kiddyvax-main-titre.pdf"),
  dark = FALSE,
  width = 15, height = 7.5, units = "cm"
)


