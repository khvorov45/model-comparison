# Script to create sample Cox time data and plot it

library(tidyverse)

# Directories to be used later
curve_cox_dir <- here::here("curve-cox")

# Functions ===================================================================

pl_assumption <- function(dat, top) {
  brks <- unique(dat$subject)
  dat %>%
    ggplot(aes(x = subject, y = time_end, ymin = time_start, ymax = time_end)) +
    ggdark::dark_theme_bw(verbose = FALSE) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "lines"),
      panel.grid.minor.y = element_blank()
    ) +
    xlab("Subject") +
    ylab("Time, days") +
    scale_x_continuous(breaks = brks) +
    scale_shape_manual(values = c(18, 5)) +
    coord_flip(ylim = c(0, 300), xlim = c(0, top)) +
    geom_rect(
      aes(xmin = -10, xmax = 30, ymin = 100, ymax = 200),
      alpha = 0.15, fill = "gray50", color = NA
    ) +
    geom_linerange() +
    geom_point(aes(shape = infected))
}

save_plot <- function(pl, pl_ind) {
  ggdark::ggsave_dark(
    plot = pl, dark = FALSE,
    filename = file.path(
      curve_cox_dir, paste0("timeplot_", pl_ind - 1, ".pdf")
    ),
    width = 10, height = 5, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

dat <- tibble(
  subject = c("1", "2"),
  time_start = c(0, 0),
  time_end = c(100, 150)
)

pl0 <- ggplot(
  dat, aes(x = subject, y = time_end, ymin = time_start, ymax = time_end)
) +
  ggdark::dark_theme_bw(verbose = FALSE) +
  xlab("Subject") +
  ylab("Time, days") +
  coord_flip()

pl1 <- pl0 +
  geom_pointrange()

dat_real <- tibble(
  subject = c(rep("1", 5), "2"),
  time_start = c(0, 10, 35, 75, 95, 145),
  time_end = c(5, 15, 40, 80, 100, 150)
)

pl3 <- pl0 +
  geom_linerange(data = dat_real) +
  geom_point(data = dat)

dat_act <- tibble(
  subject = c(1, 2),
  infected = c("Infected", "Not infected"),
  time_start = c(100, 100),
  time_end = c(150, 200)
)

pl_act <- pl_assumption(dat_act, 3)

dat_act2 <- bind_rows(
  dat_act,
  tibble(
    subject = c(3, 4),
    infected = c("Infected", "Not infected"),
    time_start = c(50, 25),
    time_end = c(150, 200)
  )
)

pl_act2 <- pl_assumption(dat_act2, 5)

pls <- list(pl0, pl1, pl3, pl_act, pl_act2)

iwalk(pls, save_plot)
