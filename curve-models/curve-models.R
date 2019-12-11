# Script for showing probability curves
# Arseniy Khvorov
# Created 2019/09/06
# Last edit 2019/12/05

library(tidyverse)
library(ggdark) # devtools:: install_github("khvorov45/ggdark")

# Directories to be used later
curve_moels_dir <- "curve-models"

# Functions ===================================================================

logistic_curve_pp <- function(x, top = 1, b0dbx = log(40), bx = -2.5) {
  (top) / (1 + exp(bx * (b0dbx - x)))
}

exponential_curve <- function(x, bx = -2, top = 0.8) {
  top * exp(bx * (x - log(5)))
}

add_labels <- function(pl, model_lbl, top_y, top_lbl, col, 
                       dashend = log(160)) {
  pl +
    geom_text(
      x = log(2.5), y = top_y + 0.01, label = model_lbl,
      hjust = 0, vjust = 0, col = col
    ) +
    geom_segment(
      x = log(5), y = top_y, xend = dashend, yend = top_y,
      lty = "1313", col = col
    ) +
    geom_text(
      x = dashend + 0.1, y = top_y, label = top_lbl,
      hjust = 0, vjust = 0.5, col = col
    )
}

save_plot <- function(pl, ind) {
  ggsave(
    pl,
    filename = file.path(
      curve_moels_dir,
      paste0("curve_", ind - 1, "_dark.pdf")
    ),
    device = "pdf", width = 10, height = 7.5, units = "cm"
  )
}

# Script ======================================================================

dat <- tibble(x = seq(1, 7, length.out = 101))
breaks <- c(1, 5 * 2^(0:8))
logistic_col <- "mediumorchid"
sclr_col <- "violetred"
sclr_top <- 0.6
cox_col <- "slateblue1"
cox_top <- 0.8

curve_0 <- ggplot(dat, aes(x, seq(0, 1, length.out = nrow(dat)))) +
  dark_theme_bw(verbose = FALSE) +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous("Infection probability", breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(
    "Antibody titre",
    breaks = log(breaks), labels = breaks
  ) +
  coord_cartesian(ylim = c(0, 1.05), xlim = c(1, 7))
curve_0

curve_1 <- curve_0 +
  stat_function(fun = logistic_curve_pp, col = logistic_col)
curve_1

curve_2 <- curve_1 %>%
  add_labels("Logistic regression", 1, "Fixed at 1", logistic_col)
curve_2

curve_3 <- curve_2 +
  stat_function(
    fun = logistic_curve_pp, args = list(top = sclr_top), col = sclr_col
  )
curve_3

curve_4 <- curve_3 %>%
  add_labels("Scaled logistic regression", sclr_top, "Estimated", sclr_col)
curve_4

curve_5 <- curve_4 +
  stat_function(
    fun = exponential_curve, args = list(bx = -0.5, top = cox_top),
    col = cox_col, xlim = c(log(5), log(1280))
  )
curve_5

curve_6 <- curve_5 %>%
  add_labels("Cox PH", cox_top, "Not estimated", cox_col)
curve_6

curve_sclr <- curve_0 +
  stat_function(
    fun = logistic_curve_pp, args = list(top = sclr_top), col = sclr_col
  )
curve_sclr <- add_labels(
  curve_sclr,
  "Scaled logistic regression", sclr_top, "Estimated", sclr_col
)
curve_sclr

curves_all <- list(
  curve_0, curve_1, curve_2, curve_3, curve_4, curve_5, curve_6, curve_sclr
)

iwalk(curves_all, save_plot)
