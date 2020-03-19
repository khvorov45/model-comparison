# Graphing fitting results
# Arseniy Khvorov
# Created 2019/09/20
# Last edit 2019/12/05

library(tidyverse)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(latex2exp)

# Directories to be used later
logistic_plot_dir <- "logistic-plot"
logistic_summ_dir <- "logistic-summary"

# Functions ===================================================================

sclr_curve <- function(x, l, b0, b1) l / (1 + exp(b0 + b1 * x))

calc_fit_prob <- function(data) {
  l <- data$est_mean[data$term == "theta"]
  if (identical(l, numeric(0))) {
    l <- 1
  } else {
    l <- 1 - 1 / (1 + exp(l))
  }
  b0 <- data$est_mean[data$term == "beta_0"]
  b1 <- data$est_mean[data$term == "beta_logTitre"]
  tibble(.rows = 101) %>%
    mutate(
      x = seq(0, 7, length.out = 101),
      fit_prob = sclr_curve(x, l, b0, b1)
    )
}

ex_plot_1 <- function(data, x_name = "x", y_name = "fit_prob") {
  xbreaks <- c(1, 2, 5 * 2^(0:8))
  colshlab <- as_labeller(
    c("logistic" = "Logistic", "scaled_logit" = "Scaled logit")
  )
  data %>%
    ggplot(aes(!!rlang::sym(x_name), !!rlang::sym(y_name))) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      legend.position = "bottom",
      legend.box.spacing = unit(0, "lines")
    ) +
    xlab("Titre") +
    ylab("Infection probability") +
    labs(lty = "Model", col = "Model", fill = "Model") +
    scale_color_discrete(labels = colshlab) +
    scale_linetype_manual(labels = colshlab, values = c("31", "11")) +
    scale_fill_discrete(labels = colshlab) +
    scale_x_continuous(
      labels = xbreaks, breaks = log(xbreaks),
      minor_breaks = log(c(1.5, 3.5, 7.5 * 2^(0:8)))
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_cartesian(ylim = c(0, 1)) +
    stat_function(fun = sclr_curve, args = list(0.5, -5, 1.5)) +
    geom_line(aes(lty = model, col = model))
}

vary_par_conv <- function(x_name, xlims, x_lab, data) {
  model_labs <- as_labeller(
    c("logistic" = "Logistic", "scaled_logit" = "Scaled logit")
  )
  data %>%
    filter(par_varied == x_name) %>%
    ggplot(
      aes(!!sym(x_name), converged, shape = model, lty = model, col = model)
    ) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      legend.margin = margin(0, 0, 0, 0, "null")
    ) +
    labs(shape = "Model", lty = "Model", col = "Model") +
    scale_x_continuous(x_lab, breaks = unique(data[[x_name]])) +
    scale_y_continuous(
      "Proportion converged",
      labels = scales::percent_format(1)
    ) +
    scale_linetype_manual(labels = model_labs, values = c("31", "11")) +
    scale_color_discrete(labels = model_labs) +
    scale_shape_manual(labels = model_labs, values = c(17, 18)) +
    coord_cartesian(xlim = xlims) +
    geom_line() +
    geom_point()
}

vary_par_se <- function(x_name, xlims, x_lab, data, y_name = "se_mean",
                        y_lab = "Expected standard error") {
  term_labeller <- function(term_names) {
    term_names[, 1] <- recode(
      term_names[, 1],
      "beta_0" = TeX("$\\beta_0$"),
      "beta_logTitre" = TeX("$\\beta_X$")
    )
    label_parsed(term_names)
  }
  model_labs <- as_labeller(
    c("logistic" = "Logistic", "scaled_logit" = "Scaled logit")
  )
  data %>%
    filter(par_varied == x_name, term != "theta") %>%
    ggplot(
      aes(!!sym(x_name), !!sym(y_name), shape = model, lty = model, col = model)
    ) +
    dark_theme_bw(verbose = FALSE) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box.spacing = unit(0, "null"),
      panel.spacing = unit(0, "null"),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.5, "lines")
    ) +
    labs(shape = "Model", lty = "Model", col = "Model") +
    scale_x_continuous(x_lab, breaks = unique(data[[x_name]])) +
    scale_y_continuous(y_lab) +
    scale_linetype_manual(labels = model_labs, values = c("31", "11")) +
    scale_shape_manual(labels = model_labs, values = c(17, 18)) +
    scale_color_discrete(labels = model_labs) +
    coord_cartesian(xlim = xlims) +
    facet_wrap(
      ~term,
      ncol = 1, labeller = term_labeller, strip.position = "right",
      scales = "free_y"
    ) +
    geom_line() +
    geom_point()
}

plot_pres_series <- function(b0, b1) {
  xbreaks <- c(1, 2, 5 * 2^(0:8))
  sclr_curve_rel <- function(x, l, b0, b1) {
    sclr_curve(x, l, b0, b1) / sclr_curve(log(5), l, b0, b1)
  }
  sclr_curve_rel_1m <- function(x, l, b0, b1) {
    1 - sclr_curve_rel(x, l, b0, b1)
  }
  pl0 <- tibble(
    x = seq(0, 7, length.out = 101), y = seq(0, 1, length.out = 101)
  ) %>%
    ggplot(aes(x, y)) +
    dark_theme_bw(verbose = FALSE) +
    xlab("Titre") +
    ylab("Infection probability") +
    scale_x_continuous(
      labels = xbreaks, breaks = log(xbreaks),
      minor_breaks = log(c(1.5, 3.5, 7.5 * 2^(0:8)))
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    coord_cartesian(ylim = c(0, 1))
  pl1 <- pl0 +
    stat_function(fun = sclr_curve, args = list(0.5, -5, 1.5), lty = "1111")
  pl2 <- pl1 +
    stat_function(fun = sclr_curve, args = list(1, b0, b1))
  pl3 <- pl0 +
    stat_function(fun = sclr_curve, args = list(1, -5, 1.5), lty = "1111") +
    ylab("Relative infection probability")
  pl4 <- pl3 +
    stat_function(fun = sclr_curve, args = list(1, b0, b1))
  pl5 <- pl0 +
    stat_function(fun = sclr_curve, args = list(1, 5, -1.5), lty = "1111") +
    stat_function(fun = sclr_curve, args = list(1, -b0, -b1)) +
    ylab("Protection")
  pl6 <- pl0 +
    stat_function(
      fun = sclr_curve_rel, args = list(0.5, -5, 1.5), lty = "1111"
    ) +
    ylab("Infection probability (relative to log(5))")
  pl7 <- pl6 +
    stat_function(fun = sclr_curve_rel, args = list(1, b0, b1))
  pl8 <- pl0 +
    stat_function(
      fun = sclr_curve_rel_1m, args = list(0.5, -5, 1.5), lty = "1111"
    ) +
    stat_function(fun = sclr_curve_rel_1m, args = list(1, b0, b1)) +
    ylab("Protection")
  list(pl0, pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8)
}

save_plot <- function(pl, name, dark, height = 7.5) {
  ggsave_dark(
    file.path(logistic_plot_dir, paste0(name, ".pdf")), pl,
    dark = dark, width = 10, height = height, units = "cm", device = "pdf"
  )
}

# Script ======================================================================

summ <- read_csv(file.path(logistic_summ_dir, "summary-10000sims.csv"))

# Example of poor logistic fit
summ_ex <- summ %>%
  filter(theta == 0, nsam == 1e4) %>%
  group_by(model) %>%
  group_modify(~ calc_fit_prob(.x))
lrex1 <- ex_plot_1(summ_ex)
# save_plot(lrex1, "lrex", FALSE)

# For presentation
b_lr <- summ %>% filter(model == "logistic", theta == 0, nsam == 1e4)
b0 <- b_lr$est_mean[b_lr$term == "beta_0"]
b1 <- b_lr$est_mean[b_lr$term == "beta_logTitre"]
pres_series <- plot_pres_series(b0, b1)
# iwalk(pres_series, ~ save_plot(.x, paste0("lrex_pres", .y), TRUE))

# Convergence with sample size
vary_nsam <- vary_par_conv("nsam", c(0, 600), "Sample size", summ)
# save_plot(vary_nsam, "vary_nsam", FALSE, 6)

# Convergence at varying lambdas
vary_theta <- vary_par_conv(
  "lambda", c(0, 1), "Lambda",
  summ %>% mutate(lambda = 1 - 1 / (1 + exp(theta)))
)
# save_plot(vary_theta, "vary_lambda", FALSE)

# SE with sample size
nsam_se_plot <- vary_par_se("nsam", c(0, 800), "Sample size", summ)
# save_plot(nsam_se_plot, "vary_nsam_se", FALSE)

# Means with sample size
nsam_mean_plot <- vary_par_se(
  "nsam", c(0, 800), "Sample size", summ, "est_mean", "Estimate mean"
) +
  geom_hline(aes(yintercept = true_value))
# save_plot(nsam_mean_plot, "vary_nsam_mean", FALSE)

# Poor logistic fit with SEs
preds <- read_csv(file.path(logistic_summ_dir, "preds-10000sims.csv"))

preds_plot <- preds %>%
  filter(nsam == 500, theta_true == 0) %>%
  ex_plot_1("logTitre", "50%") +
  geom_ribbon(
    aes(x = logTitre, ymin = `2.5%`, ymax = `97.5%`, fill = model),
    inherit.aes = FALSE, alpha = 0.5
  )

# save_plot(preds_plot, "predsplot", FALSE)
