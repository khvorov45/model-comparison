# Plotting protections curves

library(tidyverse)

# Directories used
fit_dir <- here::here("fit")
preds_plot_dir <- here::here("preds-plot")
data_dir <- here::here("data")

# Functions ===================================================================

source(file.path(data_dir, "read_data.R"))

read_pred <- function(name) {
  read_csv(
    file.path(fit_dir, glue::glue("{name}.csv")),
    col_types = cols()
  )
}

recode_viruses <- function(dat) {
  dat %>%
    mutate(
      virus_lbl = factor(
        virus,
        levels = c(
          "H1N1pdm", "A H1pdm", "h1pdm", "H3N2", "h3", "bvic", "B Vic"
        ),
        labels = c(
          "H1N1pdm", "H1N1pdm", "H1N1pdm", "H3N2", "H3N2", "B Vic", "B Vic"
        )
      )
    )
}

common_plot_els <- function() {
  xbreaks <- c(5 * 2^(0:10))
  ybreaks <- seq(0, 1, 0.1)
  list(
    ggdark::dark_theme_bw(verbose = FALSE),
    theme(
      strip.background = element_blank(),
      panel.spacing = unit(0, "null"),
      axis.text.x = element_text(angle = 45, hjust = 1),
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

plot_pred <- function(dat, facets = "vir",
                      xmin = 5, xmax = 5120, ylab = "Relative protection") {
  facets <- rlang::arg_match(facets, c("vir", "virpop", "pop", "none"))
  if (facets == "virpop") {
    faceting <- facet_grid(population ~ virus_lbl)
  } else if (facets == "vir") {
    faceting <- facet_wrap(~virus_lbl, nrow = 1)
  } else if (facets == "pop") {
    faceting <- facet_wrap(~population, nrow = 1)
  } else {
    faceting <- NULL
  }
  dat %>%
    ggplot(aes(loghi, prot, ymin = prot_low, ymax = prot_high)) +
    common_plot_els() +
    ylab(ylab) +
    faceting +
    coord_cartesian(xlim = c(log(xmin), log(xmax)), ylim = c(0, 1))
}

plot_pred_inf <- function(outsum, data, facets = "virpop",
                          xmin = 5, xmax = 5120, ymin = 0, ymax = 0.5) {
  facets <- rlang::arg_match(facets, c("vir", "virpop", "pop", "none"))
  if (facets == "virpop") {
    facets_spec <- list(
      facetscales::facet_grid_sc(
        population ~ virus_lbl,
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
      ),
      coord_cartesian(xlim = c(log(xmin), log(xmax))),
      theme(panel.spacing.y = unit(5, "points"))
    )
  } else if (facets == "vir") {
    facets_spec <- list(
      facet_wrap(~virus_lbl, nrow = 1),
      coord_cartesian(xlim = c(log(xmin), log(xmax)), ylim = c(ymin, ymax))
    )
  } else if (facets == "pop") {
    facets_spec <- list(
      facet_wrap(~population, nrow = 1),
      coord_cartesian(xlim = c(log(xmin), log(xmax)), ylim = c(ymin, ymax))
    )
  } else {
    facets_spec <- list()
  }
  data <- filter(data, virus_lbl %in% unique(outsum$virus_lbl))
  outsum %>%
    ggplot(aes(loghi, fit, ymin = fit_low, ymax = fit_high)) +
    ylab("Infection probability") +
    common_plot_els() +
    facets_spec +
    geom_point(
      data = data, mapping = aes(loghimid, inf_prop), inherit.aes = FALSE,
      shape = 18
    ) +
    ggrepel::geom_text_repel(
      data = data, mapping = aes(loghimid, inf_prop, label = ntot),
      inherit.aes = FALSE
    )
}

widen_preds <- function(preds, ...) {
  preds %>%
    pivot_wider(
      names_from = "prob_type",
      values_from = c("prob_lb", "prob_med", "prob_ub", ...)
    ) %>%
    rename(
      fit_low = prob_lb_inf, fit = prob_med_inf, fit_high = prob_ub_inf,
      prot_low = prob_lb_prot, prot = prob_med_prot, prot_high = prob_ub_prot
    ) %>%
    recode_viruses()
}

add_priors <- function(plot) {
  plot +
    geom_line(aes(y = prob_lb_prior), lty = "3333", col = "gray50") +
    geom_line(aes(y = prob_ub_prior), lty = "3333", col = "gray50")
}

save_plot <- function(pl, name, width = 12, height = 7.5) {
  ggdark::ggsave_dark(
    file.path(preds_plot_dir, glue::glue("{name}.pdf")), pl,
    device = "pdf", width = width, height = height, units = "cm"
  )
}

# Script ======================================================================

kv_main_summ <- read_data("kiddyvaxmain-summ") %>% recode_viruses()
han_hi_summ <- read_data("hanam-hi-summ") %>%
  filter(virus != "H1N1seas") %>%
  recode_viruses()

kv_cox_preds <- read_pred("kiddyvaxmain-preds-cox")
sophia_preds <- read_pred("sophia-preds-cox")
han_hi_lr <- read_pred("hanam-hi-preds-lr") %>% recode_viruses()
han_hi_sclr <- read_pred("hanam-hi-preds-sclr") %>%
  recode_viruses() %>%
  rename(
    prot = prot_point, prot_low = prot_l, prot_high = prot_u, loghi = loghimid
  )
han_hi_lr_boot <- read_pred("hanam-hi-preds-lr-boot") %>% widen_preds()
han_hi_sclr_boot <- read_pred("hanam-hi-preds-sclr-boot") %>%
  widen_preds()
han_hi_sclr_bayes <- read_pred("hanam-hi-preds-sclr-bayesian") %>%
  widen_preds("prob_lb_prior", "prob_ub_prior", "prob_med_prior") %>%
  mutate(
    prob_lb_prior = prob_lb_prior_prot, prob_ub_prior = prob_ub_prior_prot,
    prob_med_prior = prob_med_prior_prot
  )

kvm_lr <- read_pred("kiddyvaxmain-preds-lr") %>% recode_viruses()
kvm_sclr <- read_pred("kiddyvaxmain-preds-sclr") %>%
  recode_viruses() %>%
  rename(
    prot = prot_point, prot_low = prot_l, prot_high = prot_u, loghi = loghimid
  )
kvm_lr_boot <- read_pred("kiddyvaxmain-preds-lr-boot") %>% widen_preds()
kvm_sclr_boot <- read_pred("kiddyvaxmain-preds-sclr-boot") %>%
  widen_preds()
kvm_sclr_bayes <- read_pred("kiddyvaxmain-preds-sclr-bayesian") %>%
  widen_preds("prob_lb_prior", "prob_ub_prior", "prob_med_prior") %>%
  mutate(
    prob_lb_prior = prob_lb_prior_prot, prob_ub_prior = prob_ub_prior_prot,
    prob_med_prior = prob_med_prior_prot
  )

all_plots <- list(
  "kiddyvaxmain-cox" = plot_pred(kv_cox_preds),
  "kiddyvaxmain-lr" = plot_pred(kvm_lr, ylab = "Protection"),
  "kiddyvaxmain-lr-boot" = plot_pred(kvm_lr_boot, ylab = "Protection"),
  "kiddyvaxmain-lr-boot-inf" = plot_pred_inf(
    kvm_lr_boot, kv_main_summ, "vir",
    ymax = 0.2
  ),
  "kiddyvaxmain-lr-boot-rel" = plot_pred(
    mutate(
      kvm_lr_boot,
      prot = prob_med_prot_rel, prot_low = prob_lb_prot_rel,
      prot_high = prob_ub_prot_rel
    )
  ),
  "kiddyvaxmain-lr-boot-rel" = plot_pred(
    mutate(
      filter(kvm_lr_boot, virus_lbl == "B Vic"),
      prot = prob_med_prot_rel, prot_low = prob_lb_prot_rel,
      prot_high = prob_ub_prot_rel
    ),
    "none"
  ),
  "kiddyvaxmain-sclr" = plot_pred(kvm_sclr, ylab = "Protection"),
  "kiddyvaxmain-sclr-boot" = plot_pred(
    kvm_sclr_boot,
    ylab = "Protection"
  ),
  "kiddyvaxmain-sclr-boot-inf" = plot_pred_inf(
    kvm_sclr_boot, kv_main_summ, "vir",
    ymax = 0.15
  ),
  "kiddyvaxmain-sclr-bayesian" =
    plot_pred(kvm_sclr_bayes, ylab = "Protection") %>%
      add_priors(),
  "kiddyvaxmain-sclr-bayesian-inf" = plot_pred_inf(
    mutate(
      kvm_sclr_bayes,
      prob_lb_prior = prob_lb_prior_inf, prob_ub_prior = prob_ub_prior_inf,
      prob_med_prior = prob_med_prior_inf
    ),
    kv_main_summ, "vir",
    ymax = 0.2
  ) %>% add_priors(),
  "kiddyvaxmain-lr-bvic" = plot_pred(
    filter(kvm_lr, virus_lbl == "B Vic"), "none",
    ylab = "Protection"
  ),
  "kiddyvaxmain-lr-inf" = plot_pred_inf(
    kvm_lr, kv_main_summ, "vir",
    ymax = 0.2
  ),
  "kiddyvaxmain-cox-bvic" = plot_pred(
    filter(kv_cox_preds, virus_lbl == "B Vic"),
    facets = "none",
  ),
  "sophia-cox-og" = plot_pred(
    filter(sophia_preds, model == "sophia") %>%
      mutate(prot_low = prot_low_wrong, prot_high = prot_high_wrong)
  ),
  "sophia-cox-fixci" = plot_pred(
    filter(sophia_preds, model == "sophia")
  ),
  "sophia-cox-fixci-fixmod" = plot_pred(
    filter(sophia_preds, model == "me")
  ),
  "hanam-hi-lr" = plot_pred(han_hi_lr, "virpop", ylab = "Protection"),
  "hanam-hi-lr-boot" = plot_pred(han_hi_lr_boot, "virpop", ylab = "Protection"),
  "hanam-hi-lr-boot-inf" = plot_pred_inf(han_hi_lr_boot, han_hi_summ),
  "hanam-hi-lr-boot-rel" = plot_pred(
    mutate(
      han_hi_lr_boot,
      prot = prob_med_prot_rel, prot_low = prob_lb_prot_rel,
      prot_high = prob_ub_prot_rel
    ),
    "virpop"
  ),
  "hanam-hi-lr-boot-rel-h3" = plot_pred(
    mutate(
      filter(han_hi_lr_boot, virus_lbl == "H3N2"),
      prot = prob_med_prot_rel, prot_low = prob_lb_prot_rel,
      prot_high = prob_ub_prot_rel
    ),
    "pop"
  ),
  "hanam-hi-sclr" = plot_pred(han_hi_sclr, "virpop", ylab = "Protection"),
  "hanam-hi-sclr-boot" = plot_pred(
    han_hi_sclr_boot, "virpop",
    ylab = "Protection"
  ),
  "hanam-hi-sclr-boot-inf" = plot_pred_inf(
    han_hi_sclr_boot, han_hi_summ, "virpop",
    ymax = 0.3
  ),
  "hanam-hi-sclr-bayesian" =
    plot_pred(han_hi_sclr_bayes, "virpop", ylab = "Protection") %>%
      add_priors(),
  "hanam-hi-sclr-bayesian-inf" = plot_pred_inf(
    mutate(
      han_hi_sclr_bayes,
      prob_lb_prior = prob_lb_prior_inf, prob_ub_prior = prob_ub_prior_inf,
      prob_med_prior = prob_med_prior_inf
    ),
    han_hi_summ
  ) %>% add_priors(),
  "hanam-hi-lr-h3" = plot_pred(
    filter(han_hi_lr, virus == "H3N2"), "pop",
    ylab = "Protection"
  ),
  "hanam-hi-lr-inf" = plot_pred_inf(han_hi_lr, han_hi_summ)
)

iwalk(all_plots, save_plot)
