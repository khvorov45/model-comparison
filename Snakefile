# Subworkflows ================================================================

subworkflow cox_tarprop:
    workdir:
        "cox-tarprop"

subworkflow sclr_lowbase:
    workdir:
        "sclr-lowbase"

subworkflow logit_sim:
    workdir:
        "logit-sim"

# Rule to generate everything =================================================

rule all:
    input:
        ".deps-installed",
        cox_tarprop("plot/long.pdf"),
        cox_tarprop("plot/risk.pdf"),
        sclr_lowbase("sim-plot/plot2.pdf"),
        sclr_lowbase("sim-plot/plot3.pdf"),
        logit_sim("plot/vary_nsam.pdf"),
        logit_sim("plot/vary_lambda.pdf"),
        logit_sim("plot/vary_nsam_mean.pdf"),
        logit_sim("plot/vary_nsam_se.pdf"),
        logit_sim("plot/predsplot.pdf"),
        expand("curve-cox/timeplot_{i}.pdf", i = range(0, 5)),
        expand("plausible/sclr/pl_{i}.png", i = range(1, 21)),
        "plausible/plausible-titres.pdf",
        "fig-studies/fig-studies.pdf",
        "data-plot/hanam-hi-gen-scatter.pdf",
        "data-plot/hanam-hi-scatter.pdf",
        "data-plot/hanam-hi-summ-gen.pdf",
        "data-plot/hanam-hi-summ-h3n2.pdf",
        "data-plot/hanam-hi-summ.pdf",
        "data-plot/kiddyvax-main-summ.pdf",
        "data-plot/kiddyvax-main-swab.pdf",
        "data-table/hanam-hi-tbl1.csv",
        "fit/sophia-preds-cox.csv",
        "fit/kiddyvaxmain-preds-cox.csv",
        "fit/hanam-hi-preds-lr.csv",
        "fit/kiddyvaxmain-preds-lr.csv",
        "fit/hanam-hi-preds-sclr.csv",
        "fit/kiddyvaxmain-preds-sclr.csv",
        "fit/out-sclr-bayesian/hanam-hi-exp-h1pdm.csv",
        "fit/out-sclr-bayesian/hanam-hi-exp-h3.csv",
        "fit/out-sclr-bayesian/hanam-hi-gen-h1pdm.csv",
        "fit/out-sclr-bayesian/hanam-hi-gen-h3.csv",
        "fit/out-sclr-bayesian/kiddyvaxmain-bvic.csv",
        "fit/out-sclr-bayesian/kiddyvaxmain-h1pdm.csv",
        "fit/out-logistic-boot/hanam-hi.csv",
        "fit/out-logistic-boot/kiddyvaxmain.csv",
        "fit/out-sclr-boot/hanam-hi.csv",
        "fit/out-sclr-boot/kiddyvaxmain.csv"

# Dependencies ================================================================

rule install_deps:
    input:
        "renv.lock"
    output:
        ".deps-installed"
    shell:
        """Rscript -e 'renv::restore();file.create(".deps-installed")'"""

# Misc plots ==================================================================

rule curve_cox:
    input:
        ".deps-installed",
        "curve-cox/curve-cox.R"
    output:
        expand("curve-cox/timeplot_{i}.pdf", i = range(0, 5))
    shell:
        "Rscript curve-cox/curve-cox.R"

rule plausible:
    input:
        ".deps-installed",
        "plausible/plausible.R"
    output:
        expand("plausible/sclr/pl_{i}.png", i = range(1, 21)),
        "plausible/plausible-titres.pdf"
    shell:
        "Rscript plausible/plausible.R"

rule fig_studies:
    input:
        ".deps-installed",
        "fig-studies/fig-studies.R"
    output:
        "fig-studies/fig-studies.pdf"
    shell:
        "Rscript fig-studies/fig-studies.R"

# Raw data process ============================================================

rule hanam:
    input:
        ".deps-installed",
        "data/hanam.R",
        "data-raw/hanam.csv"
    output:
        "data/hanam-hi-exp.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-summ.csv"
    shell:
        "Rscript data/hanam.R"

rule sophia:
    input:
        ".deps-installed",
        "data/sophia.R"
    output:
        "data/sophia.csv"
    shell:
        "Rscript data/sophia.R"

rule kiddyvax:
    input:
        ".deps-installed",
        "data/kiddyvaxmain.R",
        "data-raw/kiddyvaxmain-serology.csv",
        "data-raw/kiddyvaxmain-swab.csv"
    output:
        "data/kiddyvaxmain.csv",
        "data/kiddyvaxmain-summ.csv",
        "data/kiddyvaxmain-swab.csv"
    shell:
        "Rscript data/kiddyvaxmain.R"

# Data plot ===================================================================

rule data_plot:
    input:
        ".deps-installed",
        "data-plot/data-plot.R",
        "data/read_data.R",
        "data/hanam-hi-exp.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-summ.csv"
    output:
        "data-plot/hanam-hi-gen-scatter.pdf",
        "data-plot/hanam-hi-scatter.pdf",
        "data-plot/hanam-hi-summ-gen.pdf",
        "data-plot/hanam-hi-summ-h3n2.pdf",
        "data-plot/hanam-hi-summ.pdf",
        "data-plot/kiddyvax-main-summ.pdf",
        "data-plot/kiddyvax-main-swab.pdf"
    shell:
        "Rscript data-plot/data-plot.R"

# Data table ==================================================================

rule data_table:
    input:
        ".deps-installed",
        "data-table/data-table.R",
        "data/read_data.R",
        "data/hanam-hi-exp.csv",
        "data/hanam-hi-gen.csv"
    output:
        "data-table/hanam-hi-tbl1.csv"
    shell:
        "Rscript data-table/data-table.R"

# Model fitting ===============================================================

rule fit_cox:
    input:
        ".deps-installed",
        "fit/fit-cox.R",
        "data/read_data.R",
        "data/kiddyvaxmain.csv",
        "data/sophia.csv"
    output:
        "fit/kiddyvaxmain-preds-cox.csv",
        "fit/sophia-preds-cox.csv"
    shell:
        "Rscript fit/fit-cox.R"

rule fit_lr:
    input:
        ".deps-installed",
        "fit/fit-lr.R",
        "data/read_data.R",
        "data/kiddyvaxmain.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-exp.csv"
    output:
        "fit/kiddyvaxmain-preds-lr.csv",
        "fit/hanam-hi-preds-lr.csv"
    shell:
        "Rscript fit/fit-lr.R"

rule fit_sclr:
    input:
        ".deps-installed",
        "fit/fit-sclr.R",
        "data/read_data.R",
        "data/kiddyvaxmain.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-exp.csv"
    output:
        "fit/kiddyvaxmain-preds-sclr.csv",
        "fit/hanam-hi-preds-sclr.csv"
    shell:
        "Rscript fit/fit-sclr.R"

rule fit_sclr_bayesian:
    input:
        ".deps-installed",
        "fit/fit-sclr-bayesian.R",
        "data/read_data.R",
        "data/kiddyvaxmain.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-exp.csv"
    output:
        protected("fit/out-sclr-bayesian/hanam-hi-exp-h1pdm.csv"),
        protected("fit/out-sclr-bayesian/hanam-hi-exp-h3.csv"),
        protected("fit/out-sclr-bayesian/hanam-hi-gen-h1pdm.csv"),
        protected("fit/out-sclr-bayesian/hanam-hi-gen-h3.csv"),
        protected("fit/out-sclr-bayesian/kiddyvaxmain-bvic.csv"),
        protected("fit/out-sclr-bayesian/kiddyvaxmain-h1pdm.csv")
    shell:
        "Rscript fit/fit-sclr-bayesian.R"

rule fit_logistic_boot:
    input:
        ".deps-installed",
        "fit/fit-logistic-boot.R",
        "data/read_data.R",
        "data/kiddyvaxmain.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-exp.csv"
    output:
        protected("fit/out-logistic-boot/hanam-hi.csv"),
        protected("fit/out-logistic-boot/kiddyvaxmain.csv")
    shell:
        "Rscript fit/fit-logistic-boot.R"

rule fit_sclr_boot:
    input:
        ".deps-installed",
        "fit/fit-sclr-boot.R",
        "data/read_data.R",
        "data/kiddyvaxmain.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-exp.csv"
    output:
        protected("fit/out-sclr-boot/hanam-hi.csv"),
        protected("fit/out-sclr-boot/kiddyvaxmain.csv")
    shell:
        "Rscript fit/fit-sclr-boot.R"

# Data plots ==================================================================

rule preds_plot:
    input:
        ".deps-installed",
        "preds-plot/preds-plot.R",
        "fit/kiddyvaxmain-preds-cox.csv",
        "fit/sophia-preds-cox.csv"
    output:
        "preds-plot/kiddyvaxmain-cox-bvic.pdf",
        "preds-plot/kiddyvaxmain-cox.pdf",
        "preds-plot/sophia-cox-og.pdf",
        "preds-plot/sophia-cox-fixci.pdf",
        "preds-plot/sophia-cox-fixci-fixmod.pdf"
    shell:
        "Rscript preds-plot/preds-plot.R"
