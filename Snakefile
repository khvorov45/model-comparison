subworkflow cox_tarprop:
    workdir:
        "cox-tarprop"

subworkflow sclr_lowbase:
    workdir:
        "sclr-lowbase"

subworkflow logit_sim:
    workdir:
        "logit-sim"

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
        "data-plot/hanam-hi-gen-scatter.pdf",
        "data-plot/hanam-hi-scatter.pdf",
        "data-plot/hanam-hi-summ-gen.pdf",
        "data-plot/hanam-hi-summ-h3n2.pdf",
        "data-plot/hanam-hi-summ.pdf",
        "data/sophia.csv",
        "data/kiddyvaxmain.csv",
        "data/kiddyvaxmain-summ.csv",
        "data/kiddyvaxmain-swab.csv"

rule install_deps:
    input:
        "renv.lock"
    output:
        ".deps-installed"
    shell:
        """Rscript -e 'renv::restore();file.create(".deps-installed")'"""

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

rule data_plot:
    input:
        ".deps-installed",
        "data-plot/hanam-plot.R",
        "data/read_data.R",
        "data/hanam-hi-exp.csv",
        "data/hanam-hi-gen.csv",
        "data/hanam-hi-summ.csv"
    output:
        "data-plot/hanam-hi-gen-scatter.pdf",
        "data-plot/hanam-hi-scatter.pdf",
        "data-plot/hanam-hi-summ-gen.pdf",
        "data-plot/hanam-hi-summ-h3n2.pdf",
        "data-plot/hanam-hi-summ.pdf"
    shell:
        "Rscript data-plot/hanam-plot.R"
