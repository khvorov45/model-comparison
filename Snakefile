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
        "plausible/plausible-titres.pdf"

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
