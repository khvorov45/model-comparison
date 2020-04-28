subworkflow cox_tarprop:
    workdir:
        "cox-tarprop"

rule all:
    input:
        cox_tarprop("plot/long.pdf"),
        cox_tarprop("plot/risk.pdf")
