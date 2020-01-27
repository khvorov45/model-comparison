# model-comparison

Comparison of models for antibody data

## Directories

* `cox-tarprop` - simulations to see if the Cox model produces reliable results when time at risk is proportional to time of follow-up.

* `cox-tarprop-summary` - summary of `cox-tarprop` simulations.

* `cox-tarprop-plot` - plot of `cox-tarprop-summary`.

* `curve-cox` - plot showing cox model assumptions.

* `curve-models` - plot showing model curves (cox, logit and scaled logit).

* `data` - contains all real data used.

* `data-raw` - contains data that needed modifying before moving to `data`. Conventions:
  * File names have the format `[study]-[dataset]` where `study` is the name of the study (e.g. `kiddyvaxmain`) and `dataset` is the name of the dataset that comes from the study (e.g. `serology`). If `dataset` part is not present, all study data is in `[study].csv` file.

* `data-plot` - plot of hanam data.

* `edu-article` - example of an educational article.

* `fit-bayesian` - bayseian fit to the Hanam data. Used to compare to maximum
likelihood fit to the midpoints.

* `fit-bayesian-summary` - summary of fit (calculated infection/protection distributions)

* `fit-bayesian-diag` - trace plots and density plots.

* `fit-bayesian-plot` - infection and protection curves.

* `fit-sclr` - scaled logit model fit to hanam data.

* `fit-sclr-boot` - bootstrap samples of the scaled logit model fit to hanam
data.

* `logistic` - logistic and scaled logit simulations.

* `logistic-summary` - summary of simulations in `logistic`.

* `logistic-plot` - plots of summaries in `logistic-summary`.

* `plausible-titres` - Plausible values for the titre distribution.

* `plausible-sclr` - Plausible values for the scaled logit model parameters.

* `presentation` - presentation I made on this.

* `tests` - function tests. Expect the functions to be in the global envirnment. Rstudio's "Run tests" will not work because it will not look there.

* `writeup` - all the results and discussion in writing.
