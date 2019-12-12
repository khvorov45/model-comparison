# Notes for model comparison

Arseniy Khvorov.
Created 2019/09/12
Last edit 2019/12/12

## Dependencies

All `.R` files (except those in `tests`) are (should be) independent but expect the working directory to be that of the `.Rproj` file.

The only local external dependency is `hanam.csv` in the `data` folder.


## Comparing to log(5)

- Cox would be the same as shifting X by log(5) since e^(bX) would evaluate
to 1 at log(5), so the hazard would correspond to the baseline hazard

- Logistic and scaled logit are not like Cox. 
Shifting X by log(5) would change b0 which
would now correspond to protection at log(5) and not log(1). Have to
actually compute protection at log(5) and divide any protection level by it.

## Effect of censoring

- Introduces bias (betas move toward zero) but not a lot of it 
(<10% difference as compared to fitting to true values). They also inflate 
errors (by ~50% for betas and ~5-fold for lambda).

- Does multiple imputation/bayesian fitting solve these?

  - Not using multiple imputation for now

## Plausible parameter ranges

- lambda: can be anything

- Titre intercept and covariate coefficients. For protective titres, in order
for the probability of infection to drop within the expected range, the 
intercept term has to be negative. The more negative, the more positive the
titre coefficient should be. From previous research we could probably expect the
coefficient to be in the single digits so the intercept is likely somewhere
between -30 and 0.

- HI titres. We are assuming a normal distribution for the underlying 
unobservable log titres. Mean around 0 - almost everyone will have 
undetectable titres. Mean around 4 - majority will have fairly high titres.
SD around 5 - more people will have titres over 1280 than in any other category
(except undetectable). True mean and sd are likely to lie somwhere in early
single digits.

## Stan

Coupling the scaled logit model with censored data makes stan not produce usable results
