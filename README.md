# Categorical Regression Misc Functions

**Miscellaneous functions for interpreting categorical regression models in R.**

This script contains miscellaneous functions for interpreting the output of categorical regression models in R. The functions implemented so far work for ordinal and multinomial logit models from the [MASS](https://cran.r-project.org/web/packages/MASS/MASS.pdf) and [mlogit](https://cran.r-project.org/web/packages/mlogit/mlogit.pdf) R packages respectively. When it comes to implementing and interpreting categorical regression models Stata has a much richer and more easily accessible set of tools for the task than R. The functions in this project are written to help a class of relatively new R users more easily make sense of the output of these models.

### Summary of functions

So far these only work for the logit models (not their probit counterparts). The functions rely on the **MASS** & **mlogit** packages and some of their dependencies.

1. `polrLogitStdCoef()` _# X standardized coefficients for ordinal logit models_
    * `polrLogitStdCoef(polrObject, help = FALSE, digits = 4)`
2. `mlogTestR()` _# does Wald and Likelihood ratio tests of joint significance for multinomial logit models from the mlogit package_
    * `mlogTestR(mlogitObject, reflevel, type = "wald", digits = 4)`
