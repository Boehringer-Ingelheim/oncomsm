<!-- badges: start -->

[![R-CMD-check](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/check.yml/badge.svg?branch=main)](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/check.yml)
[![R-CMD-check-windows](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/rcmd-check-windows.yml/badge.svg?branch=main)](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/rcmd-check-windows.yml)
[![pages-build-deployment](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/pages/pages-build-deployment)
[![metacran version](https://www.r-pkg.org/badges/version/oncomsm)](https://cran.r-project.org/package=oncomsm)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](https://github.com/Boehringer-Ingelheim/oncomsm)

<!-- badges: end -->



# Bayesian multi-state models for oncology

**Tl;dr:** This package implements methods to dynamically predict response 
and progression of individuals in early oncology trials using parametric
[multi-state models](https://boehringer-ingelheim.github.io/oncomsm/articles/multi-state-model-for-early-oncology.html) and Bayesian inference. 
This allows the dynamic computation of "probability of success" for a wide 
range of success criteria. 


## Installation

```{r}
# install.packages("remotes")
remotes::install_github("https://github.com/Boehringer-Ingelheim/oncomsm")
```


## Documentation

The package documentation is hosted [here](https://boehringer-ingelheim.github.io/oncomsm/).


## Contributing

### Giving feedback

Your feedback and error reports are highly valued!


### Updating stan models

The stan models contained in `inst/stan` are not automatically updated to avoid
taking a dependency on the `rstantools` package.  
After modifying or adding new models, run 
```{r}
rstantools::rstan_config()
```

and silence `R/stanmodels.R` via `capture.output()`.
