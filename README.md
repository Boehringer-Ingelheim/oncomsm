<!-- badges: start -->

[![R-CMD-check](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/check.yml/badge.svg?branch=main)](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/check.yml)
[![metacran version](https://www.r-pkg.org/badges/version/oncomsm)](https://cran.r-project.org/package=oncomsm)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](https://github.com/Boehringer-Ingelheim/oncomsm)

<!-- badges: end -->



# Bayesian multi-state models for oncology

Tl;dr: This package implements methods to dynamically predict response 
and progression of individuals in early oncology trials using parametric
multi-state models and Bayesian inference. 
This allows the dynamic computation of Probability of Success for a wide 
range of success criteria.


## Installation

```{r}
# install.packages("remotes")
remotes::install_github("https://github.com/Boehringer-Ingelheim/oncomsm")
```


## Background

...

## Compiling from source

The stan models contained in `inst/stan` are not automatically updated to avoid
taking a dependency on the `rstantools` package. 
After modifying or adding new models, run 
```{r}
rstantools::rstan_config()
```

and silence `R/stanmodels.R` via `capture.output()`.
