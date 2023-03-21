<!-- badges: start -->
[![Linux](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/linux.yml/badge.svg?branch=main)](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/linux.yml)
[![Windows](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/windows.yml/badge.svg?branch=main)](https://github.com/Boehringer-Ingelheim/oncomsm/actions/workflows/windows.yml)
[![metacran version](https://www.r-pkg.org/badges/version-last-release/oncomsm)](https://cran.r-project.org/package=oncomsm)
[![metacran version](https://cranlogs.r-pkg.org/badges/grand-total/oncomsm)](https://cran.r-project.org/package=oncomsm)
[![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)](https://github.com/Boehringer-Ingelheim/oncomsm)
<!-- badges: end -->

test-cz

# Bayesian multi-state models for early oncology

The R package **oncomsm** implements methods to dynamically predict response 
and progression of individuals in early oncology trials using parametric
[multi-state models](https://boehringer-ingelheim.github.io/oncomsm/articles/oncomsm.html) and Bayesian inference. 
This allows the [dynamic computation of "probability of success"](https://boehringer-ingelheim.github.io/oncomsm/articles/web_only/application-to-probability-of-success.html) for a wide 
range of success (or "go") criteria. 
For instance, the [bhmbasket](https://cran.r-project.org/package=bhmbasket) R package [can be used to define study success based on Bayesian hierarchical models](https://boehringer-ingelheim.github.io/oncomsm/articles/bhmbasket-integration.html).


## Installation

The development version can be installed from this repository.

```{r}
install.packages("oncomsm")
```


The development version can be installed from this repository.

```{r}
# install.packages("remotes")
remotes::install_github("https://github.com/Boehringer-Ingelheim/oncomsm")
```


## Documentation

The package documentation is hosted [here](https://boehringer-ingelheim.github.io/oncomsm/).


## Contributing

See the [contributing guidelines](https://boehringer-ingelheim.github.io/oncomsm/CONTRIBUTING.html).
