# Fully Bayesian generative models for delayed event rates with time-to-first-event endpoints

...

## Compiling from source

The stan models contained in `inst/stan` are not automatically updated to avoid a
dependency on the `rstantools` package. 
After modifying or adding new models, run 
```{r}
rstantools::rstan_config()
```
and silence `R/stanmodels.R` via `capture.output()`.
