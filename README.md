# Fully Bayesian generative models for delayed event rates with time-to-first-event endpoints

...

## Compiling from source

The stan models contained in `inst/stan` are not automatically updated to avoid a
dependency on the `rstantools` package. 
After modifying or adding new models, run 
```{r}
rstantools::rstan_config()
```
add
```
CXX14FLAGS += -mtune=native -O3 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2
```
to `src/Makevars.win` and silence `R/stanmodels.R` via `capture.output()`.
