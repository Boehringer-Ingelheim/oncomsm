.tbl_to_stan_data <- function(model, data, n_additional) {
  UseMethod(".tbl_to_stan_data")
}

.prior_to_stan_data <- function(model) {
  UseMethod(".prior_to_stan_data")
}
