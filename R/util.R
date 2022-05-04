.tbl_to_stan_data <- function(model, data) {
  UseMethod(".tbl_to_stan_data")
}

.prior_to_stan_data <- function(model) {
  UseMethod(".prior_to_stan_data")
}

extract_col_from_stan <- function(res, name) {
  rstan::extract(res, name)[[1]] %>%
    t() %>%
    {colnames(.) <- 1:ncol(.); .} %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(everything(), names_to = "iter", values_to = name) %>%
    dplyr::mutate(iter = as.integer(.data$iter))
}
