# get unique identifiers that are reproducible
get_identifier <- function(n = 1, exclude = NULL) {
  res <- character(n)
  duplicate <- rep(TRUE, n)
  while (any(duplicate)) {
    res[duplicate] <- sprintf(
      "ID%08i",
      sample.int(1e7, sum(duplicate), replace = FALSE)
    )
    duplicate <- res %in% exclude
  }
  res
}



# get plotting grid
get_dt_grid <- function(model,
                        parameter_sample, dt_interval, dt_n_grid, dt_expand,
                        seed) {
  if (is.null(dt_interval)) {
    dt_max <- sample_predictive(
        model,
        n_per_group = rep(1e3L, length(model$group_id)),
        sample = parameter_sample, nsim = 1, seed = seed, as_mstate = TRUE
      ) %>%
      arrange(.data$t_sot, .data$subject_id, .data$t_min) %>%
      mutate(
        dt = pmax(.data$t_max - .data$t_sot)
      ) %>%
      pull("dt") %>%
      {stats::quantile(.[is.finite(.)], probs = 0.9)} %>% # nolint cut outliers
      as.numeric()
    dt_interval <- c(0, dt_max * dt_expand)
  }
  dt_grid <- seq(dt_interval[1] + 0.1, dt_interval[2], length.out = dt_n_grid)
  return(dt_grid)
}



get_mu_sigma <- function(q05, q95) {
  f <- function(x) {
    mu <- x[1]
    sigma <- x[2]
    stats::qlnorm(c(0.05, 0.95), mu, sigma)
  }
  res <- stats::optim(
    c(1, 0.5),
    function(x) sum((f(x) - c(q05, q95))^2),
    lower = c(-Inf, 0.001),
    method = "L-BFGS-B"
  )
  return(tibble(
    mu = res$par[1],
    sigma = res$par[2]
  ))
}
