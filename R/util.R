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

.nodata <- function(model) { # nolint
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  tibble(
    subject_id = integer(),
    group_id = integer(),
    from = integer(),
    to = integer(),
    t_min = numeric(),
    t_max = numeric(),
    t_sot = numeric()
  )
}



# create a data set with no observed data
.emptydata <- function(model, n_per_group, seed = NULL) { # nolint
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (!is.null(seed)) {
    set.seed(seed) # nocov
  }
  n <- sum(n_per_group)
  group_ids <- model$group_id
  rr <- model$recruitment_rate
  res <- tibble()
  for (i in seq_along(group_ids)) {
    if (n_per_group[i] < 1) {
      next
    }
    subject_ids <- get_identifier(n = n_per_group[i])
    recruitment_times <- cumsum(stats::rexp(n_per_group[i], rate = rr[i]))
    res <- bind_rows(res, tibble(
      subject_id = subject_ids,
      group_id = group_ids[i],
      t = recruitment_times,
      state = "stable" # first visits are always stable
    ))
  }
  res <- arrange(res, t)
  attr(res, "isemptydata") <- TRUE
  return(res)
}



# generic for digesting data into rstan-ready list
data2standata <- function(data, model) { # nolint
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  # first prepare any data (if available)
  lst_stan_data <- data %>%
    transmute(
      group_id = as.integer(factor(.data$group_id,
                                   levels = model$group_id
      )),
      subject_id = as.integer(factor(.data$subject_id)),
      from = as.integer(factor(.data$from, levels = model$states)),
      to = if_else(
          is.na(.data$to),
          "unknown",
          as.character(.data$to)
        ) %>%
        factor(levels = c(model$states, "unknown")) %>%
        as.integer(),
      t_min = .data$t_min - .data$t_sot,
      t_max = .data$t_max - .data$t_sot,
    )
  # adjust data for numerical stability
  lst_stan_data <- lst_stan_data %>%
    group_by(.data$subject_id) %>%
    mutate(
      # within transition spacing
      t_min = pmax(.data$t_min, model$visit_spacing[.data$group_id] / 2),
      t_max = pmax(.data$t_max, .data$t_min + 1 / 30), # 1 day
      # align within subject spacing
      t_min = pmax(.data$t_min, lag(.data$t_max, default = 0)),
      t_max = pmax(.data$t_max, .data$t_min + 1 / 30)
    ) %>%
    arrange(.data$subject_id, .data$from) %>%
    as.list()
  # make sure everything is an array
  for (i in seq_along(lst_stan_data)) {
    lst_stan_data[[i]] <- as.array(lst_stan_data[[i]])
  }
  # data meta information
  lst_stan_data$N <- nrow(data)
  lst_stan_data$N_subjects <- length(unique(data$subject_id))
  # model meta information and priors
  lst_stan_data$M_groups <- length(model$group_id)
  lst_stan_data$maximal_time <- model$maximal_time
  lst_prior_parameters <- list(
    p_mean = model$prior$p[, "mean"] %>% as.array(),
    p_n = model$prior$p[, "n"] %>% as.array(),
    p_eta = model$prior$p[, "eta"] %>% as.array(),
    p_min = model$prior$p[, "min"] %>% as.array(),
    p_max = model$prior$p[, "max"] %>% as.array(),
    median_t_mu = matrix(NA_real_, nrow = lst_stan_data$M_groups, ncol = 3) %>%
      matrix(ncol = 3),
    median_t_sigma = matrix(NA_real_, nrow = lst_stan_data$M_groups, ncol = 3) %>%  # nolint
      matrix(ncol = 3),
    shape_mu = matrix(NA_real_, nrow = lst_stan_data$M_groups, ncol = 3),
    shape_sigma = matrix(NA_real_, nrow = lst_stan_data$M_groups, ncol = 3)
  )
  # calculate mu/sigma for log-normal priors based on modes/90% quantiles
  for (g in 1:lst_stan_data$M_groups) {
    for (trans in 1:3) {
      res <- get_mu_sigma(
        model$prior$median_t[g, trans, "q05"],
        model$prior$median_t[g, trans, "q95"]
      )
      lst_prior_parameters$median_t_mu[g, trans] <- res$mu
      lst_prior_parameters$median_t_sigma[g, trans] <- res$sigma
      res <- get_mu_sigma(
        model$prior$shape[g, trans, "q05"],
        model$prior$shape[g, trans, "q95"]
      )
      lst_prior_parameters$shape_mu[g, trans] <- res$mu
      lst_prior_parameters$shape_sigma[g, trans] <- res$sigma
    }
  }
  lst_stan_data <- c(lst_stan_data, lst_prior_parameters)
  return(lst_stan_data)
}
