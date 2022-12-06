#' A stable-response-progression model
#'
#' `create_model()` takes one or more prior-specifications for an
#' SRP multi-state model and combines them into a single model object.
#' Groups are still treated as independent.
#'
#' `define_srp_prior()` specifies a prior distribution for a
#' three state model (stable, response, progression) for
#' a single group.
#'
#' @name srpmodel
#' @aliases srp-model
NULL

#' @param p_mean numeric, mean of the beta prior for the response probability
#' @param p_n numeric, beta prior equivalent sample size (a + b)
#' @param p_eta numeric, robustification parameter for beta prior; actual
#' prior is (1 - eta) beta + eta; i.e., eta is the non-informative weight.
#' @param p_min numeric, minimal response probability
#' @param p_max numeric, maximal response probability
#' @param median_t_q05 numeric of length three,
#' 5% quantiles of the log-normal distributions for the
#' median time-to-next-event for the three transitions s->r, s->p, r->p.
#' @param median_t_q95 numeric of length three,
#' 95% quantiles of the log-normal distributions for the
#' median time-to-next-event for the three transitions s->r, s->p, r->p.
#' @param shape_q05 numeric of length three,
#' 5% quantiles of the log-normal distributions for the shapes of the
#' time-to-next-event distributions for the three transitions s->r, s->p, r->p.
#' @param shape_q95 numeric of length three,
#' 95% quantiles of the log-normal distributions for the shapes of the
#' time-to-next-event distributions for the three transitions s->r, s->p, r->p.
#' @param visit_spacing numeric, fixed duration between visits
#' @param recruitment_rate numeric, constant recruitment rate
#'
#' @examples
#' # a model with prior 25% response rate and variance equivalent to
#' # 10 data points (i.e. a Beta(2.5, 7.5) distribution).
#' grp <- define_srp_prior(p_mean = 0.25, p_n = 10)
#'
#' @rdname srpmodel
#' @export
define_srp_prior <- function(
  p_mean = 0.5,
  p_n = 3,
  p_eta = 0.0,
  p_min = 0.0,
  p_max = 1.0,
  median_t_q05 = c(1, 1, 3),
  median_t_q95 = c(12, 12, 24),
  shape_q05 = rep(0.99, 3),
  shape_q95 = rep(1.01, 3),
  visit_spacing = 1, # months
  recruitment_rate = 1
) {
  checkmate::assert_vector(median_t_q05, len = 3, any.missing = FALSE)
  checkmate::assert_vector(median_t_q95, len = 3, any.missing = FALSE)
  checkmate::assert_vector(shape_q05, len = 3, any.missing = FALSE)
  checkmate::assert_vector(shape_q95, len = 3, any.missing = FALSE)
  if (visit_spacing <= 0) stop("visit spacing must be positive") # nocov
  if (recruitment_rate <= 0) stop("recruitment_rate must be positive") # nocov
  params <- as.list(environment())
  params$visit_spacing <- NULL
  params$recrutiment_rate <- NULL
  res <- structure(
    as.list(environment()),
    visit_spacing = visit_spacing,
    recruitment_rate = recruitment_rate,
    class = "srp_prior"
  )
  return(res)
}



#' @param ... named `srp_prior` objects; the argument names serve as
#' group labels
#' @param maximal_time the maximal overall runtime of the trial as measured from
#' the first visit of any group. No visits past this point are sampled.
#'
#' @examples
#' # a model with two groups and different priors on the respective response
#' # probabilities
#' create_srpmodel(
#'   A = define_srp_prior(),
#'   B = define_srp_prior(p_mean = 0.33, p_n = 10)
#' )
#'
#' @rdname srpmodel
#' @export
create_srpmodel <- function(
  ...,
  maximal_time = 10 * 12
) {
  group_priors <- list(...)
  group_id <- names(group_priors)
  k <- length(group_id)
  if (k == 0) {
    stop("at least one group prior must be specified") # nocov
  }
  if (any(group_id == "")) {
    stop("All arguments passed to ... must be named.") # nocov
  }
  if (length(unique(group_id)) < k) {
    stop("All names of arguments passed to ... must be unique.") # nocov
  }
  transition_labels <- c("s->r", "s->p", "r->p")
  p <- matrix(NA_real_, nrow = k, ncol = 5,
              dimnames = list(group_id, c("mean", "n", "eta", "min", "max")))
  median_t <- array(NA_real_, dim = c(k, 3, 2),
                    dimnames = list(group_id, transition_labels,
                                    c("q05", "q95")))
  shape <- array(NA_real_, dim = c(k, 3, 2),
                 dimnames = list(group_id, transition_labels, c("q05", "q95")))
  visit_spacing <- numeric(k)
  names(visit_spacing) <- group_id
  recruitment_rate <- numeric(k)
  names(recruitment_rate) <- group_id
  for (g in group_id) {
    group <- group_priors[[g]]
    visit_spacing[g] <- attr(group, "visit_spacing")
    recruitment_rate[g] <- attr(group, "recruitment_rate")
    p[g, ] <- c(group$p_mean, group$p_n, group$p_eta, group$p_min, group$p_max)
    median_t[g, , ] <- cbind(group$median_t_q05, group$median_t_q95)
    shape[g, , ] <- cbind(group$shape_mode, group$shape_q05, group$shape_q95)
  }
  res <- structure(
    list(
      group_id = group_id,
      maximal_time = maximal_time,
      visit_spacing = visit_spacing,
      recruitment_rate = recruitment_rate,
      stan_model = stanmodels$srp_model_simple,
      states = c("stable", "response", "progression"),
      prior = list(
        p = p, median_t = median_t, shape = shape
      )
    ),
    class = c("srpmodel", "list"),
    parameter_names = c("p", "median_t", "shape", "scale")
  )
  check_valid(res)
  return(res)
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
  return(arrange(res, t))
}



# generic for digesting data into rstan-ready list
data2standata <- function(data, model) { # nolint
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  # first prepare any data (if available)
  lst_stan_data <- data %>%
    mutate(
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
      t_max = .data$t_max - .data$t_sot
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


# sample from stan model (hidden)
.sample <- function( # nolint
  model,
  data = NULL,
  now = NULL,
  warmup = 250L,
  nsim = 1000L,
  seed = NULL,
  pars = attr(model, "parameter_names"),
  refresh = 0L,
  nuts_control = list(),
  ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (is.null(seed)) # generate seed if none was specified
    seed <- sample.int(.Machine$integer.max, 1)
  if (is.null(data)) {
    data <- .nodata(model)
  } else {
    if (is.null(now)) {
      now <- max(data$t)
    }
    data <- visits_to_mstate(data, model, now)
  }
  # combine prior information with data for stan
  stan_data <- data2standata(data, model)
  # global seed affects permutation of extracted parameters if not set
  set.seed(seed)
  # sample
  res <- rstan::sampling(
    model$stan_model,
    data = stan_data,
    chains = 1L, cores = 1L,
    iter = warmup + nsim, warmup = warmup,
    seed = seed, pars = pars, refresh = refresh,
    control = nuts_control,
    ...
  )
  return(res)
}
