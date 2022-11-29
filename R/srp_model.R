#' Combine SRP group priors to model object
#'
#' This function takes one or more prior-specifications for an SRP multi-state
#' model and combines them into a joint model.
#' Groups are still treated as independent.
#'
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
#' @return An object of type srp_group_prior holding all prior information in a
#' list-like structure; visit_spacing and recruitment_rate are accessible as
#' attributes.
#'
#' @name srp_model
#' @export
srp_group_prior <- function(
  p_mean = 0.5,
  p_n = 3,
  p_eta = 0.0,
  p_min = 0.0,
  p_max = 1.0,
  median_t_q05 = c(1, 1, 3),
  median_t_q95 = c(12, 12, 24),
  shape_q05 = rep(0.98, 3),
  shape_q95 = rep(1.46, 3),
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
    class = "srp_group_prior"
  )
  return(res)
}



#' Combine SRP group priors to model object
#'
#' This function takes one or more prior-specifications for an SRP multi-state
#' model and combines them into a joint model.
#' Groups are still treated as independent.
#'
#' @param ... named srp_group_prior objects; the argument names serve as
#' group labels
#' @param maximal_time the maximal overall runtime of the trial as measured from
#' the first visit of any group. No visits past this point are sampled.
#'
#' @return an object of class srp_model; a named list of the group_ids,
#' the maximal time, the visit spacing, the recruitment rate, and a list of
#' prior parameters for the response probability, the median transition times,
#' and the shape of the transition time distribution (Weibull)
#'
#' @name srp_model
#' @export
create_srp_model <- function(
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
    class = c("srp_model", "Model", "list"),
    parameter_names = c("p", "median_t", "shape", "scale")
  )
  check_valid(res)
  return(res)
}



#' @param x SRP model to format
#' @template param-dotdotdot
#' @rdname srp_model
#' @export
format.srp_model <- function(x, ...) {
  sprintf("srp_model<%s>", paste(x$group_id, collapse = ","))
}


check_valid.srp_model <- function(model) { # nolint # nocov start
  checkmate::assert_character(model$group_id)
  group_ids <- model$group_id
  k <- length(group_ids)
  if (model$maximal_time <= 0) stop("maximal time must be positive")
  checkmate::assert_vector(model$visit_spacing, len = k)
  if (any(model$visit_spacing <= 0)) stop("visit spacing must be positive")
  checkmate::assert_vector(model$recruitment_rate, len = k)
  if (any(model$recruitment_rate <= 0)) stop("recruitment rate must be positive")
  checkmate::assert_class(model$stan_model, "stanmodel")
  checkmate::assert_true(all(model$prior$p[, "mean"] > 0))
  checkmate::assert_true(all(model$prior$p[, "n"] > 0))
  checkmate::assert_true(all(model$prior$p[, "min"] < model$prior$p[, "mean"]))
  checkmate::assert_true(all(model$prior$p[, "mean"] < model$prior$p[, "max"]))
  checkmate::assert_true(all(model$prior$median_t[, , "q05"] > 0))
  checkmate::assert_true(all(model$prior$median_t[, , "q05"] <
                               model$prior$median_t[, , "q95"]))
  checkmate::assert_true(all(model$prior$shape[, , "q05"] > 0))
  checkmate::assert_true(all(model$prior$shape[, , "q05"] <
                               model$prior$shape[, , "q95"]))
  return(TRUE)
} # nocov end

# see Model.R
.impute.srp_model <- function(model, data, nsim, parameter_sample = NULL, # nolint
                              seed = NULL, p = NULL, shape = NULL,
                              scale = NULL, as_mstate = FALSE, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_groups <- length(model$group_id)
  # make sure that either parameter sample or p, scale, shape are given
  if (!is.null(parameter_sample)) {
    stopifnot(isa(parameter_sample, "stanfit"))
    n_params_sample <- parameter_sample@sim$iter - parameter_sample@sim$warmup
  } else {
    # p, scale, and shape must be given
    if (is.null(p) || is.null(shape) || is.null(scale)) { # nocov start
      stop("if no parameter sample is given all of p, scale, shape must be given") # nolint
    }
    if (length(p) != n_groups) stop()
    if (all(dim(shape) != c(n_groups, 3))) stop()
    if (all(dim(scale) != c(n_groups, 3))) stop()
  } # nocov end
  # extract subject and group id levels for conversion to and back from integer
  subject_id_levels <- unique(as.character(data$subject_id))
  group_id_levels <- model$group_id # important to maintain ordering
  if (is.null(p)) {
    # extract parameter arrays from stanfit object
    # p[i,j] is probability for ith sample for jth group
    p <- rstan::extract(parameter_sample, "p")[[1]]
  } else {
    # expand fixed parameters to same format as rstan parameters
    p <- t(array(p, dim = c(n_groups, n_params_sample)))
  }
  if (is.null(scale)) {
    # scale[i,j,k] is the scale value for ith sample in jth group for
    # k-th transition (k being 1. stable-> response, 2. stable -> progression,
    # 3. response -> progression)
    scale <- rstan::extract(parameter_sample, "scale")[[1]]
  } else {
    # need to add one dimension (iterations to fit format)
    scale <- aperm(
      array(scale, dim = c(n_groups, 3, n_params_sample)),
      c(3, 1, 2)
    )
  }
  if (is.null(shape)) {
    shape <- rstan::extract(parameter_sample, "shape")[[1]]
  } else {
    # need to add one dimension (iterations to fit format)
    shape <- aperm(
      array(shape, dim = c(n_groups, 3, n_params_sample)),
      c(3, 1, 2)
    )
  }
  # sorting the samples and changing type to integer for groups and subj id
  data <- data %>%
    arrange(.data$subject_id, .data$t) %>%
    mutate( # convert group_id to properly ordered factor
      group_id = factor(.data$group_id, levels = group_id_levels)
    )
  idx <- sample(seq_len(n_params_sample), size = nsim, replace = TRUE)
  sample_once <- function(iter) {
    # extract a set of parameters
    response_probabilities <- p[idx[iter], , drop = FALSE]
    shapes <- matrix(shape[idx[iter], , ], ncol = 3)
    scales <- matrix(scale[idx[iter], , ], ncol = 3)
    # sample using C++ implementation
    res <- impute_srp_model(data, response_probabilities, shapes, scales,
        visit_spacing = model$visit_spacing,
        max_time = model$maximal_time
      ) %>%
      as_tibble() %>%
      mutate(
        group_id = as.character(.data$group_id)
      )
    if (as_mstate) {
      res <- visits_to_mstate(res, model)
    }
    res <- mutate(res, iter = as.integer(iter))
    return(res)
  }
  res <- purrr::map_df(1:nsim, sample_once)
  return(res)
}



# helper to create empty standata for model
.nodata.srp_model <- function(model) { # nolint
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



# helper to create all-missing standata for model
.emptydata.srp_model <- function(model, n_per_group, seed = NULL) { # nolint
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



# convert time to event data to stan data list
data2standata.srp_model <- function(data, model) { # nolint
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

#' @importFrom stringr str_extract
#' @export
parameter_sample_to_tibble.srp_model <- function(model, sample, ...) { # nolint
  stopifnot(isa(sample, "stanfit"))
  as.matrix(sample) %>%
    as_tibble() %>%
    mutate(
      iter = row_number()
    ) %>%
    tidyr::pivot_longer(-"iter") %>%
    filter(.data$name != "lp__") %>%
    tidyr::separate("name",
      into = c("parameter", "group_id"),
      sep = "\\[", fill = "right"
    ) %>%
    tidyr::separate("group_id",
      into = c("group_id", "transition"),
      sep = "[\\]|,]", fill = "right", extra = "drop"
    ) %>%
    mutate(
      group_id = model$group_id[as.integer(stringr::str_extract(
        .data$group_id, "[0-9]+"
      ))],
      transition = as.integer(.data$transition)
    )
}



#' @inheritParams visits_to_mstate
#' @name srp_model
#' @export
visits_to_mstate.srp_model <- function(tbl_visits, model, # nolint
                                       now = max(tbl_visits$t),
                                       eof_indicator = "EOF") {
  # make sure everything is sorted
  tbl_visits <- arrange(tbl_visits, .data$subject_id, .data$t)

  tbl_mstate <- list()

  subject_id_lagged <- 0L
  state_lagged <- 0L
  t_sot <- 0 # start of treatment
  for (i in seq_len(nrow(tbl_visits))) {
    if (tbl_visits$subject_id[i] != subject_id_lagged || i == 1) {
      # switch to new subject
      subject_id_lagged <- tbl_visits$subject_id[i]
      state_lagged <- "stable"
      t_sot <- tbl_visits$t[i]
      if (tbl_visits$state[i] != state_lagged) {
        # record jump
        stop(sprintf(
          "first visit must be in starting state; subject_id=%s, state=%s",
          tbl_visits$subject_id[i],
          tbl_visits$state[i]
        ))
      }
    }
    # handle jumps
    if (tbl_visits$state[i] != state_lagged) {
      if (tbl_visits$state[i] == eof_indicator) { # record eof
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
          subject_id = tbl_visits$subject_id[i],
          group_id = tbl_visits$group_id[i],
          from = state_lagged,
          to = NA,
          t_min = tbl_visits$t[i],
          # - Inf indicates censoring and end of follow up
          # (event can no longer be observed)
          t_max = -Inf,
          t_sot = t_sot
        ))
        state_lagged <- tbl_visits$state[i]
      } else { # record jump
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
          subject_id = tbl_visits$subject_id[i],
          group_id = tbl_visits$group_id[i],
          from = state_lagged,
          to = tbl_visits$state[i],
          t_min = tbl_visits$t[i - 1],
          t_max = tbl_visits$t[i],
          t_sot = t_sot
        ))
      }
      # updated current state
      state_lagged <- tbl_visits$state[i]
    }
    # handle non-eof censoring
    censored <- FALSE
    if (!(tbl_visits$state[i] %in% c("progression", eof_indicator))) {
      if (i < nrow(tbl_visits)) {
        if (tbl_visits$subject_id[i] != tbl_visits$subject_id[i + 1]) {
          censored <- TRUE
        }
      } else {
        censored <- TRUE
      }
    }
    if (censored) {
      tbl_mstate <- bind_rows(tbl_mstate, tibble(
        subject_id = tbl_visits$subject_id[i],
        group_id = tbl_visits$group_id[i],
        from = tbl_visits$state[i],
        to = NA,
        t_min = now,
        t_max = Inf,
        t_sot = t_sot
      ))
    }
  }
  return(tbl_mstate)
}



#' @inheritParams compute_pfs
#' @template param-warmup
#'
#' @rdname srp_model
#' @export
compute_pfs.srp_model <- function( # nolint
  model,
  t,
  parameter_sample = NULL,
  warmup = 500L,
  nsim = 1000L,
  seed = NULL,
  ...
) {
  if (is.null(parameter_sample)) {
    parameter_sample <- sample_prior(model,
                                     warmup = warmup, nsim = nsim, seed = seed,
                                     ...)
  }
  group_ids <- model$group_id
  n_groups <- length(group_ids)
  n_smpl <- dim(parameter_sample)[[1]]
  n_t <- length(t)
  res <- array(0.0, dim = c(n_smpl, n_t, n_groups))
  p <- rstan::extract(parameter_sample, "p")[[1]]
  shape <- rstan::extract(parameter_sample, "shape")[[1]]
  scale <- rstan::extract(parameter_sample, "scale")[[1]]
  for (i in 1:n_smpl) {
    for (j in 1:n_groups) {
      # calculate pfs for i-th sample, j-th group and add to result
      res[i, , j] <- pfs(t, p[i, j], shape[i, j, ], scale[i, j, ])
    }
  }
  # unpack 3d array into tibble
  tbl_pfs <- tibble(
    iter = integer(),
    group_id = character(),
    t = numeric(),
    pfs = numeric()
  )
  dimnames(res) <- list(1:n_smpl, 1:n_t, group_ids)
  for (i in seq_along(group_ids)) {
    tbl_pfs <- dplyr::bind_rows(
      tbl_pfs,
      as_tibble(t(res[ , , i])) %>%
        mutate(t = t) %>%
        tidyr::pivot_longer(
          -t, names_to = "iter", values_to = "pfs"
        ) %>%
        mutate(group_id = group_ids[i], iter = as.integer(.data$iter))
    )
  }
  return(tbl_pfs)
}
