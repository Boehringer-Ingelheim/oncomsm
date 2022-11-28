#' @name srp_model
#' @export
srp_group_prior <- function(
  p_mean = 0.5,
  p_n = 1,
  p_eta = 0.0,
  p_min = 0.0,
  p_max = 1.0,
  median_t_mode = c(3, 3, 6),
  median_t_90q = c(30, 30, 60),
  median_t_max = c(Inf, Inf, Inf),
  shape_mode = c(1, 1, 1),
  shape_90q = c(5, 5, 5),
  shape_min = c(0.01, 0.01, 0.01),
  shape_max = c(20, 20, 20), # about .999 quantile
  visit_spacing = 1, # months
  recruitment_rate = 1
) {
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

#'
#' @export
srp_modell <- function(
  ...,
  maximal_time = 10 * 12
) {
  group_priors <- list(...)
  group_id <- names(group_priors)
  k <- length(group_id)
  if (any(group_id == "")) {
    stop("All arguments passed to ... must be named.") # nocov
  }
  if (length(unique(group_id)) < k) {
    stop("All names of arguments passed to ... must be unique.") # nocov
  }
  transition_labels <- c("s->r", "s->p", "r->p")
  p <- matrix(NA_real_, nrow = k, ncol = 5,
              dimnames = list(group_id, c("mean", "n", "eta", "min", "max")))
  median_t <- array(NA_real_, dim = c(k, 3, 3),
                    dimnames = list(group_id, transition_labels,
                                    c("mode", "90q", "max")))
  shape <- array(NA_real_, dim = c(k, 3, 4),
                 dimnames = list(group_id, transition_labels,
                                 c("mode", "90q", "min", "max")))
  visit_spacing <- numeric(k)
  names(visit_spacing) <- group_id
  recruitment_rate <- numeric(k)
  names(recruitment_rate) <- group_id
  for (g in group_id) {
    group <- group_priors[[g]]
    visit_spacing[g] <- attr(group, "visit_spacing")
    recruitment_rate[g] <- attr(group, "recruitment_rate")
    p[g, ] <- c(group$p_mean, group$p_n, group$p_eta, group$p_min, group$p_max)
    median_t[g, , ] <- cbind(group$median_t_mode, group$median_t_90q,
                             group$median_t_max)
    shape[g, , ] <- cbind(group$shape_mode, group$shape_90q, group$shape_min,
                          group$shape_max)
  }
  res <- structure(
    list(
      group_id = group_id,
      maximal_time = maximal_time,
      visit_spacing = visit_spacing,
      recruitment_rate = recruitment_rate,
      prior = list(
        p = p, median_t = median_t, shape = shape
      )
    ),
    class = "srp_model"
  )
  return(res)
}



#' A Stable-Response-Progression Model
#'
#' Create a new instance of an SRP model
#'
#' TODO
#'
#' @param group_id a character vector with the group ids, these are used to
#'   check compatibility of data later
#' @param logodds_mean a vector with the means of the (truncated) normal priors
#'   on the log-odds of the response probability
#' @param logodds_sd a vector with the standard deviations of the (truncated)
#'   normal priors on the log-odds of the response probability
#' @param median_time_to_next_event_mean a matrix with the means of the
#'   (truncated) normal priors on the median time to next event for each of the
#'   Weibull transition probabilities, the (i,j)-th entry is the i-th group
#'   median time to next event for transition j
#'   (1=stable-response, 2=stable-progression, 3=response-progression)
#' @param median_time_to_next_event_sd a matrix with the standard deviations of
#'   the (truncated) normal priors on the median time to next event for each of
#'   the Weibull transition probabilities, the (i,j)-th entry is the i-th group
#'   median time to next event for transition j
#'   (1=stable-response, 2=stable-progression, 3=response-progression)
#' @param visit_spacing vector of time differences between visits per group,
#' only relevant for sampling from the predictive distribution
#' @param recruitment_rate vector with per-group recruitment rates,
#' only relevant for sampling from the predictive distribution
#' @param max_time maximal overall runtime from first visit,
#' only relevant for sampling from the predictive distribution
#' @param logodds_min lower boundary on the log-odds per group
#' @param logodds_max upper boundary on the log-odds per group
#' @param shape_min matrix of lower boundaries of the uniform prior of the
#'   Weibull distribution per group/transition
#' @param shape_max matrix of upper boundaries of the uniform prior of the
#'   Weibull distribution per group/transition
#'
#' @name srp_model
#' @aliases create_srp_model
#' @seealso [Model]
#'
#' @export
create_srp_model <- function(
  group_id,
  logodds_mean,
  logodds_sd = rep(0.01, length(group_id)),
  median_time_to_next_event_mean,
  median_time_to_next_event_sd = matrix(0.1, nrow = length(group_id), ncol = 3),
  visit_spacing,
  recruitment_rate = rep(1, length(group_id)),
  max_time = 10 * 12,
  logodds_min = rep(logodds(.001), length(group_id)),
  logodds_max = rep(logodds(.999), length(group_id)),
  shape_min = matrix(.001, nrow = length(group_id), ncol = 3),
  shape_max = matrix(20, nrow = length(group_id), ncol = 3), # about .999 quantile
  shape_mode = matrix(1, nrow = length(group_id), ncol = 3),
  shape_90q = matrix(5, nrow = length(group_id), ncol = 3),
  median_t_mode = matrix(6, nrow = length(group_id), ncol = 3),
  median_t_90q = matrix(12, nrow = length(group_id), ncol = 3),
  median_t_max = matrix(Inf, nrow = length(group_id), ncol = 3),
  p_mean = rep(0.5, length(group_id)),
  p_n = rep(1, length(group_id)),
  p_eta = rep(0.0, length(group_id)),
  p_min = rep(0.0, length(group_id)),
  p_max = rep(1.0, length(group_id))
) {
  # calculate mu, sigma for log-normal priors
  shape_mu <- matrix(NA_real_, nrow = length(group_id), ncol = 3)
  shape_sigma <- matrix(NA_real_, nrow = length(group_id), ncol = 3)
  median_t_mu <- matrix(NA_real_, nrow = length(group_id), ncol = 3)
  median_t_sigma <- matrix(NA_real_, nrow = length(group_id), ncol = 3)
  for (i in 1:3) {
    for (j in 1:length(group_id)) {
      tmp <- get_mu_sigma(shape_mode[j, i], shape_90q[j, i])
      shape_mu[j, i] <- tmp$mu
      shape_sigma[j, i] <- tmp$sigma
      tmp <- get_mu_sigma(median_t_mode[j, i], median_t_90q[j, i])
      median_t_mu[j, i] <- tmp$mu
      median_t_sigma[j, i] <- tmp$sigma
    }
  }
  mdl <- list(
    logodds_mean = logodds_mean,
    logodds_sd = logodds_sd,
    median_time_to_next_event_mean = median_time_to_next_event_mean,
    median_time_to_next_event_sd = median_time_to_next_event_sd,
    logodds_min = logodds_min,
    logodds_max = logodds_max,
    shape_min = shape_min,
    shape_max = shape_max,
    shape_mu = shape_mu,
    shape_sigma = shape_sigma,
    median_t_mu = median_t_mu,
    median_t_sigma = median_t_sigma,
    median_t_max = median_t_max,
    p_mean = p_mean,
    p_n = p_n,
    p_eta = p_eta,
    p_min = p_min,
    p_max = p_max
  )
  mdl <- lapply(mdl, base::as.array)
  attr(mdl, "group_id") <- as.character(group_id) # assert type
  attr(mdl, "states") <- c("stable", "response", "progression")
  attr(mdl, "visit_spacing") <- as.array(visit_spacing)
  attr(mdl, "recruitment_rate") <- as.array(recruitment_rate)
  attr(mdl, "max_time") <- as.array(max_time)
  attr(mdl, "stanmodel") <- stanmodels[["srp_model"]]
  attr(mdl, "parameter_names") <- c("p", "shape", "scale",
                                    "median_time_to_next_event")
  class(mdl) <- c("srp_model", "Model", class(mdl))
  is_valid(mdl)
  return(mdl)
}


#' @param x SRP model to format
#' @template param-dotdotdot
#' @rdname srp_model
#' @export
format.srp_model <- function(x, ...) {
  sprintf("srp_model<%s>", paste(attr(x, "group_id"), collapse = ","))
}


is_valid.srp_model <- function(mdl) { # nolint
  with(mdl, {
    checkmate::assert_vector(logodds_mean,
      len = length(attr(mdl, "group_id")),
      any.missing = FALSE, .var.name = "logodds_mean"
    )
    checkmate::assert_vector(logodds_sd,
      len = length(attr(mdl, "group_id")),
      any.missing = FALSE, .var.name = "logodds_mean"
    )
    checkmate::assert_array(median_time_to_next_event_mean,
      d = 2,
      any.missing = FALSE,
      .var.name = "median_time_to_next_event_mean"
    )
    checkmate::assert_array(median_time_to_next_event_sd,
      d = 2,
      any.missing = FALSE,
      .var.name = "median_time_to_next_event_sd"
    )
    checkmate::assert_vector(median_time_to_next_event_sd,
      len = length(attr(mdl, "group_id")) *
        length(attr(mdl, "states")), any.missing = FALSE,
      .var.name = "median_time_to_next_event_sd"
    )
    checkmate::assert_vector(median_time_to_next_event_mean,
      len = length(attr(mdl, "group_id")) *
        length(attr(mdl, "states")), any.missing = FALSE,
      .var.name = "median_time_to_next_event_mean"
    )
    checkmate::assertTRUE(all(logodds_mean < logodds_max),
      .var.name = "logodds_mean < logodds_max"
    )
    with(
      mdl,
      checkmate::assert_numeric(
        median_time_to_next_event_mean[median_time_to_next_event_mean < 0],
        lower = 0, upper = 0,
        .var.name = "median_time_to_next_event_mean"
      )
    )
  })
  return(TRUE)
}

# see Model.R
.impute.srp_model <- function(model, data, nsim, parameter_sample = NULL, # nolint
                              seed = NULL, p = NULL, shape = NULL,
                              scale = NULL, as_mstate = FALSE, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_groups <- length(attr(model, "group_id"))
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
  group_id_levels <- attr(model, "group_id") # important to maintain ordering
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
        visit_spacing = attr(model, "visit_spacing"),
        max_time = attr(model, "max_time")
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
  group_ids <- attr(model, "group_id")
  rr <- attr(model, "recruitment_rate")
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
  lst_stan_data <- data %>%
    mutate(
      group_id = as.integer(factor(.data$group_id,
        levels = attr(model, "group_id")
      )),
      subject_id = as.integer(factor(.data$subject_id)),
      from = as.integer(factor(.data$from, levels = attr(model, "states"))),
      to = if_else(
        is.na(.data$to),
        "unknown",
        as.character(.data$to)
      ) %>%
      factor(levels = c(attr(model, "states"), "unknown")) %>%
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
  lst_stan_data$M_groups <- length(attr(model, "group_id"))
  lst_stan_data$N <- nrow(data)
  lst_stan_data$N_subjects <- length(unique(data$subject_id))
  lst_stan_data <- c(lst_stan_data, as.list(model))
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
      group_id = attr(model, "group_id")[as.integer(stringr::str_extract(
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
                                     pars = attr(model, "parameter_names"), ...)
  }
  group_ids <- attr(model, "group_id")
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
