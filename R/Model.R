#' An abstract multi-state model class
#'
#' This is abstract class defining a standard set of methods for any implemented
#' multi-state model.
#' Objects of class 'Model' cannot be instantiated directly, only objects of
#' the respective sub-classes can be.
#'
#' @seealso [srp_model]
#'
#' @name Model
NULL


#' @param x Model to format
#' @template param-dotdotdot
#' @rdname Model
#' @export
format.Model <- function(x, ...) class(x)[1]

#' @param x Model to print
#' @template param-dotdotdot
#' @rdname Model
#' @export
print.Model <- function(x, ...) cat(format(x, ...), "\n")

# Checks for internal consistency of object "Model"
# Returns TRUE if consistent, otherwise throws an error
is_valid <- function(model) {
  UseMethod("is_valid")
}

is_valid.Model <- function(model) {
  stop("Not implemented")
}

#' Sample model prior parameters
#'
#' Sample model parameters from the model prior distribution.
#'
#' @template param-model
#' @template param-warmup
#' @template param-nsim
#' @template param-seed
#' @template param-rstan_output
#' @template param-pars
#' @template param-dotdotdot
#'
#' @return A tibble with columns iter (integer), parameter (character),
#' group_id (character), and value (numeric) with the parameter samples in long
#' format or (if rstan_output == TRUE) a rstanfit object.
#'
#' @export
sample_prior <- function(model, warmup, nsim, seed, rstan_output, pars, ...) {
  UseMethod("sample_prior")
}

#' @inheritParams sample_prior
#' @rdname Model
#' @export
sample_prior.Model <- function(
  model,
  warmup = 500L,
  nsim = 2000L,
  seed = NULL,
  rstan_output = TRUE,
  pars = attr(model, "parameter_names"),
  ...
) {
  res <- .sample(
    model, data = NULL,
    warmup = warmup, nsim = nsim, seed = seed, pars = pars, ...
  )
  if (rstan_output == FALSE) { # convert to tibble representation
    res <- parameter_sample_to_tibble(model, res)
  }
  return(res)
}





#' Sample model posterior parameters
#'
#' Posterior sample model of the parameters conditional on a data set.
#'
#' @template param-model
#' @template param-data-condition
#' @template param-warmup
#' @template param-nsim
#' @template param-seed
#' @template param-rstan_output
#' @template param-pars
#' @template param-dotdotdot
#'
#' @return A tibble with columns iter (integer), parameter (character),
#' group_id (character), and value (numeric) with the parameter samples in long
#' format or (if rstan_output == TRUE) a rstanfit object.
#'
#' @export
sample_posterior <- function(model, data, warmup, nsim, seed, rstan_output,
                             pars, ...) {
  UseMethod("sample_posterior")
}

#' @rdname Model
#' @export
sample_posterior.Model <- function(
  model,
  data,
  warmup = 500L,
  nsim = 2000L,
  seed = NULL,
  rstan_output = TRUE,
  pars = attr(model, "parameter_names"),
  ...
) {
  res <- .sample(
    model, data = data,
    warmup = warmup, nsim = nsim, seed = seed, pars = pars, ...
  )
  if (rstan_output == FALSE) { # convert to tibble representation
    res <- parameter_sample_to_tibble(model, res)
  }
  return(res)
}





#' Sample data from predictive distribution of a model
#'
#' @template param-model
#' @template param-n_per_group
#' @template param-sample
#' @template param-nsim
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return TODO:
#'
#' @export
sample_predictive <- function(model, n_per_group, sample, nsim,
                              nsim_parameters, warmup_parameters, seed, ...) {
  UseMethod("sample_predictive")
}

#'
#' @rdname Model
#' @export
sample_predictive.Model <- function(
  model,
  n_per_group,
  sample = NULL,
  nsim = 100L,
  nsim_parameters = 1000L,
  warmup_parameters = 250,
  seed = NULL,
  ...
) {
  if (is.null(sample)) {
    sample <- sample_prior(
      model, rstan_output = TRUE, seed = seed,
      warmup = warmup_parameters, nsim = nsim_parameters
    )
  }
  .impute(model = model, data = .emptydata(model, n_per_group),
          parameter_sample = sample, now = 0, nsim = nsim, seed = seed, ...)
}



#' Impute data from predictive distribution
#'
#' If no parameter sample is provided, sample from posterior predictive
#'
#' @template param-model
#' @param data the (multi-state) data frame to impute further trajectories for.
#' @param n_per_group the number of individuals per group to be recruited.
#' @param recruitment_rates the per-group recruitment rates.
#' @param now exact time point relative to start of the trial
#' @param sample a stanfit object containing samples. These parameter samples
#'   represent the parameter distribution over which the predictive distribution
#'   averages. Technically, the parameters are resampled with replacement from
#'   this sample to match the desired number of imputations.
#' @template param-nsim
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return a data frame with imputed version of the input data.
#'
#' @export
impute <- function(model, data, n_per_group, recruitment_rates, now,
                   sample, nsim, nsim_parameters, warmup_parameters,
                   seed, ...) {
  UseMethod("impute")
}

#' @inheritParams impute
#' @rdname Model
#' @export
impute.Model <- function(
  model,
  data,
  n_per_group = NULL,
  recruitment_rates = NULL,
  now = NULL,
  sample = NULL,
  nsim = 250L,
  nsim_parameters = 1000L,
  warmup_parameters = 250L,
  seed = NULL,
  ...
) {
  group_ids <- attr(model, "group_id")
  if (is.null(sample)) {
    sample <- sample_posterior(
      model, data = data, rstan_output = TRUE, seed = seed,
      warmup = warmup_parameters, nsim = nsim_parameters
    )
  }
  if (is.null(now)) {
    # convert to visits and take the last time point
    now <- max(mstate_to_visits(model, data)$t)
  }
  if (is.null(n_per_group)) {
    # no new individuals
    n_per_group <- data %>%
      select("group_id", "subject_id") %>%
      distinct() %>%
      pull("group_id") %>%
      table() %>%
      .[group_ids] %>%
      as.numeric()
  } else {
    if (is.null(recruitment_rates)) {
      stop("recruitment_rates must be specified")
    }
  }
  tbl_to_be_recruited <- list()
  for (i in seq_along(group_ids)) {
    n_recruited <- data %>%
      filter(.data$group_id == group_ids[i]) %>%
      pull("subject_id") %>%
      unique() %>%
      length()
    n_to_be_recruited <- n_per_group[i] - n_recruited
    if (n_to_be_recruited < 0) {
      stop("data contains more individuals than specified in n_per_group")
    }
    if (n_to_be_recruited > 0) {
      tbl_to_be_recruited <- rbind(tbl_to_be_recruited, tibble(
          group_id = group_ids[i],
          subject_id = uuid::UUIDgenerate(n = n_to_be_recruited),
          t_sot = now +
            cumsum(stats::rexp(n_to_be_recruited, rate = recruitment_rates[i])),
          from = "stable",
          to = NA_character_,
          t_min = .data$t_sot + 1 / 30,
          t_max = Inf # right censored
        )) %>%
        arrange(.data$t_sot)
    }
  }
  tbl_tmp <- bind_rows(data, tbl_to_be_recruited)
  res <- .impute(model = model, data = tbl_tmp, parameter_sample = sample,
                 nsim = nsim, seed = seed, ...)
  return(res)
}

# must be implemented by "Model" subclass
.impute <- function(model, data, sim, now, parameter_sample, ...) {
  UseMethod(".impute")
}





.sample <- function(model, data, ...) {
  UseMethod(".sample")
}

# sample from stan model (hidden)
.sample.Model <- function( # nolint
  model,
  data = NULL,
  warmup = 250L,
  nsim = 1000L,
  seed = NULL,
  pars = attr(model, "parameter_names"),
  refresh = 0L,
  ...
) {
  if (is.null(seed)) # generate seed if none was specified
    seed <- sample.int(.Machine$integer.max, 1)
  # combine prior information with data for stan
  stan_data <- c(
    as.list(model), # hyperparameters
    data2standata(model, if (is.null(data)) {
      .nodata(model)
      }else {
        data
        }
    ) # data
  )
  # global seed affects permutation of extracted parameters if not set
  set.seed(seed)
  # sample
  res <- rstan::sampling(
    attr(model, "stanmodel"),
    data = stan_data,
    chains = 1L, cores = 1L,
    iter = warmup + nsim, warmup = warmup,
    seed = seed, pars = pars, refresh = refresh, ...
  )
  return(res)
}





#' Convert stanfit sample to data table
#'
#' @template param-model
#' @template param-sample
#' @template param-dotdotdot
#'
#' @return a tibble with the sampled parameters in long format
#'
#' @export
parameter_sample_to_tibble <- function(model, sample, ...) {
  UseMethod("parameter_sample_to_tibble")
}

#' @inheritParams parameter_sample_to_tibble
#' @rdname Model
#' @export
parameter_sample_to_tibble.Model <- function(model, sample, ...) {
  stop("not implemented")
}





# create a data set with no observations
.nodata <- function(model) {
  UseMethod(".nodata")
}

# helper to create empty standata for model
.nodata.Model <- function(model) { # nolint
  stop("not implemented")
}





# create a data set with no observed data
.emptydata <- function(model, n_per_group) {
  UseMethod(".emptydata")
}

# helper to create all-missing standata for model
.emptydata.Model <- function(model, n_per_group) { #nolint
  stop("not implemented")
}





# generic for digesting data into rstan-ready list
data2standata <- function(model, data, ...) {
  UseMethod("data2standata")
}

# convert time to event data to stan data list
data2standata.Model <- function(model, data) {
  stop("not implemented")
}




#' Swimmer-like plot of multi-state data
#'
#' @template param-model
#' @param data a data table with multi-state data
#' @param now the current time relative to the start of the trial (sot)
#' @param relative_to_sot Boolean, should the timeline be relative to the start
#'   of trial or the start of treatment for each individual
#' @template param-dotdotdot
#'
#' @export
plot_mstate <- function(model, data, now, relative_to_sot, ...) {
  UseMethod("plot_mstate")
}

#' @inheritParams plot_mstate
#' @name Model
#' @export
plot_mstate.Model <- function(model, data, now, relative_to_sot, ...) {
  stop("not implemented")
}



#' Generate Visit data from a multi-state model
#'
#' @template param-model
#' @template param-n_per_group
#' @param recruitment_rate numeric vector with the monthly recruitment rates
#' per group
#' @template param-dotdotdot
#'
#' @export
generate_visit_data <- function(model, n_per_group, recruitment_rate, ...) {
  UseMethod("generate_visit_data")
}

#' @inheritParams generate_visit_data
#' @template param-seed
#' @rdname Model
#' @export
generate_visit_data.Model <- function(model, n_per_group, recruitment_rate,
                                      seed = NULL, ...) {
  tbl_data <- sample_predictive(model, n_per_group = n_per_group, nsim = 1,
                                seed = seed) %>% select(-"iter", -"t_sot")
  tbl_recruitment_times <- tbl_data %>%
    select("subject_id", "group_id") %>%
    distinct() %>%
    group_by(.data$group_id) %>%
    mutate( # poisson recruitment process
      rate = purrr::map_dbl(
          .data$group_id,
          ~recruitment_rate[which(attr(model, "group_id") == .)]
        ),
      t_sot = cumsum(stats::rexp(n = n(), rate = .data$rate))
    ) %>%
    ungroup() %>%
    select(-"rate")
  # join and shift transition times accordingly
  res <- left_join(
      tbl_data,
      tbl_recruitment_times,
      by = c("subject_id", "group_id")
    ) %>%
    mutate(
      t_min = .data$t_min + .data$t_sot,
      t_max = .data$t_max + .data$t_sot
    ) %>%
    mstate_to_visits(model, .)
  return(res)
}



#' Sample from the progression-free-survival rate
#'
#' Progression-free-survival rate at time t (PFS-t rate) is a function of the
#' parameters of a given multi-state model. Hence any prior or posterior sample
#' from such a model gives rise to a sample of the corresponding PFS t rate.
#'
#' @template param-model
#' @param t a vector of time-points at which the PFS rate should be computed
#' @template param-sample
#' @template param-warmup
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return a data frame with samples of PFS rates at each of the time points
#' in the vector t.
#'
#' @export
sample_pfs_rate <- function(model, t, sample, warmup, nsim, seed, ...) {
  UseMethod("sample_pfs_rate")
}

#' @inheritParams sample_pfs_rate
#' @rdname Model
#' @export
sample_pfs_rate.Model <- function(
  model,
  t, # PFS_r is 1 - Pr[progression or death before time t]
  sample = NULL,
  warmup = 500L,
  nsim = 2000L,
  seed = NULL,
  ...
) {
  stop("not implemented")
}



#' Convert cross sectional visit data to time-to-event data
#'
#' This function assumes that the visit density is high enough to not miss any
#' transient state jumps.
#'
#' @param tbl_visits visit data in long format
#' @param model a multi-state model object
#' @param now time point since start of trial (might be later than last
#'   recorded visit)
#' @param eof_indicator state name indicating (exactly observed) end of
#'   follow up.
#'
#' @return A data frame
#'
#' @export
visits_to_mstate <- function(tbl_visits, model, now = max(tbl_visits$t),
                             eof_indicator = "EOF") {
  if (!inherits(tbl_visits, "data.frame")) {
    stop("'tbl_visits' must be a data.frame")
  } else {
    assertthat::assert_that(
      inherits(tbl_visits$subject_id, "character"),
      inherits(tbl_visits$group_id, "character"),
      inherits(tbl_visits$t, "numeric"),
      inherits(tbl_visits$state, "character")
    )
  }
  UseMethod("visits_to_mstate", object = model)
}

#' @inheritParams visits_to_mstate
#' @rdname Model
#' @export
visits_to_mstate.Model <- function(tbl_visits, model, now = max(tbl_visits$t),
                                   eof_indicator = "EOF") {
  stop("not implemented")
}
