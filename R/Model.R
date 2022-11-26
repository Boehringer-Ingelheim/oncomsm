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
  stop("Not implemented") # nocov
}

#' Sample model prior parameters
#'
#' Sample model parameters from the model prior distribution.
#'
#' @template param-model
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return a rstanfit object with the sampled prior parameters
#'
#' @seealso [parameter_sample_to_tibble()]
#'
#' @export
sample_prior <- function(model, nsim, seed, ...) {
  UseMethod("sample_prior")
}

#' @inheritParams sample_prior
#' @template param-warmup
#' @template param-pars
#' @template param-nuts_control
#'
#' @rdname Model
#' @export
sample_prior.Model <- function(
  model,
  nsim = 2000L,
  seed = NULL,
  warmup = 500L,
  pars = attr(model, "parameter_names"),
  nuts_control = list(),
  ...
) {
  res <- .sample(
    model, data = NULL,
    warmup = warmup, nsim = nsim, seed = seed, pars = pars,
    nuts_control = nuts_control, ...
  )
  return(res)
}





#' Sample model posterior parameters
#'
#' Posterior sample model of the parameters conditional on a data set.
#'
#' @template param-model
#' @template param-data-condition
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return A rstanfit object with posterior samples.
#'
#' @seealso [parameter_sample_to_tibble()]
#'
#' @export
sample_posterior <- function(model, data, nsim, seed, ...) {
  UseMethod("sample_posterior")
}

#' @inheritParams sample_posterior
#' @template param-warmup
#' @template param-nuts_control
#' @template param-pars
#'
#' @rdname Model
#' @export
sample_posterior.Model <- function(
  model,
  data,
  nsim = 2000L,
  seed = NULL,
  warmup = 500L,
  nuts_control = list(),
  pars = attr(model, "parameter_names"),
  ...
) {
  res <- .sample(
    model, data = data,
    warmup = warmup, nsim = nsim, seed = seed, pars = pars,
    nuts_control = nuts_control, ...
  )
  return(res)
}





#' Sample data from predictive distribution of a model
#'
#' @template param-model
#' @template param-n_per_group
#' @template param-sample
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return TODO:
#'
#' @export
sample_predictive <- function(model, n_per_group, sample, nsim, seed, ...) {
  UseMethod("sample_predictive")
}

#' @inheritParams sample_predictive
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @param as_mstate return data in multi-state forma, see [visits_to_mstate()]
#' @template param-nuts_control
#'
#' @rdname Model
#' @export
sample_predictive.Model <- function(
  model,
  n_per_group,
  sample = NULL,
  nsim = 100L,
  seed = NULL,
  nsim_parameters = 1000L,
  warmup_parameters = 250,
  as_mstate = FALSE,
  nuts_control = list(),
  ...
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(sample)) {
    sample <- sample_prior(model,
                           warmup = warmup_parameters, nsim = nsim_parameters,
                           nuts_control = nuts_control)
  }
  # construct an empty data set
  data <- .emptydata(model, n_per_group)
  # call model-specific imputation method
  .impute(model = model, data = data, parameter_sample = sample, now = 0,
          nsim = nsim, as_mstate = as_mstate, ...)
}



#' Impute data from predictive distribution
#'
#' If no parameter sample is provided, sample from posterior predictive
#'
#' @template param-model
#' @param data the (multi-state) data frame to impute further trajectories for.
#' @param n_per_group the number of individuals per group to be recruited.
#' @param now exact time point relative to start of the trial
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return a data frame with imputed version of the input data.
#'
#' @export
impute <- function(model, data, nsim, n_per_group, now, seed, ... ) {
  UseMethod("impute")
}

#' @inheritParams impute
#' @param recruitment_rates vector of recruitment rates
#' @param sample a stanfit object containing samples. These parameter samples
#'   represent the parameter distribution over which the predictive distribution
#'   averages. Technically, the parameters are resampled with replacement from
#'   this sample to match the desired number of imputations.
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @template param-nuts_control
#'
#' @rdname Model
#' @export
impute.Model <- function(
  model,
  data,
  nsim,
  n_per_group = NULL,
  now = NULL,
  seed = NULL,
  recruitment_rates = attr(model, "recruitment_rate"),
  sample = NULL,
  nsim_parameters = 1000L,
  warmup_parameters = 250L,
  nuts_control = list(),
  ...
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  group_ids <- attr(model, "group_id")
  if (is.null(sample)) {
    sample <- sample_posterior(model,
                               data = data, seed = seed,
                               warmup = warmup_parameters,
                               nsim = nsim_parameters,
                               nuts_control = nuts_control)
  }
  if (is.null(now)) { # convert to visits and take the last time point
    now <- max(data$t)
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
      stop("recruitment_rates must be specified") # nocov
    }
  }
  tbl_to_recruit <- tibble(
    subject_id = character(0L),
    group_id = character(0L),
    t = numeric(0L),
    state = character(0L)
  )
  for (i in seq_along(group_ids)) {
    n_recruited <- data %>%
      filter(.data$group_id == group_ids[i]) %>%
      pull("subject_id") %>%
      unique() %>%
      length()
    n_to_be_recruited <- n_per_group[i] - n_recruited
    if (n_to_be_recruited < 0) {
      stop("data contains more individuals than specified in n_per_group") # nocov nolint
    }
    ids_to_exclude <- c(tbl_to_recruit$subject_id, unique(data$subject_id))
    if (n_to_be_recruited > 0) {
      subject_ids <- get_identifier(n = n_to_be_recruited,
                                    exclude = ids_to_exclude)
      recruitment_times <- now +
        cumsum(stats::rexp(n_to_be_recruited, rate = recruitment_rates[i]))
      tbl_to_recruit <- bind_rows(tbl_to_recruit, tibble(
        subject_id = subject_ids,
        group_id = group_ids[i],
        t = recruitment_times,
        state = "stable" # first visits are always stable
      ))
    }
  }
  tbl_data <- bind_rows(data, tbl_to_recruit)
  res <- .impute(model = model, data = tbl_data, parameter_sample = sample,
                 nsim = nsim, seed = seed, ...)
  return(res) # problem with forward sampling?
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
  now = NULL,
  warmup = 250L,
  nsim = 1000L,
  seed = NULL,
  pars = attr(model, "parameter_names"),
  refresh = 0L,
  nuts_control = list(),
  ...
) {
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
  stan_data <- c(as.list(model), data2standata(data, model))
  # global seed affects permutation of extracted parameters if not set
  set.seed(seed)
  # sample
  res <- rstan::sampling(
    attr(model, "stanmodel"),
    data = stan_data,
    chains = 1L, cores = 1L,
    iter = warmup + nsim, warmup = warmup,
    seed = seed, pars = pars, refresh = refresh,
    control = nuts_control,
    ...
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
  stop("not implemented") # nocov
}





# create a data set with no observations
.nodata <- function(model) {
  UseMethod(".nodata")
}

# helper to create empty standata for model
.nodata.Model <- function(model) { # nolint
  stop("not implemented") # nocov
}





# create a data set with no observed data
.emptydata <- function(model, n_per_group, seed) {
  UseMethod(".emptydata")
}

# helper to create all-missing standata for model
.emptydata.Model <- function(model, n_per_group, seed = NULL) { #nolint
  stop("not implemented") # nocov
}





# generic for digesting data into rstan-ready list
data2standata <- function(data, model, ...) {
  UseMethod("data2standata", model)
}

# convert time to event data to stan data list
data2standata.Model <- function(data, model) {
  stop("not implemented") # nocov
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
plot_mstate <- function(data, model, now, relative_to_sot, ...) {
  UseMethod("plot_mstate", object = model)
}

#' @inheritParams plot_mstate
#' @name Model
#' @export
plot_mstate.Model <- function(data, model, now, relative_to_sot, ...) {
  stop("not implemented") # nocov
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
  stop("not implemented") # nocov
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
    stop("'tbl_visits' must be a data.frame") # nocov
  } else {
    checkmate::test_true(inherits(tbl_visits$subject_id, "character"))
    checkmate::test_true(inherits(tbl_visits$group_id, "character"))
    checkmate::test_true(inherits(tbl_visits$t, "numeric"))
    checkmate::test_true(inherits(tbl_visits$state, "character"))
  }
  UseMethod("visits_to_mstate", object = model)
}

#' @inheritParams visits_to_mstate
#' @rdname Model
#' @export
visits_to_mstate.Model <- function(tbl_visits, model, now = max(tbl_visits$t),
                                   eof_indicator = "EOF") {
  stop("not implemented") # nocov
}
