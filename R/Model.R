#' @param x Model to print
#' @template param-dotdotdot
#' @rdname model
#' @export
print.model <- function(x, ...) cat(format(x, ...), "\n") # nocov

















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
#' @rdname model
#' @export
sample_predictive.model <- function(
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
impute <- function(model, data, nsim, n_per_group, now, seed, ...) {
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
#' @rdname model
#' @export
impute.model <- function(
  model,
  data,
  nsim,
  n_per_group = NULL,
  now = NULL,
  seed = NULL,
  recruitment_rates = model$recruitment_rate,
  sample = NULL,
  nsim_parameters = 1000L,
  warmup_parameters = 250L,
  nuts_control = list(),
  ...
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  group_ids <- model$group_id
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
    if (n_to_be_recruited < 0) { # nocov start
      stop("data contains more individuals than specified in n_per_group") # nolint
    } # nocov end
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
.sample.model <- function( # nolint
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




# create a data set with no observations
.nodata <- function(model) {
  UseMethod(".nodata")
}






# create a data set with no observed data
.emptydata <- function(model, n_per_group, seed) {
  UseMethod(".emptydata")
}




# generic for digesting data into rstan-ready list
data2standata <- function(data, model, ...) {
  UseMethod("data2standata", model)
}




#' Plot the transition time distributions of a model
#'
#' @template param-model
#' @template param-parameter_sample
#' @template param-dotdotdot
#'
#' @export
plot_transition_times <- function(model, parameter_sample, ...) {
  UseMethod("plot_transition_times")
}




#' Plot the response probability distributions of a model
#'
#' @template param-model
#' @template param-parameter_sample
#' @template param-dotdotdot
#'
#' @export
plot_response_probability <- function(model, parameter_sample, ...) {
  UseMethod("plot_response_probability")
}




#' Plot the progression-free-survival function of a model
#'
#' @template param-model
#' @template param-parameter_sample
#' @template param-dotdotdot
#'
#' @export
plot_pfs <- function(model, parameter_sample, ...) {
  UseMethod("plot_pfs")
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



#' Sample from the progression-free-survival rate
#'
#' Progression-free-survival rate at time t (PFS-t rate) is a function of the
#' parameters of a given multi-state model. Hence any prior or posterior sample
#' from such a model gives rise to a sample of the corresponding PFS t rate.
#'
#' @template param-model
#' @param t a vector of time-points at which the PFS rate should be computed
#' @template param-parameter_sample
#' @template param-dotdotdot
#'
#' @return a data frame with samples of PFS rates at each of the time points
#' in the vector t.
#'
#' @export
compute_pfs <- function(model, t, parameter_sample, ...) {
  UseMethod("compute_pfs")
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

