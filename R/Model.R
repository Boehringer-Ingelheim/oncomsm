# standard format/print methods
format.Model <- function(x, ...) class(x)[1]
print.Model <- function(x, ...) cat(format(x, ...), "\n")





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

#' @rdname sample_prior
#' @export
sample_prior.Model <- function(
  model, warmup = 500L, nsim = 2000L, seed = NULL, rstan_output = FALSE,
  pars = attr(model, "parameter_names"), ...
) {
  .sample(
    model, data = NULL,
    warmup = warmup, nsim = nsim, seed = seed, rstan_output = rstan_output, pars = pars, ...
  )
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
sample_posterior <- function(model, data, warmup, nsim, seed, rstan_output, pars, ...) {
  UseMethod("sample_posterior")
}

#' @rdname sample_posterior
#' @export
sample_posterior.Model <- function(
  model, data, warmup = 250L, nsim = 1000L, seed = NULL, rstan_output = FALSE,
  pars = attr(model, "parameter_names"), ...
) {
  .sample(
    model, data = data,
    warmup = warmup, nsim = nsim, seed = seed, rstan_output = rstan_output, pars = pars, ...
  )
}





#' Sample data from prior-predictive distribution
#'
#' Sample time-to-first event data from given model prior-predictive distribution.
#'
#' @template param-model
#' @template param-n_per_arm
#' @template param-nsim
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return  data frame with variables "subject_id", "group_id", "t_recruitment", "dt1" and "dt2"
#' where dt1 is the minimal and dt2 the maximal time to the event in question.
#'
#' @export
sample_prior_predictive <- function(model, n_per_arm, nsim, nsim_parameters, warmup_parameters, seed, ...) {
  UseMethod("sample_prior_predictive")
}

#' @rdname sample_prior_predictive
#' @export
sample_prior_predictive.Model <- function(model, n_per_arm, nsim = 1000L, nsim_parameters = 1000L, warmup_parameters = 250, seed = NULL, ...) {
  prior_sample <- sample_prior(
    model, rstan_output = TRUE, seed = seed,
    warmup = warmup_parameters, nsim = nsim_parameters
  )
  .impute(model = model, data = .emptydata(model, n_per_arm), parameter_sample = prior_sample, now = 0, nsim = nsim, seed = seed, ...)
}





#' Sample data from posterior-predictive distribution
#'
#' Sample time-to-first event data from given model posterior-predictive
#' distribution given data.
#'
#' @template param-model
#' @template param-data-condition
#' @template param-now
#' @template param-nsim
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return  data frame with variables "subject_id", "group_id", "t_recruitment", "dt1", "dt2", "dt_eof"
#' where dt1 is the minimal and dt2 the maximal time to the event in question.
#'
#' @export
impute_posterior_predictive <- function(model, data, now, nsim, nsim_parameters, warmup_parameters, seed, ...) {
  UseMethod("impute_posterior_predictive")
}


#' @rdname impute_posterior_predictive
#' @export
impute_posterior_predictive.Model <- function(model, data, now = NULL, nsim = 1000L, nsim_parameters = 1000L, warmup_parameters = 250L, seed = NULL, ...) {
  posterior_sample <- sample_posterior(
    model, data = data, rstan_output = TRUE, seed = seed,
    warmup = warmup_parameters, nsim = nsim_parameters, ...
  )
  .impute(model = model, data = data, now = now, parameter_sample = posterior_sample, nsim = nsim, seed = seed)
}





.sample <- function(model, data, ...) {
  UseMethod(".sample")
}

# sample from stan model (hidden)
.sample.Model <- function(
  model, data = NULL,
  warmup = 250L, nsim = 1000L, seed = NULL, rstan_output = FALSE,
  pars = attr(model, "parameter_names"), refresh = 0L, ...
) {
  if (is.null(seed)) # generate seed if none was specified
    seed <- sample.int(.Machine$integer.max, 1)
  # combine prior information with data for stan
  stan_data <- c(
    as.list(model), # hyperparameters
    data2standata(model, if (is.null(data)) { .nodata(model) } else { data }) # data
  )
  set.seed(seed) # global seed affects permutation of extracted parameters if not set
  # sample
  res <- rstan::sampling(
    attr(model, "stanmodel"),
    data = stan_data,
    chains = 1L, cores = 1L,
    iter = warmup + nsim, warmup = warmup,
    seed = seed, pars = pars, refresh = refresh, ...
  )
  if (rstan_output == TRUE)
    return(res)
  else
    return(.parameter_sample_to_tibble(model, res))
}





# must be implemented by "Model" subclass
.impute <- function(model, data, sim, now, parameter_sample, ...) {
  UseMethod(".impute")
}





# convert stanfit result to tibble
.parameter_sample_to_tibble <- function(model, sample, ...) {
  UseMethod(".parameter_sample_to_tibble")
}

.parameter_sample_to_tibble.Model <- function(model, sample) {
  stop("not implemented")
}





# create a data set with no observations
.nodata <- function(model) {
  UseMethod(".nodata")
}

# helper to create empty standata for model
.nodata.Model <- function(model) {
  stop("not implemented")
}





# create a data set with no observed data
.emptydata <- function(model, n_per_arm) {
  UseMethod(".emptydata")
}

# helper to create all-missing standata for model
.emptydata.Model <- function(model, n_per_arm) {
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




#' @export
plot_mstate <- function(model, data, ...) {
  UseMethod("plot_mstate")
}

#' @export
plot_mstate.Model <- function(model, data, ...) {
  stop("not implemented")
}



#' @export
generate_visit_data <- function(model, n_per_group, ...) {
  UseMethod("generate_visit_data")
}

#' @export
generate_visit_data.Model <- function(model, n_per_group, seed = NULL, ...) {
  sample_prior_predictive(mdl, n_per_group, 1, seed = seed) %>%
    select(-iter) %>%
    mstate_to_visits(mdl, .)
}
