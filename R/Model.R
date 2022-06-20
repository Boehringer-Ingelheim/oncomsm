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
  model, warmup = 250L, nsim = 1000L, seed = NULL, rstan_output = FALSE,
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
    warmup = warmup_parameters, nsim = nsim_parameters
  )
  .impute(model = model, data = data, now = now, nsim = nsim, seed = seed, ...)
}





.sample <- function(model, data, ...) {
  UseMethod(".sample")
}

# sample from stan model (hidden)
.sample.Model <- function(
  model, data = NULL,
  warmup = 250L, nsim = 1000L, seed = NULL, rstan_output = FALSE,
  pars = attr(model, "parameter_names"), ...
) {
  if (is.null(seed)) # generate seed if none was specified
    seed <- sample.int(.Machine$integer.max, 1)
  # combine prior information with data for stan
  stan_data <- c(
    as.list(model), # hyperparameters
    data2standata(model, if (is.null(data)) { .nodata(model) } else { data }) # data
  )
  # sample
  res <- rstan::sampling(
    attr(model, "stanmodel"),
    data = stan_data,
    pars = pars,
    chains = 1L, cores = 1L,
    iter = warmup + nsim, warmup = warmup,
    seed = seed, ...
  )
  if (rstan_output == TRUE)
    return(res)
  else
    return(.parameter_sample_to_tibble(model, res))
}





# must be implemented by "Model" subclass
.impute <- function(model, data, ...) {
  UseMethod(".impute")
}





# convert stanfit result to tibble
.parameter_sample_to_tibble <- function(model, sample, ...) {
  UseMethod(".parameter_sample_to_tibble")
}

#' @importFrom stringr str_extract
.parameter_sample_to_tibble.Model <- function(model, sample) {
  stopifnot(isa(sample, "stanfit"))
  as.matrix(sample) %>%
    as_tibble() %>%
    mutate(
      iter = row_number()
    ) %>%
    tidyr::pivot_longer(-.data$iter) %>%
    filter(.data$name != "lp__") %>%
    tidyr::separate(.data$name, into = c("parameter", "group_id"), sep = "\\[", fill = "right") %>%
    mutate(
      group_id = attr(model, "group_id")[as.integer(stringr::str_extract(.data$group_id, "[0-9]+"))]
    )
}





# create a data set with no observations
.nodata <- function(model) {
  UseMethod(".nodata")
}

# helper to create empty standata for model
.nodata.Model <- function(model) {
  tibble(
    subject_id = integer(),
    group_id = integer(),
    t_recruitment = numeric(),
    dt1 = numeric(),
    dt2 = numeric()
  )
}





# create a data set with no observed data
.emptydata <- function(model, n_per_arm) {
  UseMethod(".emptydata")
}

# helper to create all-missing standata for model
.emptydata.Model <- function(model, n_per_arm) {
  n <- sum(n_per_arm)
  tibble(
    subject_id = 1:n,
    group_id = attr(model, "group_id")[rep(1:length(n_per_arm), times = n_per_arm)],
    t_recruitment = rep(NA_real_, n),
    dt1 = rep(NA_real_, n),
    dt2 = rep(NA_real_, n),
    dt_eof = rep(NA_real_, n)
  )
}





# generic for digesting data into rstan-ready list
data2standata <- function(model, data, ...) {
  UseMethod("data2standata")
}





# convert time to event data to stan data list
data2standata.Model <- function(model, data) {
  group_id <- attr(model, "group_id")
  lst_stan_data <- list(
    M_groups = length(group_id)
  )
  # convert group_id to factor (needs to be passed as integer vector)
  data$group_id <- as.integer(factor(data$group_id, levels = group_id))
  # interval censored data
  tmp <- dplyr::filter(data, is.finite(.data$dt1) & is.finite(.data$dt2))
  lst_stan_data$N_A <- nrow(tmp)
  lst_stan_data$group_id_A <- base::as.array(tmp$group_id)
  lst_stan_data$dt1_A <- base::as.array(pmax(1/30, tmp$dt1)) # as.array makes sure we have a vector, even if it is just a single number
  lst_stan_data$dt2_A <- base::as.array(pmax(tmp$dt1 + sqrt(.Machine$double.eps), tmp$dt2))
  # right censored data
  tmp <- dplyr::filter(data, is.finite(.data$dt1) & !is.finite(.data$dt2))
  lst_stan_data$N_B <- nrow(tmp)
  lst_stan_data$group_id_B <- base::as.array(tmp$group_id)
  lst_stan_data$dt1_B <- base::as.array(pmax(1/30, tmp$dt1))
  # definite non-responders
  tmp <- dplyr::filter(data, is.infinite(.data$dt1))
  lst_stan_data$N_C <- nrow(tmp)
  lst_stan_data$group_id_C <- base::as.array(tmp$group_id)
  # new individuals
  tmp <- dplyr::filter(data, is.na(.data$dt1) & is.na(.data$dt2))
  lst_stan_data$N_D <- nrow(tmp)
  lst_stan_data$group_id_D <- base::as.array(tmp$group_id)
  # construct waiting times for recruitment
  dt_recruitment <- with(lst_stan_data,
    matrix(0, nrow = N_A + N_B + N_C, ncol = M_groups)
  )
  n_recruited_per_group <- integer(lst_stan_data$M_groups)
  for (i in 1:lst_stan_data$M_groups) {
    t_recruitment <- data %>%
      filter(group_id == i, !is.na(t_recruitment)) %>%
      pull(t_recruitment)
    # set first waiting time to small positive number (referencing to first patient in)
    # enforce minimal wait between individuals to avoid numerical issues
    n_recruited_per_group[i] <- length(t_recruitment)
    if (n_recruited_per_group[i] > 0) {
      dt_recruitment[1:(n_recruited_per_group[i]), i] <- c(
        sqrt(.Machine$double.eps),
        pmax(sqrt(.Machine$double.eps), diff(sort(t_recruitment)))
      )
    }
  }
  lst_stan_data$dt_recruitment <- base::as.array(dt_recruitment)
  lst_stan_data$n_recruited_per_group <- base::as.array(n_recruited_per_group)
  # return
  return(lst_stan_data)
}
