#' Sample parameters from a model
#'
#' @description `sample_posterior()` draws samples from the
#' posterior distribution of the specified model given a data set with
#' visit data.
#'
#' @template param-model
#' @template param-data-condition
#' @param now numeric, time from first visit in data if different form last
#' recorded visit
#' @template param-nsim
#' @template param-seed
#' @template param-warmup
#' @template param-nuts_control
#' @template param-pars
#' @template param-dotdotdot
#'
#' @return A [rstan::stanfit] object with posterior samples.
#'
#' @seealso [rstan::stan()] [parameter_sample_to_tibble()]
#' [sample_predictive()] [impute()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' tbl <- tibble::tibble(
#'   subject_id = c("A1", "A1"),
#'   group_id = c("A", "A"),
#'   t = c(0, 1.5),
#'   state = c("stable", "response")
#' )
#' sample_posterior(mdl, tbl, 500L, 42L)
#'
#' @export
sample_posterior <- function(model,
                             data,
                             now = NULL,
                             nsim = 1000L,
                             seed = NULL,
                             warmup = 500L,
                             nuts_control = list(),
                             pars = attr(model, "parameter_names"),
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
    seed = seed, pars = pars, refresh = 0L,
    control = nuts_control,
    ...
  )
  return(res)
}
