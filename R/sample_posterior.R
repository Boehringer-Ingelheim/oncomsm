#' Sample from posterior distribution of a model
#'
#' @description `sample_posterior()` draws samples from the
#' posterior distribution of the specified model given a data set with
#' visit data
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

#' @template param-warmup
#' @template param-nuts_control
#' @template param-pars
#'
#' @rdname sample_posterior
#' @export
sample_posterior.model <- function(
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
