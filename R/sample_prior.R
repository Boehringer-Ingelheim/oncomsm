#' Sample from prior distribution of a model
#'
#' @description `sample_prior()` draws samples from the
#' prior distribution of the specified model.
#'
#' @template param-model
#' @template param-nsim
#' @template param-seed
#' @template param-dotdotdot
#'
#' @return A rstanfit object with sampled prior parameters.
#'
#' @seealso [parameter_sample_to_tibble()]
#'
#' @export
sample_prior <- function(model, nsim, seed, ...) {
  UseMethod("sample_prior")
}

#' @template param-warmup
#' @template param-pars
#' @template param-nuts_control
#'
#' @examples
#' mdl <- create_srp_model(A = srp_group_prior())
#' sample_prior(mdl, 1000L, seed = 42L)
#'
#' @rdname sample_prior
#' @export
sample_prior.model <- function(
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
