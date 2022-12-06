#' Sample from prior distribution of a model
#'
#' @description `sample_prior()` draws samples from the
#' prior distribution of the specified model object.
#'
#' @template param-model
#' @template param-nsim
#' @template param-seed
#' @template param-warmup
#' @template param-pars
#' @template param-nuts_control
#' @template param-dotdotdot
#'
#' @return A [rstan::stanfit] object with sampled prior parameters.
#'
#' @seealso [rstan::stan()] [parameter_sample_to_tibble()] [sample_posterior()]
#' [sample_predictive()] [impute()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' sample_prior(mdl, 1000L, seed = 42L)
#'
#' @export
sample_prior <- function(
  model,
  nsim = 2000L,
  seed = NULL,
  warmup = 500L,
  pars = attr(model, "parameter_names"),
  nuts_control = list(),
  ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  res <- .sample(
    model, data = NULL,
    warmup = warmup, nsim = nsim, seed = seed, pars = pars,
    nuts_control = nuts_control, ...
  )
  return(res)
}
