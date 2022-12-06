#' Sample from posterior distribution of a model
#'
#' @description `sample_posterior()` draws samples from the
#' posterior distribution of the specified model given a data set with
#' visit data.
#'
#' @template param-model
#' @template param-data-condition
#' @template param-nsim
#' @template param-seed
#' @template param-warmup
#' @template param-nuts_control
#' @template param-pars
#' @template param-dotdotdot
#'
#' @return A [rstan::stanfit] object with posterior samples.
#'
#' @seealso [rstan::stan()] [parameter_sample_to_tibble()] [sample_prior()]
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
sample_posterior <- function(
  model,
  data,
  nsim = 2000L,
  seed = NULL,
  warmup = 500L,
  nuts_control = list(),
  pars = attr(model, "parameter_names"),
  ...
) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  res <- .sample(
    model, data = data,
    warmup = warmup, nsim = nsim, seed = seed, pars = pars,
    nuts_control = nuts_control, ...
  )
  return(res)
}
