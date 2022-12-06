#' @description `sample_prior()` draws samples from the
#' prior distribution of the specified model object.
#'
#' @examples
#' sample_prior(mdl, 1000L, seed = 42L)
#'
#' @rdname sample_posterior
#' @export
sample_prior <- function(model,
                         nsim = 1000L,
                         seed = NULL,
                         warmup = 500L,
                         pars = attr(model, "parameter_names"),
                         nuts_control = list(),
                         ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  # just call sample_posterior without data
  res <- sample_posterior(
    model = model,
    data = NULL,
    now = 0,
    nsim = nsim,
    seed = seed,
    warmup = warmup,
    pars = pars,
    nuts_control = nuts_control,
    ...
  )
  return(res)
}
