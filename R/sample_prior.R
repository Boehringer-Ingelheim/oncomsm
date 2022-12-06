#' @description `sample_prior()` draws samples from the
#' prior distribution of the specified model object.
#'
#' @examples
#' sample_prior(mdl, seed = 42L)
#'
#' @rdname sample_posterior
#' @export
sample_prior <- function(model,
                         nsim = 2000L,
                         seed = NULL,
                         warmup = 500L,
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
    nuts_control = nuts_control,
    ...
  )
  return(res)
}
