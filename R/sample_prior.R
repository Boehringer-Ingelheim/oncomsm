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
    data = tibble::tibble(
        group_id = character(0L),
        subject_id = character(0L),
        t = numeric(0L),
        state = character(0L)
      ),
    now = 0,
    nsim = nsim,
    seed = seed,
    warmup = warmup,
    nuts_control = nuts_control,
    ...
  )
  return(res)
}
