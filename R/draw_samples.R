#' Draw samples from a model
#'
#' simulate exact response times from the specified model,
#' optionally conditioning on additional data.
#'
#' @param model the model to simulate from
#' @param n_to_be_recruited int, number of individuals to sample for; these are
#' future individuals and additional the (optional) individuals provided via the "data" argument.
#' @param nsim int >= 10L, number of samples to draw
#' @param data **optional** data to condition on
#' @param seed **optional** random seed for simulation
#'
#' @export
draw_samples <- function(model, n_to_be_recruited, nsim, now = 0, data = NULL, seed = NULL, ...) {
  UseMethod("draw_samples")
}
