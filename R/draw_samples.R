#' Draw samples from a model
#'
#' simulate exact response times from the specified model,
#' optionally conditioning on additional data.
#'
#' @param model the model to simulate from
#'
#' @export
draw_samples <- function(model, ...) {
  UseMethod("draw_samples")
}
