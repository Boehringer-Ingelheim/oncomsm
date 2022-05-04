#' Draw samples from a model
#'
#' simulate exact response times from the specified model,
#' optionally conditioning on additional data.
#'
#' @param model the model to simulate from
#' @param data NULL (then n must be given) or data frame with columns "group_id", "subject_id",
#' "t_recruitment" (recruitment time since start of trial),
#' "dt1" (minimal time since recruitment to response, Inf if non-responder),
#' "dt2" (maximal time since recruitment to response, Inf if non-responder, NA if censored)
#' @param n NULL (then data must be provided) or an integer vector with the
#' per-arm sample sizes for each group_id as given in "model"
#'
#' @export
draw_samples <- function(model, data, n, ...) {
  UseMethod("draw_samples")
}
