#' Sample visits from predictive distribution
#'
#' `sample_predictive()` draws samples from the predictive distribution of a
#' model given a parameter sample.
#'
#' @template param-model
#' @template param-n_per_group
#' @template param-sample
#' @template param-nsim
#' @template param-seed
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @param as_mstate return data in multi-state forma, see [visits_to_mstate()]
#' @template param-nuts_control
#' @template param-dotdotdot
#'
#' @return a data frame with variables
#' `subject_id<chr>` (subject identifier),
#' `group_id<chr>` (group identifier),
#' `t<dbl>` (time of visit, relative to first visit in study),
#' `state<chr>` (state recorded at visit)
#' `iter<int>` (re-sample indicator).
#' Allowed states are "stable", "response", "progression" (or death),
#' and "EOF" (end of follow-up).
#' The EOF state marks the end of an individual's follow-up before the absorbing
#' state "progression".
#'
#' @seealso [sample_prior()] [sample_posterior()] [impute()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' sample_predictive(mdl, 1L, nsim = 1L, seed = 38L)
#'
#' @export
sample_predictive <- function(
  model,
  n_per_group,
  sample = NULL,
  p = NULL,
  shape = NULL,
  scale = NULL,
  nsim = 100L,
  seed = NULL,
  nsim_parameters = 1000L,
  warmup_parameters = 250,
  as_mstate = FALSE,
  nuts_control = list(),
  ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(sample)) {
    sample <- sample_prior(model,
                           warmup = warmup_parameters, nsim = nsim_parameters,
                           nuts_control = nuts_control)
  }
  # construct an empty data set
  data <- .emptydata(model, n_per_group)
  # call model-specific imputation method
  .impute(model = model, data = data, parameter_sample = sample, now = 0,
          nsim = nsim, as_mstate = as_mstate, shape = shape, scale = scale, p = p, ...)
}
