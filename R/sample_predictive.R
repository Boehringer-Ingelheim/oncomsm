#' Sample visits from predictive distribution
#'
#' `sample_predictive()` draws samples from the predictive distribution of a
#' model given a parameter sample.
#'
#' @template param-model
#' @template param-n_per_group
#' @template param-sample
#' @param p numeric, vector of optional fixed response probabilities to use for
#' sampling
#' @param scale numeric, matrix of optional fixed Weibull scale parameters to
#' use for sampling must be a matrix of dim c(n_groups, 3) where the second
#' dimension corresponds to the transitions between s->r, s->p, r->p
#' @param shape numeric, matrix of optional fixed Weibull shape parameters to
#' use for sampling must be a matrix of dim c(n_groups, 3) where the second
#' dimension corresponds to the transitions between s->r, s->p, r->p
#' @template param-nsim
#' @template param-seed
#' @template param-nsim_parameters
#' @template param-warmup_parameters
#' @template param-nuts_control
#' @param as_mstate logical, return data in mstate format?
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
#' @seealso [sample_prior()] [sample_posterior()]
#'
#' @examples
#' sample_predictive(mdl, 1L, 20L, seed = 38L)
#'
#' @aliases impute
#' @export
sample_predictive <- function(model,
                              nsim,
                              n_per_group,
                              sample = NULL,
                              p = NULL,
                              shape = NULL,
                              scale = NULL,
                              seed = NULL,
                              nsim_parameters = 1000L,
                              warmup_parameters = 250,
                              nuts_control = list(),
                              as_mstate = FALSE,
                              ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  # just call impute with empty data set
  impute(
    model,
    data = .emptydata(model, n_per_group), # construct an empty data set
    nsim = nsim,
    n_per_group = n_per_group,
    now = NULL,
    seed = seed,
    recruitment_rate = model$recruitment_rate,
    p = p,
    shape = shape,
    scale = scale,
    sample = sample,
    nsim_parameters = nsim_parameters,
    warmup_parameters = warmup_parameters,
    nuts_control = nuts_control,
    as_mstate = as_mstate,
    ...
  )
}
