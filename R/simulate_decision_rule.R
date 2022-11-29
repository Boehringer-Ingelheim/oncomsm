#' Simulate a custom decision rule
#'
#' @param model model to use for sampling
#' @param decision_rule a function with signature `rule(mdl, data, ...)`
#' returning a
#' data frame with results from a applying the decision rule to data set`data`,
#' typically contains a column `group_id` and a one column per decision/result.
#' @param n_per_group group size
#' @param data a data frame with visit data to condition on
#' @param parameter_sample an optional parameter sample to reuse
#' @param seed optional fixed seed
#' @param nsim the number of resamples to draw from the predictive distribution
#'
#' @return A data frame with columns `iter` (the resample index) and any columns
#' returned by `decision_rule` applied to each of the `nsim` datasets sampled
#' from the predictive distribution.
#'
#' @export
simulate_decision_rule <- function(model,
                                   n_per_group,
                                   decision_rule,
                                   data = NULL,
                                   parameter_sample = NULL,
                                   seed = NULL,
                                   nsim = 1L) {
  # sample parameters
  if (!is.null(seed)) {
    set.seed(seed) # nocov
  }
  smpl <- if (is.null(parameter_sample)) {
    if (is.null(data)) {
      sample_prior(model)
    } else {
      sample_posterior(model, data = data)
    }
  }
  # simulate from predictive distribution
  tbl_complete <- if (is.null(data)) {
    sample_predictive(model, n_per_group = n_per_group, nsim = nsim, sample = smpl)
  } else {
    impute(model, data, n_per_group = n_per_group, nsim = nsim, sample = smpl)
  }
  # apply decision rule
  res <- tbl_complete %>%
    group_by(.data$iter) %>%
    tidyr::nest() %>%
    ungroup() %>%
    transmute(
      iter,
      res = furrr::future_map(
        data, ~decision_rule(model, data = .),
        .options = furrr::furrr_options(seed = TRUE)
      )
    ) %>%
    tidyr::unnest(res)
  return(res)
}