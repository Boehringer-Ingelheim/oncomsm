# Checks for internal consistency of object "Model"
# Returns TRUE if consistent, otherwise throws an error
check_valid <- function(model){ # nolint # nocov start
  checkmate::assert_character(model$group_id)
  group_ids <- model$group_id
  k <- length(group_ids)
  if (model$maximal_time <= 0) stop("maximal time must be positive")
  checkmate::assert_vector(model$visit_spacing, len = k)
  if (any(model$visit_spacing <= 0)) stop("visit spacing must be positive")
  checkmate::assert_vector(model$recruitment_rate, len = k)
  if (any(model$recruitment_rate <= 0)) {
    stop("recruitment rate must be positive")
  }
  checkmate::assert_class(model$stan_model, "stanmodel")
  checkmate::assert_true(all(model$prior$p[, "mean"] > 0))
  checkmate::assert_true(all(model$prior$p[, "n"] > 0))
  checkmate::assert_true(all(model$prior$p[, "min"] < model$prior$p[, "mean"]))
  checkmate::assert_true(all(model$prior$p[, "mean"] < model$prior$p[, "max"]))
  checkmate::assert_true(all(model$prior$median_t[, , "q05"] > 0))
  checkmate::assert_true(all(model$prior$median_t[, , "q05"] <
                               model$prior$median_t[, , "q95"]))
  checkmate::assert_true(all(model$prior$shape[, , "q05"] > 0))
  checkmate::assert_true(all(model$prior$shape[, , "q05"] <
                               model$prior$shape[, , "q95"]))
  return(TRUE)
} # nocov end
