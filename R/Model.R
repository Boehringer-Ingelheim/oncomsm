#' @export
model <- function(tte_model, recruitement_model) {
  assertthat::assert_that(
    inherits(tte_model, "TTEModel"),
    inherits(recruitement_model, "RecruitmentModel")
  )
  res <- list(
    tte_model = tte_model,
    recruitement_model = recruitement_model
  )
  attr(res, "class") <- c("OverallModel", class(res))
  res
}

#' @include draw_samples.R
#'
#' @export
draw_samples.OverallModel <- function(model, data, nsim, now = max(c(0, data$t), na.rm = TRUE), seed = NULL, ...) {
  tbl_responses <- draw_samples(
    model$tte_model,
    data = dplyr::select(data, subject_id, dt1, dt2),
    nsim = nsim,
    seed = seed,
    ...
  )
  tbl_recruitment <- draw_samples(
    model$recruitement_model,
    data = select(data, subject_id, t),
    now = now,
    nsim = nsim,
    seed = seed,
    ...
  )
  dplyr::full_join(
    tbl_responses,
    tbl_recruitment,
    by = c("subject_id", "iter")
  )
}
