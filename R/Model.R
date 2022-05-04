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
draw_samples.OverallModel <- function(model, data = NULL, n = NULL, nsim = 1000L, now = NULL, seed = NULL, ...) {
  if (is.null(now) & is.null(data))
    now <- 0
  if (is.null(now) & !is.null(data))
    now <- max(data$t_recruitment)
  tbl_responses <- draw_samples(
    model$tte_model,
    data = data,
    n = n,
    nsim = nsim,
    seed = seed,
    ...
  )
  tbl_recruitment <- draw_samples(
    model$recruitement_model,
    data = data,
    n = n,
    now = now,
    nsim = nsim,
    seed = seed,
    ...
  )
  dplyr::full_join(
    tbl_responses,
    tbl_recruitment,
    by = c("group_id", "subject_id", "iter")
  )
}


#' @import patchwork ggplot2
#'
#' @export
plot.OverallModel <- function(x, n, nsim = 1e4, ...) {
  p1 <- plot(x$tte_model, n = n, nsim = nsim, ...)
  p2 <- plot(x$recruitement_model, n = n, nsim = nsim, now = 0, ...)
  p1 / p2
}
