#' Create an instance of the mixture cure-rate model for binary tumor response
#'
#' This model assumes that the response rates across groups are independent.
#'
#' @param group_id character vector of group ids
#' @param logodds_mean numeric vector with the means of the normal prior on the logodds of a response per group
#' @param logodds_sd numeric vector with the standard deviations of the normal prior on the logodds of a response per group
#' @param logodds_min numeric vector with the minimal logodds per group
#' @param logodds_max numeric vector with the maximal logodds per group
#' @param log_shape_mean numeric vector with the means of the normal prior on the log-shape parameter of the Weibull distribution for time to response
#' @param log_shape_sd numeric vector with the standard deviations of the normal prior on the log-shape parameter of the Weibull distribution for time to response
#' @param median_time_to_response_mean numeric vector with the means of the normal prior on the median time to response (truncated at 0)
#' @param median_time_to_response_sd numeric vector with the standard deviations of the normal prior on the median for time to response
#' @param max_time_to_response numeric vector with the maximal time to response per group
#' @param visit_spacing vector of deterministic spacing between future visits (in months)
#' @param monthly_rate_mean expected monthly recruitment rate (prior mean)
#' @param monthly_rate_sd standard deviation of the expected recruitment rate prior
#'
#' @return An object of class "IndependentMixtureCureRateModel" holding all relevant
#' prior information.
#'
#' @rdname IndependentMixtureCureRateModel
#'
#' @export
new_IndependentMixtureCureRateModel <- function(
  group_id,
  logodds_mean, logodds_sd,
  logodds_min = rep(logodds(.001), length(group_id)),
  logodds_max = rep(logodds(.999), length(group_id)),
  log_shape_mean, log_shape_sd,
  median_time_to_response_mean, median_time_to_response_sd,
  max_time_to_response,
  visit_spacing,
  monthly_rate_mean,
  monthly_rate_sd
) {
  mdl <- as.list(environment()) # store all input parameters
  mdl$group_id <- NULL
  mdl <- lapply(mdl, base::as.array)
  attr(mdl, "group_id") <- group_id
  attr(mdl, "stanmodel") <- stanmodels[["simplified_model"]]
  attr(mdl, "parameter_names") <- c("p", "shape", "median_time_to_response", "monthly_rate")
  class(mdl) <- c("IndependentMixtureCureRateModel", "Model", class(mdl))
  return(mdl)
}



#' @import patchwork
#'
#' @export
plot.IndependentMixtureCureRateModel <- function(x, sample = NULL, data = NULL, warmup = 250L, nsim = 1e4, seed = NULL, ...) {
  if (is.null(sample)) {
    if (is.null(data)) {
      sample <- sample_prior(x, warmup = warmup, nsim = nsim, seed = seed, ...)
    } else {
      sample <- sample_posterior(x, data = data, warmup = warmup, nsim = nsim, seed = seed, ...)
    }
  }

  p1 <- sample %>%
    filter(parameter == "p") %>%
    ggplot() +
    aes(value) +
    geom_histogram(
      bins = 25,
      alpha = 0.5,
      aes(fill = group_id), position = position_identity()
    ) +
    scale_x_continuous("event probability") +
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )

  p2 <- sample %>%
    filter(parameter %in% c("shape", "median_time_to_response")) %>%
    tidyr::pivot_wider(names_from = "parameter") %>%
    mutate(
      scale = median_time_to_response/(log(2)^(1/shape)),
      t = rweibull(n(), shape, scale)
    ) %>%
    ggplot() +
    aes(t) +
    geom_histogram(
      bins = 25,
      alpha = 0.5,
      aes(fill = group_id), position = position_identity()
    ) +
    scale_x_continuous("months to event") +
    coord_cartesian(xlim = c(0, NA_real_)) +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )

  p3 <- sample %>%
    filter(parameter == "monthly_rate") %>%
    ggplot() +
    aes(value) +
    geom_histogram(
      bins = 25,
      alpha = 0.5,
      aes(fill = group_id), position = position_identity()
    ) +
    scale_x_continuous("monthly recruitment rate") +
    coord_cartesian(xlim = c(0, NA_real_)) +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )

  (p1 + p2 + p3) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
}
