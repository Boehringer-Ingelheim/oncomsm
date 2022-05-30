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
  attr(mdl, "parameter_names") <- c("p", "shape", "scale", "median_time_to_response", "monthly_rate")
  class(mdl) <- c("IndependentMixtureCureRateModel", "Model", class(mdl))
  return(mdl)
}

.impute.IndependentMixtureCureRateModel <- function(model, data, nsim, now = NULL, parameter_sample = NULL, warmup_parameters = 250L, nsim_parameters = 1000L, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # sample from prior/posterior if not given
  if (is.null(parameter_sample)) {
    parameter_sample <- .sample(
      model, data = data, warmup = warmup_parameters, nsim = nsim_parameters,
      seed = seed, rstan_output = TRUE, pars = c("p", "shape", "scale", "monthly_rate"), ...
    )
  }
  stopifnot(isa(parameter_sample, "stanfit"))
  # extract parameter matrices
  p <- rstan::extract(parameter_sample, "p")[[1]]
  shape <- rstan::extract(parameter_sample, "shape")[[1]]
  scale <- rstan::extract(parameter_sample, "scale")[[1]]
  monthly_rate <- rstan::extract(parameter_sample, "monthly_rate")[[1]]
  # figure out latest time entry from data
  if (all(is.na(data$t_recruitment))) {
    t_max_data <- 0
  } else {
    t_max_data <- data %>%
      mutate(t_max = t_recruitment + if_else(is.finite(dt2), dt2, ifelse(is.finite(dt1), dt1, 0))) %>%
      filter(is.finite(t_max)) %>%
      pull(t_max) %>%
      max()
  }
  # if "now" is not specified, use heuristic from data
  if (!is.null(now) ) {
    stopifnot(now >= t_max_data)
  } else {
    now <- t_max_data
  }
  # figure out last recruitment time per group
  t_last <- data %>%
    group_by(group_id) %>%
    summarize(
      t_last = max(c(0, t_recruitment), na.rm = TRUE)
    ) %>%
    pull(t_last)
  # convert groups to integers
  data$group_id <- as.integer(factor(data$group_id, levels = attr(model, "group_id")))

  ff <- function(data, p, shape, scale, monthly_rate, t_last, now) {

    for (i in 1:nrow(data)) {
      if (is.na(data$t_recruitment[i])) {
        group <- data$group_id[i]
        data$t_recruitment[i] <- t_last[group] + rtruncexp(monthly_rate[group], max(0, now - t_last[group]), Inf)
        t_last[group] <- data$t_recruitment[i]
      }
      if (!is.infinite(data$dt1[i]) & !is.finite(data$dt2[i])) {
        group <- data$group_id[i]
        if (rbinom(1, n = 1, prob = p[group]) == 1) {
          # event
          dtmin <- max(data$dt1[i], now - data$t_recruitment[i], na.rm = TRUE)
          dt <- rtruncweibull(shape[group], scale[group], dtmin, 999)
          dt2 <- 0
          while (dt2 < dt) {
            dt2 <- dt2 + model$visit_spacing[group]
          }
          data$dt2[i] <- dt2
          if (is.na(data$dt1[i]))
            data$dt1[i] <- dt2 - model$visit_spacing[group]
        } else {
          # no event
          data$dt1[i] <- Inf
          data$dt2[i] <- Inf
        }
      }
    }
    return(data)
  }

  res <- tibble(
      iter = 1:nsim,
      idx = sample.int(nrow(p), size = nsim),
      p = purrr::map(idx, ~p[., ]),
      shape = purrr::map(idx, ~shape[., ]),
      scale = purrr::map(idx, ~scale[., ]),
      monthly_rate = purrr::map(idx, ~monthly_rate[., ]),
    ) %>%
    mutate(
      data = purrr::pmap(
        list(p, shape, scale, monthly_rate),
        ~ff(data, ..1, ..2, ..3, ..4, t_last, now)
      )
    ) %>%
    select(iter, data) %>%
    tidyr::unnest(data)

  # convert groups back
  res$group_id <- as.character(factor(res$group_id, levels = 1:length(attr(model, "group_id")), labels = attr(model, "group_id")))

  return(res)
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
