#' Plot the transition times of a model
#'
#' `plot_transition_times()` plots a the survival functions for the transition
#' times in a multi-state model.
#'
#' @template param-model
#' @template param-parameter_sample
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#' @template param-dt-params
#' @template param-nuts_control
#' @param confidence numeric in (0, 1) confidence level for point-wise
#' confidence bands around mean; none plotted if NULL.
#' @template param-dotdotdot
#'
#' @return a [ggplot2::ggplot] object
#'
#' @seealso [plot_pfs()] [plot_response_probability()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' plot_transition_times(mdl)
#'
#' @export
plot_transition_times <- function(model, # nolint
                                  parameter_sample = NULL,
                                  seed = 42L,
                                  nsim = 500L,
                                  warmup = 250,
                                  nuts_control = list(),
                                  dt_interval = NULL,
                                  dt_n_grid = 25,
                                  dt_expand = 1.1,
                                  dt_grid = NULL,
                                  confidence = NULL,
                                  ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (is.null(parameter_sample)) { # sample parameters from prior if none given
    parameter_sample <- sample_prior(model,
                                     seed = seed, nsim = nsim,
                                     warmup = warmup,
                                     nuts_control = nuts_control, ...)
  }
  if (is.null(dt_grid)) {
    # determine plotting grid
    dt_grid <- get_dt_grid(model, parameter_sample, dt_interval,
                           dt_n_grid, dt_expand, seed)
  }
  # convert parameters to tibble and prepare for plotting
  if (!is.null(confidence)) {
    alpha_half <- (1 - confidence) / 2
  } else {
    alpha_half <- 0
  }
  tbl <- parameter_sample_to_tibble(model, parameter_sample) %>%
    filter(.data$parameter %in% c("shape", "scale")) %>%
    tidyr::pivot_wider(
      names_from = "parameter",
      values_from = "value"
    ) %>%
    tidyr::expand_grid(dt = dt_grid) %>%
    mutate(
      survival = 1 - stats::pweibull(.data$dt,
                                     shape = .data$shape, scale = .data$scale)
    ) %>%
    group_by(.data$group_id, .data$transition, .data$dt) %>%
    summarize(
      lo = stats::quantile(.data$survival, probs = alpha_half) %>%
        as.numeric(),
      hi = stats::quantile(.data$survival, probs = 1 - alpha_half) %>%
        as.numeric(),
      survival = mean(.data$survival),
      .groups = "drop"
    ) %>%
    filter(
      is.finite(.data$survival)
    ) %>%
    mutate(
      transition = case_when(
        .data$transition == 1 ~ "stable to response",
        .data$transition == 2 ~ "stable to progression",
        .data$transition == 3 ~ "response to progression"
      ) %>%
        factor(levels = c(
          "stable to response", "stable to progression",
          "response to progression"
        ))
    )
  plt <- ggplot2::ggplot(tbl)
  if (!is.null(confidence)) {
    plt <- plt + ggplot2::geom_ribbon(
      ggplot2::aes(.data$dt, ymin = .data$lo,
                   ymax = .data$hi, fill = .data$group_id),
      alpha = 0.2
    )
  }
  plt <- plt + ggplot2::geom_line(ggplot2::aes(.data$dt, .data$survival,
                                               color = .data$group_id)) +
    ggplot2::labs(x = "time to next event", y = "'Survival' fraction") +
    ggplot2::scale_color_discrete("") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = .1)
    ) +
    ggplot2::facet_wrap(~ .data$transition, nrow = 1) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(1.5, "lines")
    )
  # add plotting data to return value
  attr(plt, "data") <- tbl
  return(plt)
}
