#' Plot the response probability distributions
#'
#' `plot_response_probability()` plots the distribution over the response
#' probability parameter in the specified model.
#'
#' @template param-model
#' @template param-parameter_sample
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#' @template param-nuts_control
#' @template param-dotdotdot
#'
#' @return a [ggplot2::ggplot] object
#'
#' @seealso [plot_transition_times()] [plot_pfs()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' plot_response_probability(mdl)
#'
#' @export
plot_response_probability <- function(model, # nolint
                                      parameter_sample = NULL,
                                      seed = 42L,
                                      nsim = 500L,
                                      warmup = 250,
                                      nuts_control = list(),
                                      ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (is.null(parameter_sample)) { # sample parameters from prior if none given
    parameter_sample <- sample_prior(model,
                                     seed = seed, nsim = nsim,
                                     warmup = warmup,
                                     nuts_control = nuts_control, ...)
  }
  tbl <- parameter_sample_to_tibble(model, parameter_sample) %>%
    filter(.data$parameter == "p") %>%
    select(-"transition")
  plt <- ggplot2::ggplot(tbl) +
    ggplot2::stat_ecdf(ggplot2::aes(.data$value, color = .data$group_id),
                       geom = "line") +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(x = "response probability", y = "CDF") +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = .1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = .1)
    ) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  attr(plt, "data") <- tbl
  return(plt)
}
