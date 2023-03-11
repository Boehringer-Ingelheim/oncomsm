#' Plot progression-free-survival function
#'
#' `plot_pfs()` plots the progression-free-survival function of a model.
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
#' @seealso [plot_transition_times()] [plot_response_probability()]
#'
#' @examples
#' \dontrun{
#' mdl <- create_srpmodel(A = define_srp_prior())
#' plot_pfs(mdl)
#' }
#' @export
plot_pfs <- function(model, # nolint
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
  if (!is.null(confidence)) {
    alpha_half <- (1 - confidence) / 2
  } else {
    alpha_half <- 0
  }
  tbl <- compute_pfs(
    model,
    t = dt_grid,
    parameter_sample = parameter_sample
  ) %>%
    # integrate over sample
    group_by(.data$group_id, .data$t) %>%
    summarize(
      lo = stats::quantile(.data$pfs, probs = alpha_half) %>% as.numeric(),
      hi = stats::quantile(.data$pfs, probs = 1 - alpha_half) %>% as.numeric(),
      pfs = mean(.data$pfs),
      .groups = "drop"
    )
  plt <- ggplot2::ggplot(tbl)
  if (!is.null(confidence)) {
    plt <- plt + ggplot2::geom_ribbon(ggplot2::aes(.data$t, ymin = .data$lo,
                                                   ymax = .data$hi,
                                                   fill = .data$group_id),
                                      alpha = 0.2)
  }
  plt <- plt + ggplot2::geom_line(ggplot2::aes(
    x = .data$t,
    y = .data$pfs, color = .data$group_id
  )) +
    ggplot2::labs(x = "time", y = "PFS") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .1)) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  attr(plt, "data") <- tbl
  return(plt)
}
