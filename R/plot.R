#' Summary plot of model prior
#'
#' @param x the model to plot
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#' @template param-nuts_control
#' @template param-dt-params
#' @template param-dotdotdot
#' @param confidence numeric in (0, 1) confidence level for point-wise
#' confidence bands around mean; none plotted if NULL.
#'
#' @return A `patchwork` object, see [patchwork::patchwork]
#'
#' @examples
#' \dontrun{
#' mdl <- create_srpmodel(A = define_srp_prior())
#' plot(mdl)
#' }
#' @seealso [plot_pfs()] [plot_transition_times()]
#' [plot_response_probability()]
#'
#' @export
plot.srpmodel <- function(x,
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
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("the patchwork package is required to plot SRP models") # nocov
  }
  if (is.null(parameter_sample)) { # sample parameters from prior if none given
    parameter_sample <- sample_prior(x,
                                     seed = seed, nsim = nsim,
                                     warmup = warmup,
                                     nuts_control = nuts_control, ...)
  }
  if (is.null(dt_grid)) {
    # determine plotting grid
    dt_grid <- get_dt_grid(x, parameter_sample, dt_interval,
                           dt_n_grid, dt_expand, seed)
  }
  plt_trans <- plot_transition_times(x, parameter_sample, dt_grid = dt_grid,
                                     confidence = confidence)
  plt_pfs <- plot_pfs(x, parameter_sample, dt_grid = dt_grid,
                      confidence = confidence)
  plt_pr <- plot_response_probability(x, parameter_sample)
  design <- "
  111
  234
  "
  plt_trans + plt_pfs + plt_pr + patchwork::guide_area() +
    patchwork::plot_layout(design = design, guides = "collect")
}
