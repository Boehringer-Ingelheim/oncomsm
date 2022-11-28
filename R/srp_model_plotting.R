#' @inheritParams plot_mstate
#' @name srp_model
#' @export
plot_mstate.srp_model <- function(data, model, now = max(tbl_mstate$t_max), # nolint
                                  relative_to_sot = TRUE, ...) {
  starting_state <- attr(model, "states")[1]

  tbl_mstate <- data %>%
    rename(`Group ID` = "group_id")

  if (relative_to_sot) {
    tbl_mstate <- tbl_mstate %>%
      mutate(
        t_min = .data$t_min - .data$t_sot,
        t_max = .data$t_max - .data$t_sot,
        t_sot = 0
      )
  }

  tbl_points <- tbl_mstate %>%
    mutate(
      tmp = purrr::pmap(
        list(.data$from, .data$to, .data$t_min, .data$t_max, .data$t_sot),
        ~ if (is.na(..2)) {
          tibble(t = c(..3, ..5), state = c(..1, starting_state))
        } else {
          tibble(t = c(..3, ..4, ..5), state = c(..1, ..2, starting_state))
        }
      )
    ) %>%
    select(all_of(c("subject_id", "Group ID", "tmp"))) %>%
    tidyr::unnest("tmp") %>%
    filter(is.finite(.data$t), .data$t < now) %>%
    distinct() %>%
    arrange(.data$subject_id, .data$t)

  tbl_intervals <- tbl_mstate %>%
    bind_rows(
      select(tbl_mstate, all_of(c("subject_id", "Group ID", "t_sot"))) %>%
        distinct() %>%
        mutate(
          from = starting_state,
          to = starting_state,
          t_min = .data$t_sot,
          t_max = .data$t_sot
        )
    ) %>%
    arrange(.data$subject_id, .data$t_min, .data$t_max) %>%
    distinct() %>%
    group_by(.data$subject_id) %>%
    transmute(
      .data$subject_id,
      .data$`Group ID`,
      state = if_else(.data$to == lead(.data$from), lead(.data$from),
        NA_character_
      ),
      tmp1 = .data$t_max,
      tmp2 = lead(.data$t_min)
    ) %>%
    ungroup() %>%
    filter(
      !is.na(.data$state), is.finite(.data$tmp1), is.finite(.data$tmp2),
      .data$tmp2 > .data$tmp1
    )

  tbl_at_risk <- tbl_mstate %>%
    filter(.data$t_max == Inf) %>%
    transmute(
      .data$subject_id,
      .data$`Group ID`,
      t = .data$t_min,
      state = .data$from
    )

  tbl_censored <- tbl_mstate %>%
    filter(.data$t_max == -Inf) %>%
    transmute(
      .data$subject_id,
      .data$`Group ID`,
      t = .data$t_min,
      state = .data$from
    )

  scale <- max(tbl_points$t)

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data$tmp1, xend = .data$tmp2,
        y = .data$subject_id, yend = .data$subject_id,
        color = .data$state
      ),
      data = tbl_intervals
    ) +
    ggplot2::geom_point(ggplot2::aes(.data$t, .data$subject_id,
      color = .data$state
    ), data = tbl_points) +
    ggplot2::geom_segment(
      ggplot2::aes(.data$t, .data$subject_id,
        xend = .data$t + scale / 33,
        yend = .data$subject_id, color = .data$state
      ),
      arrow = ggplot2::arrow(
        type = "closed", angle = 10,
        length = ggplot2::unit(0.05, "npc")
      ),
      data = tbl_at_risk
    ) +
    ggplot2::geom_point(ggplot2::aes(.data$t, .data$subject_id,
      color = .data$state
    ),
    shape = "x",
    size = 5, data = tbl_censored
    ) +
    ggplot2::geom_vline(xintercept = now) +
    ggplot2::labs(x = if (relative_to_sot) {
      "Time since SoT"
    } else {
      "Time since first SoT"
    }, y = "Subject ID") +
    ggplot2::scale_color_discrete("") +
    ggplot2::facet_wrap(~ .data$`Group ID`,
      ncol = 1,
      labeller = ggplot2::label_both,
      strip.position = "right", scales = "free_y"
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 5),
      legend.position = "right"
    )
}



#' @inheritParams plot_transition_times
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#' @template param-dt-params
#'
#' @rdname srp_model
#' @export
plot_transition_times.srp_model <- function(model, # nolint
                                            parameter_sample = NULL,
                                            seed = 42L,
                                            nsim = 500L,
                                            warmup = 250,
                                            nuts_control = list(),
                                            dt_interval = NULL,
                                            dt_n_grid = 25,
                                            dt_expand = 1.1,
                                            dt_grid = NULL,
                                            ...) {
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
      survival = mean(.data$survival),
      lo = quantile(.data$survival, probs = 0.1) %>% as.numeric(),
      hi = quantile(.data$survival, probs = 0.9) %>% as.numeric(),
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
  plt <- ggplot2::ggplot(tbl) +
    ggplot2::geom_ribbon(ggplot2::aes(.data$dt, ymin = .data$lo,
                                      ymax = .data$hi)) +
    ggplot2::geom_line(ggplot2::aes(.data$dt, .data$survival,
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



#' @inheritParams plot_response_probability
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#'
#' @rdname srp_model
#' @export
plot_response_probability.srp_model <- function(model, # nolint
                                               parameter_sample = NULL,
                                               seed = 42L,
                                               nsim = 500L,
                                               warmup = 250,
                                               nuts_control = list(),
                                               ...) {
  if (is.null(parameter_sample)) { # sample parameters from prior if none given
    parameter_sample <- sample_prior(model,
                                     seed = seed, nsim = nsim,
                                     warmup = warmup,
                                     nuts_control = nuts_control, ...)
  }
  tbl <- parameter_sample_to_tibble(model, parameter_sample) %>%
    filter(.data$parameter == "p")
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



#' @inheritParams plot_pfs
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#' @template param-dt-params
#'
#' @rdname srp_model
#' @export
plot_pfs.srp_model <- function(model, # nolint
                               parameter_sample = NULL,
                               seed = 42L,
                               nsim = 500L,
                               warmup = 250,
                               nuts_control = list(),
                               dt_interval = NULL,
                               dt_n_grid = 25,
                               dt_expand = 1.1,
                               dt_grid = NULL,
                               ...) {
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
  # TODO make sure this works from 0
  tbl <- compute_pfs(
      model,
      t = dt_grid,
      parameter_sample = parameter_sample
    ) %>%
    # integrate over sample
    group_by(.data$group_id, .data$t) %>%
    summarize(
      pfs = mean(.data$pfs),
      .groups = "drop"
    )
  plt <- ggplot2::ggplot(tbl) +
    ggplot2::geom_line(ggplot2::aes(
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



#' @param x the model to plot
#' @template param-parameter_sample
#' @template param-seed
#' @template param-nsim
#' @template param-warmup
#' @template param-nuts_control
#' @template param-dt-params
#' @template param-dotdotdot
#'
#' @rdname srp_model
#' @export
plot.srp_model <- function(x,
                           parameter_sample = NULL,
                           seed = 42L,
                           nsim = 500L,
                           warmup = 250,
                           nuts_control = list(),
                           dt_interval = NULL,
                           dt_n_grid = 25,
                           dt_expand = 1.1,
                           dt_grid = NULL,
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
  plt_trans <- plot_transition_times(x, parameter_sample, dt_grid = dt_grid)
  plt_pfs <- plot_pfs(x, parameter_sample, dt_grid = dt_grid)
  plt_pr <- plot_response_probability(x, parameter_sample)
  design <- "
  111
  234
  "
  plt_trans + plt_pfs + plt_pr + patchwork::guide_area() +
    patchwork::plot_layout(design = design, guides = "collect")
}
