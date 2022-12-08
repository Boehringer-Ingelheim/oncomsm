#' Swimmer plot of multi-state data
#'
#' `plot_mstate()` plots data in 'multi-state-format' as swimmer plot.
#'
#' @template param-model
#' @param data a data frame with multi-state data; variables are
#' `subject_id<chr>`, `group_id<chr>`, `subject_id<chr>`, `from<chr>`,
#' `to<chr>`, `t_min<dbl>`, `t_max<dbl>`, `t_sot<dbl>`, where
#' `to` and `from` indicate the state from which and into which the transitions
#' occurs (stable, response, progression), `t_max` and `t_min` specify the
#' interval in which the transition occurred relative to `t_sot`
#' (start of treatment).
#' @param now the current time relative to the start of the trial
#' @param relative_to_sot logical, should the timeline be relative to the start
#'   of trial or the start of treatment for each individual
#' @template param-dotdotdot
#'
#' @return a [ggplot2::ggplot] object
#'
#' @seealso [visits_to_mstate()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' tbl_visits <- sample_predictive(mdl, n_per_group = 5L, nsim = 1, seed = 468L)
#' tbl_mstate <- visits_to_mstate(tbl_visits, mdl)
#' plot_mstate(tbl_mstate, mdl)
#'
#' @export
plot_mstate <- function(data, model, now = max(tbl_mstate$t_max), # nolint
                              relative_to_sot = TRUE, ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  starting_state <- model$states[1]
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
