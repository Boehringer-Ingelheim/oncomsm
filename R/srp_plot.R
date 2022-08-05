#' @export
srp_plot <- function(tbl_mstate, starting_state = "stable") {

  tbl_mstate <- tbl_mstate %>%
    mutate(subject_id = factor(subject_id))

  tbl_points <- tbl_mstate %>%
    filter(t_max != -Inf) %>%
    mutate(
      tmp = purrr::pmap(
          list(from, to, t_min, t_max, t_sot),
          ~tibble(t = c(..3, ..4, ..5), state = c(..1, ..2, starting_state))
        )
    ) %>%
    select(subject_id, group_id, tmp) %>%
    tidyr::unnest(tmp) %>%
    filter(is.finite(t)) %>%
    distinct() %>%
    arrange(subject_id, t)

  tbl_intervals <- tbl_mstate %>%
    bind_rows(
      select(tbl_mstate, subject_id, group_id, t_sot) %>%
        distinct() %>%
        mutate(from = starting_state, to = starting_state, t_min = t_sot, t_max = t_sot)
    ) %>%
    arrange(subject_id, t_min, t_max) %>%
    distinct() %>%
    group_by(subject_id) %>%
    transmute(
      subject_id,
      group_id,
      state = if_else(to == lead(from), lead(from), NA_character_),
      tmp1 = t_max,
      tmp2 = lead(t_min)
    ) %>%
    ungroup() %>%
    filter(!is.na(state), is.finite(tmp1), is.finite(tmp2), tmp2 > tmp1)

  tbl_at_risk <- tbl_mstate %>%
    filter(t_max == Inf) %>%
    transmute(
      subject_id,
      group_id,
      t = t_min,
      state = from
    )

  tbl_censored <- tbl_mstate %>%
    filter(t_max == -Inf) %>%
    transmute(
      subject_id,
      group_id,
      t = t_min,
      state = from
    )

  ggplot() +
    geom_segment(aes(x = tmp1, xend = tmp2, y = subject_id, yend = subject_id, color = state), data = tbl_intervals) +
    geom_point(aes(t, subject_id, color = state), data = tbl_points) +
    geom_segment(aes(t, subject_id, xend = t + .25, yend = subject_id, color = state), arrow = arrow(type = "closed", angle = 10), data = tbl_at_risk) +
    geom_point(aes(t, subject_id, color = state), shape = "x", size = 5, data = tbl_censored) +
    labs(x = "Time since SoT", y = "Subject ID") +
    scale_color_discrete("") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "top")

}
