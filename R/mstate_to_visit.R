#' @export
mstate_to_visits <- function(model, tbl_mstate, ...) {
  UseMethod("mstate_to_visits")
}

#' @export
mstate_to_visits.Model <- function(model, tbl_mstate, ...) {
  stop("not implemented")
}


#' @export
mstate_to_visits.srp_model <- function(model, tbl_mstate, ...) {

  message("'mstate_to_visits.srp_model' is experimental")

  starting_state <- attr(model, "states")[1]
  visit_spacing <- attr(model, "visit_spacing")[1] # TODO: respect groups

  tbl_mstate %>% select(
      .data$subject_id,
      .data$group_id,
      .data$from,
      .data$to,
      .data$t_min,
      .data$t_max,
      .data$t_sot
    ) %>%
    bind_rows(
      tbl_mstate %>%
        select(.data$subject_id, .data$group_id, .data$t_sot) %>%
        distinct() %>%
        mutate(from = starting_state, to = starting_state, t_min = .data$t_sot, t_max = .data$t_sot)
    ) %>%
    arrange(.data$subject_id, .data$t_min, .data$t_max) %>%
    select(-.data$t_sot) %>%
    distinct() %>%
    group_by(.data$subject_id) %>%
    mutate(
      a = .data$t_max,
      b = lead(.data$t_min) %>% if_else(is.na(.), .data$a, .),
      c = .data$to,
      d = lead(.data$from)
    ) %>%
    filter(is.finite(.data$a), is.finite(.data$b), !is.na(.data$to)) %>%
    rowwise() %>%
    mutate(
      tmp = purrr::pmap(
        list(.data$c, .data$d, .data$a, .data$b),
        ~{
          if (a == b) {
            res <- tibble(t = ..3, state = ..1)
          } else {
            res <- tibble(t = seq(..3, max(..3, ..4), length.out = max(2, (..4 - ..3)/visit_spacing)), state = rep(..1, length(t)))
            res$state[nrow(res)] <- ..2
          }
          res
        }
      )
    ) %>%
    tidyr::unnest(.data$tmp) %>%
    select(.data$subject_id, .data$group_id, .data$t, .data$state) %>%
    ungroup()

}
