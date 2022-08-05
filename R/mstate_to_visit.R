#' @export
mstate_to_visits <- function(tbl_mstate, starting_state = "stable", visit_spacing = 1) {

  select(tbl_mstate, subject_id, group_id, from, to, t_min, t_max, t_sot) %>%
    bind_rows(
      select(tbl_mstate, subject_id, group_id, t_sot) %>%
        distinct() %>%
        mutate(from = starting_state, to = starting_state, t_min = t_sot, t_max = t_sot)
    ) %>%
    arrange(subject_id, t_min, t_max) %>%
    select(-t_sot) %>%
    distinct() %>%
    group_by(subject_id) %>%
    mutate(
      a = t_max,
      b = lead(t_min) %>% if_else(is.na(.), a, .),
      c = to, d = lead(from)
    ) %>%
    filter(is.finite(a), is.finite(b), !is.na(to)) %>%
    rowwise() %>%
    mutate(
      tmp = purrr::pmap(
        list(c, d, a, b),
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
    unnest(tmp) %>%
    select(subject_id, group_id, t, state)

}
