#' Convert longitudinal visit data to time-to-event data
#'
#' @param tbl_visits visit data in long format
#' @param start_state staring state
#' @param absorbing_states character vector of absorbing states
#' @param now time point since start of trial (might be later than last recorded visit) # nolint
#' @param eof_indicator state name indicating (exactly observed) eond of follow up. # nolint
#'
#' @return A data frame
#'
#' @importFrom dplyr group_by summarize left_join
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#'
#' @export
visits_to_mstate <- function(tbl_visits, start_state, absorbing_states,
                             now = max(tbl_visits$t), eof_indicator = "EOF") {

  tbl_visits <- arrange(tbl_visits, .data$subject_id, .data$t) # make sure everything is sorted # nolint

  tbl_mstate <- list()

  subject_id_lagged <- 0L
  state_lagged <- 0L
  t_sot <- 0 # start of treatment
  for (i in seq_len(nrow(tbl_visits))) {
    if (tbl_visits$subject_id[i] != subject_id_lagged || i == 1) {
      # switch to new subject
      subject_id_lagged <- tbl_visits$subject_id[i]
      state_lagged <- start_state
      t_sot <- tbl_visits$t[i]
      if (tbl_visits$state[i] != state_lagged) {
        # record jump
        stop(sprintf(
          "first visit must be in starting state; subject_id=%s, state=%s",
          tbl_visits$subject_id[i],
          tbl_visits$state[i]
        ))
      }
    }
    # handle jumps
    if (tbl_visits$state[i] != state_lagged) {
      if (tbl_visits$state[i] == eof_indicator) { # record eof
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
          subject_id = tbl_visits$subject_id[i],
          group_id = tbl_visits$group_id[i],
          from = state_lagged,
          to = NA,
          t_min = tbl_visits$t[i],
          t_max = -Inf, # - Inf indicates censoring and end of follow up (event can no longer be observed) # nolint
          t_sot = t_sot
        ))
        state_lagged <- tbl_visits$state[i]
      } else { # record jump
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
            subject_id = tbl_visits$subject_id[i],
            group_id = tbl_visits$group_id[i],
            from = state_lagged,
            to = tbl_visits$state[i],
            t_min = tbl_visits$t[i - 1],
            t_max = tbl_visits$t[i],
            t_sot = t_sot
          ))
      }
      # updated current state
      state_lagged <- tbl_visits$state[i]
    }
    # handle non-eof censoring
    censored <- FALSE
    if (!(tbl_visits$state[i] %in% c(absorbing_states, eof_indicator))) {
      if (i < nrow(tbl_visits)) {
        if (tbl_visits$subject_id[i] != tbl_visits$subject_id[i + 1]) {
          censored <- TRUE
        }
      } else {
        censored <- TRUE
      }
    }
    if (censored) {
      tbl_mstate <- bind_rows(tbl_mstate, tibble(
        subject_id = tbl_visits$subject_id[i],
        group_id = tbl_visits$group_id[i],
        from = tbl_visits$state[i],
        to = NA,
        t_min = now,
        t_max = Inf, # Inf indicates censoring while still at risk (event can still be observed) # nolint
        t_sot = t_sot
      ))
    }
  }
  return(tbl_mstate)
}
