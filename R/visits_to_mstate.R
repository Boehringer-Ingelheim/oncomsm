#' Convert longitudinal visit data to time-to-event data
#'
#' @param tbl_visits visit data in long format; columns "group_id" (character),
#' "subject_id" (character), "t" (numeric, visit time point in months since start of trial),
#' "status" (character, S/P/R), and "eof"
#' (logical, end of follow-up, is the visit the last visit of the individual?)
#' @param cutoff maximal visit times to consider - allows retrospective analyses.
#' @param event character, string encoding an event of interest
#' @param nonevent character, string encoding a competing event
#'
#' @return A data frame with one row per subject in tbl_visits. columns:
#' "group_id", "subject_id", "t_recruitment" (recruitment time since start of trial),
#' "dt1" (minimal time since recruitment to event), "dt2" (maximal time since recruitment to event)
#' "dt_eof" (time since recruitment to end of follow up)
#'
#' @importFrom dplyr group_by summarize left_join
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#'
#' @export
visits_to_mstate <- function(tbl_visits, start_state, absorbing_states,
                             now = max(tbl_visits$t), eof_indicator = "EOF") {

  tbl_mstate <- list()

  subject_id_lagged <- 0L
  state_lagged <- 0L
  t_sot <- 0 # start of treatment
  for (i in 1:nrow(tbl_visits)) {
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
          dt_min = tbl_visits$t[i] - t_sot,
          dt_max = Inf,
          t_sot = t_sot
        ))
        state_lagged <- tbl_visits$state[i]
      } else { # record jump
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
            subject_id = tbl_visits$subject_id[i],
            group_id = tbl_visits$group_id[i],
            from = state_lagged,
            to = tbl_visits$state[i],
            dt_min = tbl_visits$t[i - 1] - t_sot,
            dt_max = tbl_visits$t[i] - t_sot,
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
        dt_min = now - t_sot,
        dt_max = Inf,
        t_sot = t_sot
      ))
    }
  }
  return(tbl_mstate)
}
