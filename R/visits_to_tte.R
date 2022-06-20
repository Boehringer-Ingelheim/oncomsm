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
visits_to_tte <- function(tbl_visits, cutoff = Inf, event = "E", nonevent = "N") {

  # save all individuals for later
  res <- dplyr::select(tbl_visits, .data$group_id, .data$subject_id) %>%
    dplyr::distinct()

  # restrict visit data to <= cutoff
  tbl_visits <- dplyr::filter(tbl_visits, .data$t <= cutoff)

  if (nrow(tbl_visits) > 0) {
    # extract recruitment times (first visit)
    tbl_recruitment_times <- tbl_visits %>%
      dplyr::group_by(.data$group_id, .data$subject_id) %>%
      dplyr::summarize(t_recruitment = min(.data$t))

    # computed lower and upper bound on time since recruitment to event
    f <- function(status, dt, eof) {
      idx_event <- which(status %in% c("N", "E"))
      idx_eof <- which(eof)
      if (length(idx_eof) == 0) {
        dt_eof <- NA_real_
      } else {
        dt_eof <- dt[idx_eof]
      }
      if (length(idx_event) == 0) {
        # no event, no non-event
        return(tibble(dt1 = utils::tail(dt, 1), dt2 = Inf, dt_eof = dt_eof))
      }
      if (status[idx_event[1]] == "E") {
        # event before nonevent
        return(tibble(dt1 = dt[idx_event[1] - 1], dt2 = dt[idx_event[1]], dt_eof = dt_eof))
      }
      if (status[idx_event[1]] == "N") {
        return(tibble(dt1 = Inf, dt2 = Inf, dt_eof = dt_eof))
      }
    }
    tbl_interval_censored_event_times_since_sot <- tbl_visits %>%
      dplyr::group_by(.data$group_id, .data$subject_id) %>%
      dplyr::mutate(dt = .data$t - min(.data$t)) %>%
      dplyr::summarize(
        tmp = f(.data$status, .data$dt, .data$eof),
        .groups = "drop"
      ) %>%
      tidyr::unnest(.data$tmp)

    # combine
    res %>%
      dplyr::left_join(tbl_recruitment_times, by = c("group_id", "subject_id")) %>%
      dplyr::left_join(tbl_interval_censored_event_times_since_sot, by = c("group_id", "subject_id")) %>%
      dplyr::select("group_id", "subject_id", "t_recruitment", "dt1", "dt2", "dt_eof") %>%
      return()
  } else {
    res %>%
      mutate(t_recruitment = NA_real_, dt1 = NA_real_, dt2 = NA_real_, dt_eof = NA_real_) %>%
      return()
  }

}
