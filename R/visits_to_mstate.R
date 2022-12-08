#' Convert cross-sectional visit data to multi-state format
#'
#' `visits_to_mstate()` converts visits to interval-censored multi-state
#' data where each row corresponds to a transition between states.
#' The conversion assumes that visit spacing is tight enough to not miss any
#' transitions.
#'
#' @param tbl_visits data frame, visit data in long format
#' @template param-model
#' @param now time point since start of trial (might be later than last
#'   recorded visit)
#' @param eof_indicator state name indicating (exactly observed) end of
#'   follow up.
#'
#' @return A data frame with multi-state data; variables are
#' `subject_id<chr>`, `group_id<chr>`, `subject_id<chr>`, `from<chr>`,
#' `to<chr>`, `t_min<dbl>`, `t_max<dbl>`, `t_sot<dbl>`, where
#' `to` and `from` indicate the state from which and into which the transitions
#' occurs (stable, response, progression), `t_max` and `t_min` specify the
#' interval in which the transition occurred relative to `t_sot`
#' (start of treatment).
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' tbl_visits <- sample_predictive(mdl, n_per_group = 5L, nsim = 1, seed = 468L)
#' visits_to_mstate(tbl_visits, mdl)
#'
#' @export
visits_to_mstate <- function(tbl_visits, model, now = max(tbl_visits$t),
                             eof_indicator = "EOF") {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (!inherits(tbl_visits, "data.frame")) {
    stop("'tbl_visits' must be a data.frame") # nocov
  } else {
    checkmate::test_true(inherits(tbl_visits$subject_id, "character"))
    checkmate::test_true(inherits(tbl_visits$group_id, "character"))
    checkmate::test_true(inherits(tbl_visits$t, "numeric"))
    checkmate::test_true(inherits(tbl_visits$state, "character"))
  }
  checkmate::check_class(model, c("srpmodel", "list"))
  # make sure everything is sorted
  tbl_visits <- arrange(tbl_visits, .data$subject_id, .data$t)
  tbl_mstate <- list()
  subject_id_lagged <- 0L
  state_lagged <- 0L
  t_sot <- 0 # start of treatment
  for (i in seq_len(nrow(tbl_visits))) {
    if (tbl_visits$subject_id[i] != subject_id_lagged || i == 1) {
      # switch to new subject
      subject_id_lagged <- tbl_visits$subject_id[i]
      state_lagged <- "stable"
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
          # - Inf indicates censoring and end of follow up
          # (event can no longer be observed)
          t_max = -Inf,
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
    if (!(tbl_visits$state[i] %in% c("progression", eof_indicator))) {
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
        t_max = Inf,
        t_sot = t_sot
      ))
    }
  }
  if (!is.null(attr(tbl_visits, "isemptydata"))) {
    # this is only used in conjunction with .emptydata since the global now
    # does not make sense in that context; minimum time to first event is
    # assumed to be equal to visit spacing / 2
    tbl_mstate <- tbl_mstate %>%
      mutate(t_min = model$visit_spacing[.data$group_id] / 2 + .data$t_sot)
  }
  return(tbl_mstate)
}
