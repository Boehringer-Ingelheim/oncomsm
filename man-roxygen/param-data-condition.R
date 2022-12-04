#' @param data a data frame with variables
#' `subject_id<chr>` (subject identifier),
#' `group_id<chr>` (group identifier),
#' `t<dbl>` (time of visit, relative to first visit in study),
#' `state<chr>` (state recorded at visit).
#' Allowed states are "stable", "response", "progression" (or death),
#' and "EOF" (end of follow-up).
#' The EOF state marks the end of an individual's follow-up before the absorbing
#' state "progression".
