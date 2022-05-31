#' @param data a data frame with variables "subject_id", "group_id", "t_recruitment", "dt1" and "dt2"
#' where dt1 is the minimal and dt2 the maximal time to the event in question.
#' If both are Inf, the individual has not and will never experience the event.
#' If dt1 < Inf but dt2 == Inf the individual is still at risk.
#' "t_recruitment", "dt1" and "dt2" can also be missing if individual is not yet recruited.
#'
