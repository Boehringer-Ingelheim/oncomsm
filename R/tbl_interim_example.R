#' A Interim-Analysis Example Data Set
#'
#' This data set consists of data from a three-arm study at time point
#' 10 months after start of recruitment.
#' Jump times of this multistate data set are properly censored.
#'
#' The data was generated with a recruitment rate of 3/per month (overall) and
#' aiming at 30 individuals per group.
#'
#' @format ## `tbl_interim_example`
#' A data frame with 45 rows and 7 columns:
#' \describe{
#'   \item{subject_id}{subject identifier}
#'   \item{group_id}{group identifier}
#'   \item{from}{starting state of transition}
#'   \item{to}{end state of transition, missing if censored}
#'   \item{t_min}{lower boundary on jump time (in months)}
#'   \item{t_max}{upper boundary on jump time (in months)}
#'   \item{t_sot}{recrutiment time since start of study (in months)}
#' }
"tbl_interim_example"
