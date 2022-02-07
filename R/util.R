.tbl_to_stan_data <- function(model, data, n_additional) {
  UseMethod(".tbl_to_stan_data")
}

.prior_to_stan_data <- function(model) {
  UseMethod(".prior_to_stan_data")
}

check_input_data <- function(tbl_data) {
  if (is.null(tbl_data))
    tbl_data <- tibble::tibble(
      subject_id = integer(0L),
      t_recruitment = numeric(0L),
      dt_response_min = numeric(0L),
      dt_response_max = numeric(0L)
    )
  tmp <- dplyr::select(
    tbl_data,
    subject_id, t_recruitment, dt_response_min, dt_response_max
  )
  assertthat::assert_that(
    all(purrr::map_chr(tbl_data, class) == c("integer", "numeric", "numeric", "numeric"))
  )
  assertthat::assert_that(
    length(tbl_data$subject_id) == length(unique(tbl_data$subject_id))
  )
  assertthat::assert_that(
    all(tbl_data$t_recruitment > 0 | is.na(tbl_data$t_recruitment)),
    all(tbl_data$dt_response_min <= tbl_data$dt_response_max | !is.finite(tbl_data$dt_response_max))
  )
  tbl_data
}
