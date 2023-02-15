#' Check a visits data set for correct format
#'
#' Raises specific errors when encountering issues in the data.
#'
#' @param data data.frame to check
#' @param model [srpmodel] object used to fit data
#'
#' @return data.frame, same as input but all censoring events after terminal
#' states are removed.
#'
#' @examples
#' tbl <- data.frame(group_id = "A", subject_id = "A1", t = 0, state = "stable")
#' mdl <- create_srpmodel(A = define_srp_prior())
#' check_data(tbl, mdl)
#' @export
check_data <- function(data, model) {
  checkmate::assert_class(data, "data.frame")
  checkmate::assert_class(model, "srpmodel")
  # check column names, do we have everything?
  expected_columns <- c("group_id", "subject_id", "t", "state")
  names_missing <- expected_columns[!(expected_columns %in% names(data))]
  if (length(names_missing) > 0) {
    stop(sprintf("columns: %s; missing from data", paste(names_missing, collapse = ", "))) # nolint
  }
  checkmate::test_true(inherits(data$subject_id, "character"))
  checkmate::test_true(inherits(data$group_id, "character"))
  checkmate::test_true(inherits(data$t, "numeric"))
  checkmate::test_true(inherits(data$state, "character"))
  # check no missing data
  if (any(!stats::complete.cases(data))) {
    stop("no missing data in visits allowed")
  }
  # check sorted, by patient and t
  sid <- factor(data$subject_id, levels = unique(data$subject_id))
  properly_sorted <- all(diff(order(sid, data$t)) > 0)
  if (!properly_sorted) {
    stop("data needs to be sorted by 'subject_id' and 't'")
  }
  # assert labels
  if (!all(data$state %in% c(model$states, model$censored))) {
    stop("state must be consistent with model specification.")
  }
  # assert order of states
  states_sorted <- data %>%
    group_by(.data$subject_id) %>%
    summarize(
      sorted = factor(.data$state, levels = model$states) %>% # ignore censored here!
        as.integer() %>%
        diff() %>%
        {. >= 0} %>%
        all(na.rm = TRUE)
    ) %>%
    pull(.data$sorted) %>%
    all()
  if (!states_sorted) {
    stop("no reverse jumps in state are allowed")
  }
  # assert positive times
  if (any(data$t < 0)) {
    stop("'t' must be non-negative")
  }
  # handle eof after terminal state
  data_ <- data %>%
    group_by(.data$subject_id) %>%
    filter({
      idx <- which(.data$state %in% c(model$states[3], model$censored))
      if (length(idx) > 0) {
        row_number() <= idx[1]
      } else {
        TRUE
      }
    }) %>%
    ungroup()
  if (nrow(data_) < nrow(data)) {
    message("Removed censoring events after terminal state, make sure this is intended") # nolint
  }
  return(data_)
}
