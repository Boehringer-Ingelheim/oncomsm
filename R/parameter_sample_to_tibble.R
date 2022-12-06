#' Convert parameter sample to data table
#'
#' `parameter_sample_to_tibble()` takes a [rstan::stanfit] parameter sample of
#' a model, extracts the paramters values and returns them in a data frame.
#'
#' @template param-model
#' @template param-sample
#' @template param-dotdotdot
#'
#' @return a tibble with the sampled parameters, in long format
#'
#' @seealso [sample_prior()] [sample_posterior()]
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' smpl <- sample_prior(mdl, seed = 3647L)
#' parameter_sample_to_tibble(mdl, smpl)
#'
#' @export
parameter_sample_to_tibble <- function(model, sample, ...) { # nolint
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  stopifnot(isa(sample, "stanfit"))
  as.matrix(sample) %>%
    as_tibble() %>%
    mutate(
      iter = row_number()
    ) %>%
    tidyr::pivot_longer(-"iter") %>%
    filter(.data$name != "lp__") %>%
    tidyr::separate("name",
                    into = c("parameter", "group_id"),
                    sep = "\\[", fill = "right"
    ) %>%
    tidyr::separate("group_id",
                    into = c("group_id", "transition"),
                    sep = "[\\]|,]", fill = "right", extra = "drop"
    ) %>%
    mutate(
      group_id = model$group_id[as.integer(stringr::str_extract(
        .data$group_id, "[0-9]+"
      ))],
      transition = as.integer(.data$transition)
    )
}
