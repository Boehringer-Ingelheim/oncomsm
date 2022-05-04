#' Generate visit data
#'
#' Generate example data from mixture cure rate model with Weibull distributed
#' time-to-response and time-to-progression.
#'
#' @param group_id ids of arms
#' @param n sample size of arms
#' @param response_rate response rates of arms
#' @param delay delay (in months) of recruitment per arm
#' @param recruitment_rate recruitment rates (per month) of arms
#' @param visit_spacing distance between visits (in months)
#' @param max_duration maximal overall duration of all arms
#' @param responder_weibull_scale parameter of Weibull distribution for TTE responders
#' @param responder_weibull_shape parameter of Weibull distribution for TTE responders
#' @param nonresponder_weibull_scale parameter of Weibull distribution for TTE non-responders (progression)
#' @param nonresponder_weibull_shape parameter of Weibull distribution for TTE non-responders (progression)
#' @param seed random seed to use, NULL uses no seed
#'
#' @return A data frame with sampled data in columns "group_id", "subject_id",
#' "t" (visit time point in months since start of trial), "status" (S(table), R(esponse), P(rogression)),
#' "eof" (end of follow-up, is the visit the last visit of the individual?)
#'
#' @export
generate_visit_data <- function(
  group_id,
  n,
  response_rate,
  delay = numeric(length(group_id)),
  recruitment_rate,
  visit_spacing,
  max_duration,
  responder_weibull_scale,
  responder_weibull_shape,
  nonresponder_weibull_scale,
  nonresponder_weibull_shape,
  seed = NULL
) {

  if (!is.null(seed))
    set.seed(seed)

  # function for individual arm
  f <- function(group_id, n, response_rate, delay, recruitment_rate,
                visit_spacing, max_duration, responder_weibull_scale, responder_weibull_shape,
                nonresponder_weibull_scale, nonresponder_weibull_shape
  ) {
    tibble::tibble(
      subject_id = uuid::UUIDgenerate(n = n), # generate unique patient identifiers
      responder = stats::rbinom(n, 1, response_rate)
    ) %>%
    mutate(
      t = delay + cumsum(stats::rexp(n(), recruitment_rate)), # sample recruitment times
      dt = purrr::map_dbl(.data$responder, function(responder) {
        stats::rweibull(
          1,
          scale = if (responder) responder_weibull_scale else nonresponder_weibull_scale,
          shape = if (responder) responder_weibull_shape else nonresponder_weibull_shape
        ) %>%
          pmax(visit_spacing + 1e-4) # no events before first visit (first visit needs to be SD)
      }
      )
    ) %>%
    tidyr::expand_grid(
      dt_visit = seq(0, max_duration, by = visit_spacing),
    ) %>%
    mutate(
      status = dplyr::case_when(
        .data$dt_visit < .data$dt ~ "S", # no event yet
        .data$dt_visit >= .data$dt & .data$responder ~ "R",
        .data$dt_visit >= .data$dt & !.data$responder ~ "P"
      )
    ) %>%
    group_by(
      .data$subject_id
    ) %>%
    filter( # restrict to first non SD visit
      .data$status == "S" | lag(.data$status) == "S",
      .data$t + .data$dt_visit <= max_duration
    ) %>%
    transmute( # last visit is end of treatment (more precisely end of follow up)
      group_id = group_id,
      .data$subject_id,
      t = .data$t + .data$dt_visit,
      .data$status,
      eof = row_number() == n()
    ) %>%
    ungroup()
  }
  res <- list()
  for (i in seq_along(group_id)) {
    res <- dplyr::bind_rows(
      res,
      f(group_id[i], n[i], response_rate[i], delay[i], recruitment_rate[i],
        visit_spacing[i], max_duration,
        responder_weibull_scale[i], responder_weibull_shape[i],
        nonresponder_weibull_scale[i], nonresponder_weibull_shape[i]
      )
    )
  }
  return(res)
}
