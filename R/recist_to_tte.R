#' @export
recist_to_tte <- function(tbl_data, cutoff_date = Inf) {
  subject_ids <- levels(tbl_data$subject_id)
  dplyr::filter(tbl_data, visit_date <= cutoff_date) %>%
    dplyr::mutate(
      visit_t = difftime(visit_date, min(visit_date), units = "days") %>% {as.numeric(.)/30},
      subject_id = droplevels(subject_id)
    ) %>%
    recist_to_tte_ %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      subject_id = factor(subject_id, levels = subject_ids) # bring lost levels back in
    )
}
