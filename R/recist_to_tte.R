#' @export
recist_to_tte <- function(tbl_visits, dt_max = Inf) {
  dplyr::filter(tbl_visits, dt <= dt_max) %>%
    recist_to_tte_() %>%
    tibble::as_tibble()
}
