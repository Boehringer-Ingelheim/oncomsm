#' Log-odds function
#'
#' computed the log odds of a probability.
#'
#' @param p numeric of probabilities
#'
#' @return log(p/(1-p))
#'
#' @export
logodds <- function(p) log(p / (1 - p))


# get unique identifiers that are reproducible
get_identifier <- function(n = 1, exclude = NULL) {
  res <- character(n)
  duplicate <- rep(TRUE, n)
  while (any(duplicate)) {
    res[duplicate] <- sprintf(
      "ID%08i",
      sample.int(1e7, sum(duplicate), replace = FALSE)
    )
    duplicate <- res %in% exclude
  }
  res
}
