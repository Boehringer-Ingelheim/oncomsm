#' Log-odds function
#'
#' computed the log odds of a probability.
#'
#' @param p numeric of probabilities
#'
#' @return log(p/(1-p))
#'
#' @export
logodds <- function(p) log( p / (1 - p) )


# helper function to broadcast paramters
.repifnot <- function(x, m) {
  if (length(x) == 1)
    return(rep(x, m))
  stopifnot(length(x) == m)
  return(x)
}
