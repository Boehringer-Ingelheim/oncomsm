#' The 'bhmbasket.predict' package.
#'
#' @description This package implements methods to dynamically predict response status
#' of individuals in a basket trial with delayed binary response endpoint
#' (e.g. tumor response) at arbitrary time points.
#'
#' @docType package
#' @name bhmbasket.predict-package
#' @aliases bhmbasket.predict
#' @useDynLib bhmbasket.predict, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rlang .data
#'
#' @references
#' Stan Development Team (2021). RStan: the R interface to Stan. R package version 2.21.3. https://mc-stan.org
#'
NULL

## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
