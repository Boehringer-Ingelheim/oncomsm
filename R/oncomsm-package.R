#' The oncomsm package.
#'
#' @description This package implements methods to dynamically predict response
#' and progression of individuals in early oncology trials using parametric
#' multi-state models and Bayesian inference.
#' This allows the dynamic computation of Probability of Success for a wide
#' range of success criteria.
#' The inference is implemented using stan.
#'
#' @docType package
#' @name oncomsm-package
#' @aliases oncomsm
#' @useDynLib oncomsm, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import dplyr
#' @importFrom rstan sampling
#' @importFrom rlang .data
#'
#' @references
#' Stan Development Team (2021). RStan: the R interface to Stan.
#' R package version 2.21.3. https://mc-stan.org
#'
NULL

## usethis namespace: start
#' @importFrom magrittr %>%
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
