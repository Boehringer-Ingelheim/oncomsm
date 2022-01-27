#' Fully Exchangeable Model (FEM)
#'
#' todo
#'
#' @param mu mean of Gaussian prior for latent heterogeneity
#' @param tau standard deviation of Gaussian prior for latent heterogeneity
#'
#' @return An object of class FEM with prior information.
#'
#' @examples
#' this_function_needs_a_good_name(0, 0.05, 2, .1, 3, 1)
#'
#' @import purrr
#'
#' @export
this_function_needs_a_good_name <- function(
  logor_loc, logor_scale, # normal prior on log(OR)
  alpha_loc, alpha_scale, # normal prior for Weibull alpha parameter (shape)
  sigma_loc, sigma_scale # normal prior for Weibull scale parameter
) {
  assertthat::assert_that(
    all(purrr::map_lgl(
      c(logor_loc, alpha_loc, sigma_loc),
      is.numeric
    )),
    all(purrr::map_lgl(
      c(logor_loc, alpha_loc, sigma_loc),
      ~length(.) == 1
    )),
    msg = "location parameters must be single numeric values"
  )
  assertthat::assert_that(
    all(purrr::map_lgl(
      c(logor_scale, alpha_scale, sigma_scale),
      is.numeric
    )),
    all(purrr::map_lgl(
      c(logor_scale, alpha_scale, sigma_scale),
      ~length(.) == 1
    )),
    all(purrr::map_lgl(
      c(logor_scale, alpha_scale, sigma_scale),
      ~. > 1e-6
    )),
    msg = "scale parameters must be single, positive numeric value ( > 1e-6)"
  )
  res <- list(
    prior_params = list(
      logor = c(logor_loc, logor_scale),
      alpha = c(alpha_loc, alpha_scale),
      sigma = c(sigma_loc, sigma_scale)
    ),
    stanmodel = stanmodels$ThisModelNeedsAGoodName
  )
  # assign class information and return
  class(res) <- c("ThisModelNeedsAGoodName", class(res))
  res
}

#' @export
print.ThisModelNeedsAGoodName <- function(x, ...) {
  cat(sprintf(
    "ThisModelNeedsAGoodName<log(OR)~N(%.2f,%.2f),alpha~N(%.2f,%.2f)[1,Inf],sigma~N(%.2f,%.2f)[0,Inf]>",
    x$prior_params$logor[1], x$prior_params$logor[2],
    x$prior_params$alpha[1], x$prior_params$alpha[2],
    x$prior_params$sigma[1], x$prior_params$sigma[2]
  ))
}
