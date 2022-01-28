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

#' @include draw_samples.R
#' @export
draw_samples.ThisModelNeedsAGoodName <- function(model, n_to_be_recruited, nsim, data = NULL, seed = NULL, warmup = 1000L, ...) {
  assertthat::assert_that(nsim >= 10, msg = "nsim must be >= 10")
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  wrn <- list() # container for sampler warnings
  suppressWarnings(withCallingHandlers({
    res <- rstan::sampling(
      model$stanmodel,
      data = list(
        N_nonresponder = 0L,
        N_responder_exact = 0L,
        N_interval_censored = 0L,
        N_right_censored = 0L,
        N_to_be_recruited = n_to_be_recruited,
        t = numeric(),
        t1_right_censored = numeric(), #as.array(2),
        t1_interval_censored = numeric(),
        t2_interval_censored = numeric(),

        prior_logor_loc = model$prior_params$logor[1], prior_logor_scale = model$prior_params$logor[2],
        prior_alpha_loc = model$prior_params$alpha[1], prior_alpha_scale = model$prior_params$alpha[2],
        prior_sigma_loc = model$prior_params$sigma[1], prior_sigma_scale = model$prior_params$sigma[2]
      ),
      chains = 1L, cores = 1L,
      iter = warmup + nsim, warmup = warmup,
      seed = seed,
      verbose = FALSE, show_messages = FALSE, refresh = 0
    )
  },
  warning = function(w) wrn <<- c(wrn, list(w)) # log warnings
  ))
  return(list(sim = res, warnings = wrn))
}

