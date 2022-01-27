#' Simulate data
#'
#' simulate exact response time points from a Weibull mixture cure rate model.
#'
#' @param model the model to simulate from
#' @param data **optional** data to condition on
#'
#' @export
simulate_data <- function(model, n_to_be_recruited, nsim, data = NULL, seed = 42L, warmup = 1000L, ...) {
  assertthat::assert_that(nsim >= 10, msg = "nsim must be >= 10")
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
