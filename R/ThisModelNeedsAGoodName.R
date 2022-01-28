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
draw_samples.ThisModelNeedsAGoodName <- function(
  model, n_to_be_recruited, nsim, data = NULL, seed = NULL, warmup = 1000L,
  verbose = FALSE, show_messages = FALSE, refresh = 0, ...
) {
  # sample using stan
  assertthat::assert_that(nsim >= 10, msg = "nsim must be >= 10")
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  wrn <- list() # container for sampler warnings
  stan_data <- c(.tbl_to_stan_data(model, data, n_to_be_recruited), .prior_to_stan_data(model))
  suppressWarnings(withCallingHandlers({
      res <- rstan::sampling(
        model$stanmodel,
        data = stan_data,
        chains = 1L, cores = 1L,
        iter = warmup + nsim, warmup = warmup,
        seed = seed,
        init = function() {
          list(alpha = 1, sigma = model$prior_params$sigma[1], logor_response = 0)
        },
        verbose = verbose, show_messages = show_messages, refresh = refresh, ...
      )
    },
    warning = function(w) wrn <<- c(wrn, list(w)) # log warnings
  ))
  # post-process sampled data
  tbl_res <- tibble::tibble(id = integer(0L), iter = integer(0L), t = numeric(0L))
  for (type in c("interval_censored", "right_censored", "to_be_recruited") ) {
    var_name <- sprintf("response_time_%s", type)
    tmp <- rstan::extract(res, var_name)[[var_name]]
    if (is.null(tmp))
      next
    tbl_res <- dplyr::bind_rows(
      tbl_res,
      tmp %>%
        {colnames(.) <- stan_data[[sprintf("id_%s", type)]]; rownames(.) <- 1:nrow(.); .} %>%
        tibble::as_tibble(rownames = NA) %>%
        tibble::rownames_to_column("iter") %>%
        tidyr::pivot_longer(-iter, names_to = "id", values_to = "t") %>%
        dplyr::mutate(
          iter = as.integer(.data$iter),
          id = as.integer(.data$id)
        )
    )
  }
  # add t = Inf for definite non-responders
  tbl_res <- dplyr::bind_rows(
    tbl_res,
    tidyr::expand_grid(id = stan_data$id_nonresponder, iter = as.integer(1:nsim), t = Inf)
  )
  tbl_res <- dplyr::left_join(tbl_data, tbl_res, by = "id")
  attr(tbl_res, "stan_warnings") <- wrn
  tbl_res
}

.prior_to_stan_data.ThisModelNeedsAGoodName <- function(model) {
  list(
    prior_logor_loc = model$prior_params$logor[1], prior_logor_scale = model$prior_params$logor[2],
    prior_alpha_loc = model$prior_params$alpha[1], prior_alpha_scale = model$prior_params$alpha[2],
    prior_sigma_loc = model$prior_params$sigma[1], prior_sigma_scale = model$prior_params$sigma[2]
  )
}

.tbl_to_stan_data.ThisModelNeedsAGoodName <- function(model, tbl_data, n_additional) {
  if (is.null(tbl_data))
    tbl_data <- tibble::tibble(
      id = integer(0L),
      t1 = numeric(0L),
      t2 = numeric(0L)
    )
  assertthat::assert_that(
    ncol(tbl_data) == 3,
    all(colnames(tbl_data) == c("id", "t1", "t2")),
    all(purrr::map_chr(tbl_data, class) == c("integer", "numeric", "numeric"))
  )
  data <- list() # container for data to pass to stan
  # interval censored data
  tmp <- dplyr::filter(tbl_data, is.finite(t1) & is.finite(.data$t2))
  data$N_interval_censored <- nrow(tmp)
  data$t1_interval_censored <- base::as.array(tmp$t1) # as.array makes sure we have a vector, even if it is just a single number
  data$t2_interval_censored <- base::as.array(tmp$t2)
  data$id_interval_censored <- base::as.array(tmp$id)
  # right censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$t1) & !is.finite(.data$t2))
  data$N_right_censored <- nrow(tmp)
  data$t1_right_censored <- base::as.array(tmp$t1)
  data$id_right_censored <- base::as.array(tmp$id)
  # definite non-responders
  tmp <- dplyr::filter(tbl_data, is.infinite(.data$t1))
  data$N_nonresponder <- nrow(tmp)
  data$id_nonresponder <- base::as.array(tmp$id)
  # no data at all
  tmp <- dplyr::filter(tbl_data, is.na(.data$t1) & is.na(.data$t2))
  data$N_to_be_recruited <- nrow(tmp) + n_additional
  data$id_to_be_recruited <- base::as.array(tmp$id)
  if (n_additional > 0) {
    new_id_start <- max(c(1, tbl_data$id))
    data$id_to_be_recruited <- c(data$id_to_be_recruited, new_id_start:(new_id_start + n_additional))
  }
  return(data)
}
