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
  sigma_loc, sigma_scale, # normal prior for Weibull scale parameter
  lograte_loc, lograte_scale # normal prior for (log) exponential rate
) {
  assertthat::assert_that(
    all(purrr::map_lgl(
      c(logor_loc, alpha_loc, sigma_loc, lograte_loc),
      is.numeric
    )),
    all(purrr::map_lgl(
      c(logor_loc, alpha_loc, sigma_loc, lograte_loc),
      ~length(.) == 1
    )),
    msg = "location parameters must be single numeric values"
  )
  assertthat::assert_that(
    all(purrr::map_lgl(
      c(logor_scale, alpha_scale, sigma_scale, lograte_scale),
      is.numeric
    )),
    all(purrr::map_lgl(
      c(logor_scale, alpha_scale, sigma_scale, lograte_scale),
      ~length(.) == 1
    )),
    all(purrr::map_lgl(
      c(logor_scale, alpha_scale, sigma_scale, lograte_scale),
      ~. > 1e-6
    )),
    msg = "scale parameters must be single, positive numeric value ( > 1e-6)"
  )
  res <- list(
    prior_params = list(
      logor = c(logor_loc, logor_scale),
      alpha = c(alpha_loc, alpha_scale),
      sigma = c(sigma_loc, sigma_scale),
      lograte = c(lograte_loc, lograte_scale)
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
    "ThisModelNeedsAGoodName<log(OR)~N(%.2f,%.2f),alpha~N(%.2f,%.2f)[1,Inf],sigma~N(%.2f,%.2f)[0,Inf],log(lambda)~N(%.2f,%.2f)>",
    x$prior_params$logor[1], x$prior_params$logor[2],
    x$prior_params$alpha[1], x$prior_params$alpha[2],
    x$prior_params$sigma[1], x$prior_params$sigma[2],
    x$prior_params$lograte[1], x$prior_params$lograte[2]
  ))
}

#' @include draw_samples.R
#' @export
draw_samples.ThisModelNeedsAGoodName <- function(
  model, n_to_be_recruited, nsim, now = 0, data = NULL, seed = NULL, warmup = 1000L,
  verbose = FALSE, show_messages = FALSE, refresh = 0, ...
) {
  # sample using stan
  assertthat::assert_that(nsim >= 10, msg = "nsim must be >= 10")
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  wrn <- list() # container for sampler warnings
  stan_data <- c(
    list(now = now),
    .tbl_to_stan_data(model, data, n_to_be_recruited),
    .prior_to_stan_data(model)
  )
  suppressWarnings(withCallingHandlers({
      res <- rstan::sampling(
        model$stanmodel,
        data = stan_data,
        chains = 1L, cores = 1L,
        iter = warmup + nsim, warmup = warmup,
        seed = seed,
        init = function() {
          list(alpha = 1,
               sigma = model$prior_params$sigma[1],
               logor_response = 0,
               lograte = model$prior_params$lograte[1]
          )
        },
        verbose = verbose, show_messages = show_messages, refresh = refresh, ...
      )
    },
    warning = function(w) wrn <<- c(wrn, list(w)) # log warnings
  ))
  # post-process sampled data
  tbl_res <- tibble::tibble(
    iter = integer(0L),
    subject_id = integer(0L),
    dt_response = numeric(0L)
  )
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
        tidyr::pivot_longer(-iter, names_to = "subject_id", values_to = "dt_response") %>%
        dplyr::mutate(
          iter = as.integer(.data$iter),
          subject_id = as.integer(.data$subject_id)
        ) %>%
        select(iter, subject_id, dt_response)
    )
  }
  # add t = Inf for definite non-responders
  tbl_res <- dplyr::bind_rows(
    tbl_res,
    tidyr::expand_grid(iter = as.integer(1:nsim), subject_id = stan_data$id_nonresponder, dt_response = Inf)
  )
  # add recruitment time for to-be-recruited
  tmp <- rstan::extract(res, "t_recruitment_to_be_recruited")[["t_recruitment_to_be_recruited"]]
  if (!is.null(tmp)) {
    tbl_res <- left_join(
      tbl_res,
      tmp %>%
        {colnames(.) <- stan_data[["id_to_be_recruited"]]; rownames(.) <- 1:nrow(.); .} %>%
        tibble::as_tibble(rownames = NA) %>%
        tibble::rownames_to_column("iter") %>%
        tidyr::pivot_longer(-iter, names_to = "subject_id", values_to = "t_recruitment_sampled") %>%
        dplyr::mutate(
          iter = as.integer(.data$iter),
          subject_id = as.integer(.data$subject_id)
        ) %>%
        select(iter, subject_id, t_recruitment_sampled),
      by = c("subject_id", "iter")
    )
  }
  # merge with original data
  if (!is.null(data))  {
    tbl_res <- dplyr::full_join(data, tbl_res, by = "subject_id") %>%
      mutate(
        t_recruitment = if_else(is.na(t_recruitment), t_recruitment_sampled, t_recruitment)
      ) %>%
      select(-t_recruitment_sampled)
  } else {
    tbl_res <- rename(tbl_res, t_recruitment = t_recruitment_sampled)
  }
  attr(tbl_res, "stan_warnings") <- wrn
  tbl_res
}

.prior_to_stan_data.ThisModelNeedsAGoodName <- function(model) {
  list(
    prior_logor_loc = model$prior_params$logor[1], prior_logor_scale = model$prior_params$logor[2],
    prior_alpha_loc = model$prior_params$alpha[1], prior_alpha_scale = model$prior_params$alpha[2],
    prior_sigma_loc = model$prior_params$sigma[1], prior_sigma_scale = model$prior_params$sigma[2],
    prior_lograte_loc = model$prior_params$lograte[1], prior_lograte_scale = model$prior_params$lograte[2]
  )
}

.tbl_to_stan_data.ThisModelNeedsAGoodName <- function(model, tbl_data, n_additional) {
  tbl_data <- check_input_data(tbl_data)
  data <- list() # container for data to pass to stan
  # interval censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$dt_response_min) & is.finite(.data$dt_response_max))
  data$N_interval_censored <- nrow(tmp)
  data$t_recruitment_interval_censored <- base::as.array(tmp$t_recruitment)
  data$t1_interval_censored <- base::as.array(tmp$dt_response_min) # as.array makes sure we have a vector, even if it is just a single number
  data$t2_interval_censored <- base::as.array(tmp$dt_response_max)
  data$id_interval_censored <- base::as.array(tmp$subject_id)
  # right censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$dt_response_min) & !is.finite(.data$dt_response_max))
  data$N_right_censored <- nrow(tmp)
  data$t_recruitment_right_censored <- base::as.array(tmp$t_recruitment)
  data$t1_right_censored <- base::as.array(tmp$dt_response_min)
  data$id_right_censored <- base::as.array(tmp$subject_id)
  # definite non-responders
  tmp <- dplyr::filter(tbl_data, is.infinite(.data$dt_response_min))
  data$N_nonresponder <- nrow(tmp)
  data$t_recruitment_nonresponder <- base::as.array(tmp$t_recruitment)
  data$id_nonresponder <- base::as.array(tmp$subject_id)
  # no data at all
  tmp <- dplyr::filter(tbl_data, is.na(.data$dt_response_min) & is.na(.data$dt_response_max))
  data$N_to_be_recruited <- nrow(tmp) + n_additional
  data$id_to_be_recruited <- base::as.array(tmp$subject_id)
  if (n_additional > 0) {
    new_id_start <- max(c(1, tbl_data$subject_id))
    data$id_to_be_recruited <- c(data$id_to_be_recruited, new_id_start:(new_id_start - 1 + n_additional))
  }
  return(data)
}
