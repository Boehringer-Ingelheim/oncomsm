#' @import purrr
#'
#' @export
this_function_needs_a_good_name <- function(
  logodds_loc, logodds_scale, # normal prior on log(odds)
  alpha_loc, alpha_scale, # normal prior for Weibull alpha parameter (shape)
  sigma_loc, sigma_scale, # normal prior for Weibull scale parameter
  dvisit
) {
  assertthat::assert_that(
    all(purrr::map_lgl(
      c(logodds_loc, alpha_loc, sigma_loc),
      is.numeric
    )),
    all(purrr::map_lgl(
      c(logodds_loc, alpha_loc, sigma_loc),
      ~length(.) == 1
    )),
    msg = "location parameters must be single numeric values"
  )
  assertthat::assert_that(
    all(purrr::map_lgl(
      c(logodds_scale, alpha_scale, sigma_scale),
      is.numeric
    )),
    all(purrr::map_lgl(
      c(logodds_scale, alpha_scale, sigma_scale),
      ~length(.) == 1
    )),
    all(purrr::map_lgl(
      c(logodds_scale, alpha_scale, sigma_scale),
      ~. > 1e-6
    )),
    msg = "scale parameters must be single, positive numeric value ( > 1e-6)"
  )
  res <- list(
    prior_params = list(
      logodds = c(logodds_loc, logodds_scale),
      alpha = c(alpha_loc, alpha_scale),
      sigma = c(sigma_loc, sigma_scale)
    ),
    dvisit = dvisit,
    stanmodel = stanmodels$ThisModelNeedsAGoodName
  )
  # assign class information and return
  class(res) <- c("ThisModelNeedsAGoodName", "TTEModel",class(res))
  res
}

#' @export
print.ThisModelNeedsAGoodName <- function(x, ...) {
  cat(sprintf(
    "ThisModelNeedsAGoodName<log(odds)~N(%.2f,%.2f),alpha~N(%.2f,%.2f)[1,Inf],sigma~N(%.2f,%.2f)[0,Inf],dvisit=%.2f>",
    x$prior_params$logodds[1], x$prior_params$logodds[2],
    x$prior_params$alpha[1], x$prior_params$alpha[2],
    x$prior_params$sigma[1], x$prior_params$sigma[2],
    x$dvisit
  ))
}

#' @include draw_samples.R
#' @export
draw_samples.ThisModelNeedsAGoodName <- function(
  model,
  data,
  nsim,
  seed = NULL, warmup = 1000L, verbose = FALSE, show_messages = FALSE, refresh = 0,
  return_raw_stan_output = FALSE,
  ...
) {
  # sampling with stan requires at least 2 samples to be drawn
  nsim_ <- max(2, nsim)
  # generate seed if none was specified
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  # convert factors to integers and save levels for mapping back
  data$subject_id <- factor(data$subject_id, levels = unique(data$subject_id))
  subject_id_lvls <- levels(data$subject_id)
  data$subject_id <- as.integer(data$subject_id)
  stan_data <- c(
    dvisit = model$dvisit,
    .tbl_to_stan_data(model, data),
    .prior_to_stan_data(model)
  )
  wrn <- list() # container for sampler warnings
  suppressWarnings(withCallingHandlers({
      res <- rstan::sampling(
        model$stanmodel,
        data = stan_data,
        chains = 1L, cores = 1L,
        iter = warmup + nsim, warmup = warmup,
        seed = seed,
        init = function() {
          list(
            alpha = 1,
            sigma = model$prior_params$sigma[1],
            theta = model$prior_params$logodds[1]
          )
        },
        verbose = verbose, show_messages = show_messages, refresh = refresh, ...
      )
    },
    warning = function(w) wrn <<- c(wrn, list(w)) # log warnings
  ))
  if (return_raw_stan_output) { # return stan samples directly
    attr(res, "stan_warnings") <- wrn
    return(res)
  } else {
    tbl_res <- dplyr::bind_cols(
      extract_col_from_stan(res, "subject_id"),
      extract_col_from_stan(res, "dt1") %>% select(-iter),
      extract_col_from_stan(res, "dt2") %>% select(-iter),
      extract_col_from_stan(res, "dt") %>% select(-iter)
    ) %>%
      filter(
        iter <= nsim # prune to intended number of samples
      )
    # map discrete variables back
    tbl_res$subject_id <- tbl_res$subject_id %>%
      factor(levels = 1:length(subject_id_lvls), labels = subject_id_lvls) %>%
      as.character()
    attr(tbl_res, "stan_warnings") <- wrn
    return(tbl_res)
  }
}

.prior_to_stan_data.ThisModelNeedsAGoodName <- function(model) {
  list(
    prior_logodds_loc = model$prior_params$logodds[1], prior_logodds_scale = model$prior_params$logodds[2],
    prior_alpha_loc = model$prior_params$alpha[1], prior_alpha_scale = model$prior_params$alpha[2],
    prior_sigma_loc = model$prior_params$sigma[1], prior_sigma_scale = model$prior_params$sigma[2]
  )
}

.tbl_to_stan_data.ThisModelNeedsAGoodName <- function(model, tbl_data) {
  lst_stan_data <- list()
  # interval censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$dt1) & is.finite(.data$dt2))
  lst_stan_data$N_A <- nrow(tmp)
  lst_stan_data$subject_id_A <- base::as.array(tmp$subject_id)
  lst_stan_data$dt1_A <- base::as.array(tmp$dt1) # as.array makes sure we have a vector, even if it is just a single number
  lst_stan_data$dt2_A <- base::as.array(tmp$dt2)
  # right censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$dt1) & !is.finite(.data$dt2))
  lst_stan_data$N_B <- nrow(tmp)
  lst_stan_data$subject_id_B <- base::as.array(tmp$subject_id)
  lst_stan_data$dt1_B <- base::as.array(tmp$dt1)
  # definite non-responders
  tmp <- dplyr::filter(tbl_data, is.infinite(.data$dt1))
  lst_stan_data$N_C <- nrow(tmp)
  lst_stan_data$subject_id_C <- base::as.array(tmp$subject_id)
  # new individuals
  tmp <- dplyr::filter(tbl_data, is.na(.data$dt1) & is.na(.data$dt2))
  lst_stan_data$N_D <- nrow(tmp)
  lst_stan_data$subject_id_D <- base::as.array(tmp$subject_id)
  # return
  return(lst_stan_data)
}
