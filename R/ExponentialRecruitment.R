#' @export
exponential_recruitment_model <- function(
  theta_loc, theta_scale, dt_max
) {
  res <- list(
    prior_params = list(
      theta = c(theta_loc, theta_scale)
    ),
    dt_max = dt_max,
    stanmodel = stanmodels$ExponentialRecruitment
  )
  # assign class information and return
  class(res) <- c("ExponentialRecruitment", "RecruitmentModel", class(res))
  res
}

#' @export
print.ExponentialRecruitment <- function(x, ...) {
  cat(sprintf(
    "ExponentialRecruitment<log(rate)~N(%.2f,%.2f),wait<=%.2f>",
    x$prior_params$theta[1], x$prior_params$theta[2],
    x$dt_max
  ))
}

#' @include draw_samples.R
#' @include util.R
#' @export
draw_samples.ExponentialRecruitment <- function(
  model,
  data,
  now,
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
  wrn <- list() # container for sampler warnings
  # convert factors to integers and save levels for mapping back
  data$subject_id <- factor(data$subject_id, levels = unique(data$subject_id))
  subject_id_lvls <- levels(data$subject_id)
  data$subject_id <- as.integer(data$subject_id)
  stan_data <- c(
    list(
      now = now,
      dt_max = model$dt_max
    ),
    .tbl_to_stan_data(model, data),
    .prior_to_stan_data(model)
  )
  suppressWarnings(withCallingHandlers({
      res <- rstan::sampling(
        model$stanmodel,
        data = stan_data,
        chains = 1L, cores = 1L,
        iter = warmup + nsim, warmup = warmup,
        seed = seed,
        init = function() { # set params to mean values
          list(
            theta = model$prior_params$theta[1]
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
        extract_col_from_stan(res, "t") %>% select(-iter)
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

.prior_to_stan_data.ExponentialRecruitment <- function(model) {
  list(
    prior_theta_loc = model$prior_params$theta[1],
    prior_theta_scale = model$prior_params$theta[2]
  )
}

.tbl_to_stan_data.ExponentialRecruitment <- function(model, tbl_data) {
  lst_stan_data <- list()
  # observed data
  tmp <- dplyr::filter(tbl_data, !is.na(t))
  lst_stan_data$N_old <- nrow(tmp)
  lst_stan_data$subject_id_old <- base::as.array(tmp$subject_id)
  lst_stan_data$t_old <- base::as.array(tmp$t)
  # unobserved data
  tmp <- dplyr::filter(tbl_data, is.na(t))
  lst_stan_data$N_new <- nrow(tmp)
  lst_stan_data$subject_id_new <- base::as.array(tmp$subject_id)
  return(lst_stan_data)
}
