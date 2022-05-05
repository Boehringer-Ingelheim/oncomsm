#' Create an instance of the mixture cure-rate model for binary tumor response
#'
#' This model assumes that the response rates across groups are independent.
#'
#' @param group_id character vector of group ids
#' @param logodds_mean numeric vector with the means of the normal prior on the logodds of a response per group
#' @param logodds_sd numeric vector with the standard deviations of the normal prior on the logodds of a response per group
#' @param logodds_min numeric vector with the minimal logodds per group
#' @param logodds_max numeric vector with the maximal logodds per group
#' @param shape_mean numeric vector with the means of the normal prior on the shape parameter of the Weibull distribution for time to response
#' @param shape_sd numeric vector with the standard deviations of the normal prior on the shape parameter of the Weibull distribution for time to response
#' @param scale_mean numeric vector with the means of the normal prior on the scale parameter of the Weibull distribution for time to response
#' @param scale_sd numeric vector with the standard deviations of the normal prior on the scale parameter of the Weibull distribution for time to response
#' @param visit_spacing vector of deterministic spacing between future visits (in months)
#'
#' @return An object of class "IndependentMixtureCureRateModel" holding all relevant
#' prior information.
#'
#' @rdname IndependentMixtureCureRateModel
#'
#' @export
independent_mixture_cure_rate_model <- function(
  group_id,
  logodds_mean, logodds_sd, logodds_min = rep(-Inf, length(group_id)), logodds_max = rep(Inf, length(group_id)), # (truncated) normal prior on log(odds)
  shape_mean, shape_sd, # normal prior for Weibull alpha parameter (shape)
  scale_mean, scale_sd, # normal prior for Weibull scale parameter
  visit_spacing # (deterministic) spacing between (future) visits
) {
  res <- as.list(environment()) # store all input parameters
  res$group_id <- NULL
  attr(res, "group_id") <- group_id
  attr(res, "stanmodel") <- stanmodels$IndependentMixtureCureRateModel
  # assign class information and return
  class(res) <- c("IndependentMixtureCureRateModel", "TTEModel",class(res))
  res
}



#' @include draw_samples.R util.R
#'
#' @describeIn IndependentMixtureCureRateModel draw samples
#' @export
draw_samples.IndependentMixtureCureRateModel <- function(
  model,
  data = NULL,
  n = NULL,
  nsim,
  seed = NULL, warmup = 1000L, verbose = FALSE, show_messages = FALSE, refresh = 0,
  return_raw_stan_output = FALSE,
  ...
) {
  if (is.null(data)) {
    if (is.null(n))
      stop("either data or n must be specified")
    group_id <- attr(model, "group_id")
    if (length(n) != length(group_id))
      stop("number of group sizes (n) not consistent with length of group_id of model")
    # if only n is specified, create empty data set of appropriate sample size
    data <- tibble::tibble(
      group_id = rep(attr(model, "group_id"), times = n), # pull group ids directly from model
      subject_id = as.integer(1:sum(n)),
      dt1 = NA_real_,
      dt2 = NA_real_
    )
  }
  # for numerical stability, avoid dt1 == 0 or dt2 = dt1
  data <- data %>% dplyr::mutate(
    dt1 = pmax(1e-3, .data$dt1),
    dt2 = pmax(.data$dt1 + 1e-3, .data$dt2)
  )
  # sampling with stan requires at least 2 samples to be drawn oO
  nsim_ <- max(2, nsim)
  # generate seed if none was specified
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  # convert factors to integers and save levels for mapping back
  data$group_id <- factor(data$group_id, levels = unique(data$group_id))
  group_id_lvls <- levels(data$group_id)
  # make sure group ids in data are compatible with model
  assertthat::assert_that(all(
    group_id_lvls %in% attr(model, "group_id")
  ))
  data$group_id <- as.integer(data$group_id)
  data$subject_id <- factor(data$subject_id, levels = unique(data$subject_id))
  subject_id_lvls <- levels(data$subject_id)
  data$subject_id <- as.integer(data$subject_id)
  # combine prior information with data for stan
  stan_data <- c(
    as.list(model),
    .tbl_to_stan_data(model, data)
  )
  wrn <- list() # container for sampler warnings
  suppressWarnings(withCallingHandlers({
    res <- rstan::sampling(
      attr(model, "stanmodel"),
      data = stan_data,
      chains = 1L, cores = 1L,
      iter = warmup + nsim_, warmup = warmup,
      seed = seed,
      init = function() {
        list( # initialize to parameter means
          shape = model$shape_mean,
          scale = model$scale_mean,
          theta = model$logodds_mean
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
  }
  tbl_res <- dplyr::bind_cols(
    extract_col_from_stan(res, "group_id"),
    extract_col_from_stan(res, "subject_id") %>% dplyr::select(-.data$iter),
    extract_col_from_stan(res, "dt1") %>% dplyr::select(-.data$iter),
    extract_col_from_stan(res, "dt2") %>% dplyr::select(-.data$iter),
    extract_col_from_stan(res, "dt") %>% dplyr::select(-.data$iter)
  ) %>%
    dplyr::filter(
      .data$iter <= nsim # prune to intended number of samples
    )
  # map discrete variables back to character
  tbl_res$group_id <- tbl_res$group_id %>%
    factor(levels = 1:length(group_id_lvls), labels = group_id_lvls) %>%
    as.character()
  tbl_res$subject_id <- tbl_res$subject_id %>%
    factor(levels = 1:length(subject_id_lvls), labels = subject_id_lvls) %>%
    as.character()
  attr(tbl_res, "stan_warnings") <- wrn
  return(tbl_res)
}



.tbl_to_stan_data.IndependentMixtureCureRateModel <- function(model, tbl_data) {
  lst_stan_data <- list(
    M_groups = length(unique(tbl_data$group_id))
  )
  # interval censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$dt1) & is.finite(.data$dt2))
  lst_stan_data$N_A <- nrow(tmp)
  lst_stan_data$group_id_A <- base::as.array(tmp$group_id)
  lst_stan_data$subject_id_A <- base::as.array(tmp$subject_id)
  lst_stan_data$dt1_A <- base::as.array(tmp$dt1) # as.array makes sure we have a vector, even if it is just a single number
  lst_stan_data$dt2_A <- base::as.array(tmp$dt2)
  # right censored data
  tmp <- dplyr::filter(tbl_data, is.finite(.data$dt1) & !is.finite(.data$dt2))
  lst_stan_data$N_B <- nrow(tmp)
  lst_stan_data$group_id_B <- base::as.array(tmp$group_id)
  lst_stan_data$subject_id_B <- base::as.array(tmp$subject_id)
  lst_stan_data$dt1_B <- base::as.array(tmp$dt1)
  # definite non-responders
  tmp <- dplyr::filter(tbl_data, is.infinite(.data$dt1))
  lst_stan_data$N_C <- nrow(tmp)
  lst_stan_data$group_id_C <- base::as.array(tmp$group_id)
  lst_stan_data$subject_id_C <- base::as.array(tmp$subject_id)
  # new individuals
  tmp <- dplyr::filter(tbl_data, is.na(.data$dt1) & is.na(.data$dt2))
  lst_stan_data$N_D <- nrow(tmp)
  lst_stan_data$group_id_D <- base::as.array(tmp$group_id)
  lst_stan_data$subject_id_D <- base::as.array(tmp$subject_id)
  # return
  return(lst_stan_data)
}



#' @describeIn IndependentMixtureCureRateModel plot prior assumptions
#'
#' @export
plot.IndependentMixtureCureRateModel <- function(
  x, n = rep(1L, length(attr(model, "group_id"))), nsim = 1e4, seed = NULL, ...
) {

  # prior sample from response probability
  stan_res <- draw_samples(x, n = n, nsim = nsim, seed = seed, return_raw_stan_output = TRUE, ...) %>%
    rstan::extract("p")

  p1 <- stan_res$p %>%
    as.matrix() %>%
    {colnames(.) <- attr(x, "group_id"); .} %>%
    as_tibble(rownames = "iter") %>%
    tidyr::pivot_longer(-.data$iter, names_to = "group_id", values_to = "response probability") %>%
    ggplot() +
    aes(`response probability`) +
    geom_histogram(bins = 25) +
    scale_x_continuous(breaks = seq(0, 1, by = .2), limits = c(0, 1), oob = scales::oob_keep) +
    facet_wrap(~group_id) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  p2 <- draw_samples(x, n = n, nsim = nsim, seed = seed, ...) %>%
    filter(is.finite(dt)) %>%
    ggplot() +
    aes(dt) +
    geom_histogram(bins = 25) +
    labs(x = "time to response (since first visit)") +
    facet_wrap(~group_id, scales = "free_y") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  p1 + p2

}
