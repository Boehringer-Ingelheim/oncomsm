#' Create an instance of a Poisson recruitment model
#'
#' @rdname IndependentPoissonRecruitmentModel
#'
#' @export
independent_poisson_recruitment_model <- function(
  group_id,
  log_monthly_rate_mean, log_monthly_rate_sd,
  maximal_recruitment_interval
) {
  res <- as.list(environment()) # store all input parameters
  res$group_id <- NULL
  res <- lapply(res, base::as.array)
  attr(res, "group_id") <- group_id
  attr(res, "stanmodel") <- stanmodels$IndependentPoissonRecruitmentModel
  # assign class information and return
  class(res) <- c("IndependentPoissonRecruitmentModel", "RecruitmentModel", class(res))
  res
}



#' @include draw_samples.R util.R
#'
#' @describeIn IndependentPoissonRecruitmentModel draw samples
#' @export
draw_samples.IndependentPoissonRecruitmentModel <- function(
  model,
  data = NULL,
  n = NULL,
  now,
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
      t_recruitment = NA_real_
    )
  }
  # sampling with stan requires at least 2 samples to be drawn
  nsim_ <- max(2, nsim)
  # generate seed if none was specified
  if (is.null(seed))
    seed <- sample.int(.Machine$integer.max, 1)
  wrn <- list() # container for sampler warnings
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
  stan_data <- c(
    as.list(model),
    .tbl_to_stan_data(model, data),
    list(
      now = now
    )
  )
  suppressWarnings(withCallingHandlers({
    res <- rstan::sampling(
      attr(model, "stanmodel"),
      data = stan_data,
      chains = 1L, cores = 1L,
      iter = warmup + nsim_, warmup = warmup,
      seed = seed,
      init = function() { # set params to mean values
        list(
          theta = model$monthly_rate_mean
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
    extract_col_from_stan(res, "t") %>% dplyr::select(-.data$iter)
  ) %>%
    dplyr::filter(
      .data$iter <= nsim # prune to intended number of samples
    )
  # map discrete variables back
  tbl_res$group_id <- tbl_res$group_id %>%
    factor(levels = 1:length(group_id_lvls), labels = group_id_lvls) %>%
    as.character()
  tbl_res$subject_id <- tbl_res$subject_id %>%
    factor(levels = 1:length(subject_id_lvls), labels = subject_id_lvls) %>%
    as.character()
  attr(tbl_res, "stan_warnings") <- wrn
  return(tbl_res)
}



.tbl_to_stan_data.IndependentPoissonRecruitmentModel <- function(model, tbl_data) {
  lst_stan_data <- list(
    M_groups = length(unique(tbl_data$group_id))
  )
  # observed data
  tmp <- dplyr::filter(tbl_data, !is.na(.data$t_recruitment))
  lst_stan_data$N_old <- nrow(tmp)
  lst_stan_data$group_id_old <- base::as.array(tmp$group_id)
  lst_stan_data$subject_id_old <- base::as.array(tmp$subject_id)
  lst_stan_data$t_old <- base::as.array(tmp$t_recruitment)
  # unobserved data
  tmp <- dplyr::filter(tbl_data, is.na(.data$t_recruitment))
  lst_stan_data$N_new <- nrow(tmp)
  lst_stan_data$group_id_new <- base::as.array(tmp$group_id)
  lst_stan_data$subject_id_new <- base::as.array(tmp$subject_id)
  return(lst_stan_data)
}



#' @describeIn IndependentPoissonRecruitmentModel plot prior assumptions
#'
#' @export
plot.IndependentPoissonRecruitmentModel <- function(
  x, n, nsim = 1e4, seed = NULL, ...
) {

  # prior sample from response probability
  stan_res <- draw_samples(
    x, n = n, nsim = nsim, seed = seed, return_raw_stan_output = TRUE, ...
  ) %>%
    rstan::extract("rate")

  p1 <- stan_res$rate %>%
    as.matrix() %>%
    {colnames(.) <- attr(x, "group_id"); .} %>%
    as_tibble(rownames = "iter") %>%
    tidyr::pivot_longer(-.data$iter, names_to = "group_id", values_to = "recruitment rate [monthly]") %>%
    ggplot() +
    aes(`recruitment rate [monthly]`) +
    geom_histogram(bins = 25) +
    facet_wrap(~group_id) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )

  p2 <- draw_samples(x, n = n, nsim = nsim, seed = seed, ...) %>%
    group_by(group_id, .data$iter) %>%
    summarize(`end time [months]` = max(t)) %>%
    ggplot() +
    aes(`end time [months]`) +
    geom_histogram(bins = 25) +
    facet_wrap(~group_id, scales = "free_y") +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )

  p1 + p2

}
