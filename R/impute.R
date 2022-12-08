#' @description `impute()` samples visits for individuals in `data`
#' and potentially missing
#' individuals up to a maximum of `n_per_group` from the posterior
#' predictive distribution of the given model.
#'
#' @param now numeric, time since first visit in data if not last recorded
#' visit time
#' @template param-data-condition
#'
#' @examples
#' mdl <- create_srpmodel(A = define_srp_prior())
#' tbl <- tibble::tibble(
#'   subject_id = c("A1", "A1"),
#'   group_id = c("A", "A"),
#'   t = c(0, 1.5),
#'   state = c("stable", "stable")
#' )
#' impute(mdl, tbl, 1L, seed = 38L)
#'
#' @rdname sample_predictive
#' @export
impute <- function(model,
                   data,
                   nsim,
                   n_per_group = NULL,
                   sample = NULL,
                   p = NULL,
                   shape = NULL,
                   scale = NULL,
                   now = NULL,
                   seed = NULL,
                   nsim_parameters = 1000L,
                   warmup_parameters = 250L,
                   nuts_control = list(),
                   as_mstate = FALSE,
                   ...) {
  checkmate::check_class(model, classes = c("srpmodel", "list"))
  if (!is.null(seed)) {
    set.seed(seed)
  }
  recruitment_rate <- model$recruitment_rate
  group_ids <- model$group_id
  if (is.null(now)) { # convert to visits and take the last time point
    now <- max(data$t)
  }
  if (is.null(sample)) {
    sample <- sample_posterior(model,
                               data = data, now = now, seed = seed,
                               warmup = warmup_parameters,
                               nsim = nsim_parameters,
                               nuts_control = nuts_control)
  }
  if (is.null(n_per_group)) {
    # no new individuals
    n_per_group <- data %>%
      select("group_id", "subject_id") %>%
      distinct() %>%
      pull("group_id") %>%
      table() %>%
      .[group_ids] %>%
      as.numeric()
  } else {
    if (is.null(recruitment_rate)) {
      stop("recruitment_rates must be specified") # nocov
    }
  }
  tbl_to_recruit <- tibble(
    subject_id = character(0L),
    group_id = character(0L),
    t = numeric(0L),
    state = character(0L)
  )
  for (i in seq_along(group_ids)) {
    n_recruited <- data %>%
      filter(.data$group_id == group_ids[i]) %>%
      pull("subject_id") %>%
      unique() %>%
      length()
    n_to_be_recruited <- n_per_group[i] - n_recruited
    if (n_to_be_recruited < 0) { # nocov start
      stop("data contains more individuals than specified in n_per_group") # nolint
    } # nocov end
    ids_to_exclude <- c(tbl_to_recruit$subject_id, unique(data$subject_id))
    if (n_to_be_recruited > 0) {
      subject_ids <- get_identifier(n = n_to_be_recruited,
                                    exclude = ids_to_exclude)
      recruitment_times <- now +
        cumsum(stats::rexp(n_to_be_recruited, rate = recruitment_rate[i]))
      tbl_to_recruit <- bind_rows(tbl_to_recruit, tibble(
        subject_id = subject_ids,
        group_id = group_ids[i],
        t = recruitment_times,
        state = "stable" # first visits are always stable
      ))
    }
  }
  tbl_data <- bind_rows(data, tbl_to_recruit)
  # do actual imputation
  data <- tbl_data
  parameter_sample <- sample
  n_groups <- length(model$group_id)
  # make sure that either parameter sample or p, scale, shape are given
  if (!is.null(parameter_sample)) {
    stopifnot(isa(parameter_sample, "stanfit"))
    n_params_sample <- parameter_sample@sim$iter - parameter_sample@sim$warmup
  } else {
    # p, scale, and shape must be given
    if (is.null(p) || is.null(shape) || is.null(scale)) { # nocov start
      stop("if no parameter sample is given all of p, scale, shape must be given") # nolint
    }
    if (length(p) != n_groups) stop()
    if (all(dim(shape) != c(n_groups, 3))) stop()
    if (all(dim(scale) != c(n_groups, 3))) stop()
  } # nocov end
  # extract subject and group id levels for conversion to and back from integer
  subject_id_levels <- unique(as.character(data$subject_id))
  group_id_levels <- model$group_id # important to maintain ordering
  if (is.null(p)) {
    # extract parameter arrays from stanfit object
    # p[i,j] is probability for ith sample for jth group
    p <- rstan::extract(parameter_sample, "p")[[1]]
  } else {
    # expand fixed parameters to same format as rstan parameters
    p <- t(array(p, dim = c(n_groups, n_params_sample)))
  }
  if (is.null(shape)) {
    shape <- rstan::extract(parameter_sample, "shape")[[1]]
  } else {
    # need to add one dimension (iterations to fit format)
    shape <- aperm(
      array(shape, dim = c(n_groups, 3, n_params_sample)),
      c(3, 1, 2)
    )
  }
  if (is.null(scale)) {
    # scale[i,j,k] is the scale value for ith sample in jth group for
    # k-th transition (k being 1. stable-> response, 2. stable -> progression,
    # 3. response -> progression)
    scale <- rstan::extract(parameter_sample, "scale")[[1]]
  } else {
    # need to add one dimension (iterations to fit format)
    scale <- aperm(
      array(scale, dim = c(n_groups, 3, n_params_sample)),
      c(3, 1, 2)
    )
  }
  # sorting the samples and changing type to integer for groups and subj id
  data <- data %>%
    arrange(.data$subject_id, .data$t) %>%
    mutate( # convert group_id to properly ordered factor
      group_id = factor(.data$group_id, levels = group_id_levels)
    )
  idx <- sample(seq_len(n_params_sample), size = nsim, replace = TRUE)
  sample_once <- function(iter) {
    # extract a set of parameters
    response_probabilities <- as.array(p[idx[iter], ])
    shapes <- matrix(shape[idx[iter], , ], ncol = 3)
    scales <- matrix(scale[idx[iter], , ], ncol = 3)
    # sample using C++ implementation
    res <- impute_srp_model(data, response_probabilities, shapes, scales,
                            visit_spacing = model$visit_spacing,
                            max_time = model$maximal_time
      ) %>%
      as_tibble() %>%
      mutate(
        group_id = as.character(.data$group_id)
      )
    if (as_mstate) {
      res <- visits_to_mstate(res, model)
    }
    res <- mutate(res, iter = as.integer(iter))
    return(res)
  }
  res <- purrr::map_df(1:nsim, sample_once)
  return(res)
}
