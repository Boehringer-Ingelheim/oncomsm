#' A Stable-Response-Progression Model
#'
#' Create a new instance of an SRP model
#'
#' TODO
#'
#' @param group_id a character vector with the group ids, these are used to
#'   check compatibility of data later
#' @param logodds_mean a vector with the means of the (truncated) normal priors
#'   on the log-odds of the response probability
#' @param logodds_sd a vector with the standard deviations of the (truncated)
#'   normal priors on the log-odds of the response probability
#' @param median_time_to_next_event_mean a matrix with the means of the
#'   (truncated) normal priors on the median time to next event for each of the
#'   Weibull transition probabilities, the (i,j)-th entry is the i-th group
#'   median time to next event for transition j
#'   (1=stable-response, 2=stable-progression, 3=response-progression)
#' @param median_time_to_next_event_sd a matrix with the standard deviations of
#'   the (truncated) normal priors on the median time to next event for each of
#'   the Weibull transition probabilities, the (i,j)-th entry is the i-th group
#'   median time to next event for transition j
#'   (1=stable-response, 2=stable-progression, 3=response-progression)
#' @param visit_spacing vector of time differences between visits per group
#' @param logodds_min lower boundary on the log-odds per group
#' @param logodds_max upper boundary on the log-odds per group
#' @param shape_min matrix of lower boundaries of the uniform prior of the
#'   Weibull distribution per group/transition
#' @param shape_max matrix of upper boundaries of the uniform prior of the
#'   Weibull distribution per group/transition
#'
#' @name srp_model
#' @aliases create_srp_model
#' @seealso [Model]
#'
#' @export
create_srp_model <- function(
  group_id,
  logodds_mean,
  logodds_sd = rep(0.01, length(group_id)),
  median_time_to_next_event_mean,
  median_time_to_next_event_sd = matrix(0.01, nrow = length(group_id), ncol = 3),
  visit_spacing,
  recruitment_rate = rep(1, length(group_id)), # TODO document and remove default
  max_time = 10*12, # maximal follow up time since start of trial
  logodds_min = rep(logodds(.001), length(group_id)),
  logodds_max = rep(logodds(.999), length(group_id)),
  shape_min = matrix(.99, nrow = length(group_id), ncol = 3),
  shape_max = matrix(1.01, nrow = length(group_id), ncol = 3)
) {
  mdl <- as.list(environment()) # store all input parameters
  mdl$group_id <- NULL
  mdl$visit_spacing <- NULL
  mdl <- lapply(mdl, base::as.array)
  attr(mdl, "group_id") <- as.character(group_id) # assert type
  attr(mdl, "states") <- c("stable", "response", "progression")
  attr(mdl, "visit_spacing") <- as.array(visit_spacing)
  attr(mdl, "recruitment_rate") <- as.array(recruitment_rate)
  attr(mdl, "max_time") <- as.array(max_time)
  attr(mdl, "stanmodel") <- stanmodels[["srp_model"]]
  attr(mdl, "parameter_names") <- c("p", "shape", "scale",
                                    "median_time_to_next_event")
  class(mdl) <- c("srp_model", "Model", class(mdl))
  is_valid(mdl)
  return(mdl)
}


#' @param x SRP model to format
#' @template param-dotdotdot
#' @rdname srp_model
#' @export
format.srp_model <- function(x, ...) {
  sprintf("srp_model<%s>", paste(attr(x, "group_id"), collapse = ","))
}


is_valid.srp_model <- function(mdl) { # nolint
  with(mdl, {
    checkmate::assert_vector(logodds_mean, len = length(attr(mdl, "group_id")),
                             any.missing = FALSE, .var.name = "logodds_mean")
    checkmate::assert_vector(logodds_sd, len = length(attr(mdl, "group_id")),
                             any.missing = FALSE, .var.name = "logodds_mean")
    checkmate::assert_array(median_time_to_next_event_mean, d = 2,
                             any.missing = FALSE,
                            .var.name = "median_time_to_next_event_mean")
    checkmate::assert_array(median_time_to_next_event_sd, d = 2,
                             any.missing = FALSE,
                            .var.name = "median_time_to_next_event_sd")
    checkmate::assert_vector(median_time_to_next_event_sd,
                             len = length(attr(mdl, "group_id")) *
                               length(attr(mdl, "states")), any.missing = FALSE,
                             .var.name = "median_time_to_next_event_sd")
    checkmate::assert_vector(median_time_to_next_event_mean,
                             len = length(attr(mdl, "group_id")) *
                               length(attr(mdl, "states")), any.missing = FALSE,
                             .var.name = "median_time_to_next_event_mean")
    checkmate::assertTRUE(all(logodds_mean < logodds_max),
                          .var.name = "logodds_mean < logodds_max")
    with(mdl,
      checkmate::assert_numeric(
        median_time_to_next_event_mean[median_time_to_next_event_mean < 0],
        lower = 0, upper = 0,
        .var.name = "median_time_to_next_event_mean"
      )
    )

  })
  return(TRUE)
}

# see Model.R
.impute.srp_model <- function(model, data, nsim, parameter_sample = NULL, # nolint
                              seed = NULL, p = NULL, shape = NULL,
                              scale = NULL, as_mstate = FALSE, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_groups <- length(attr(model, "group_id"))
  # make sure that either parameter sample or p, scale, shape are given
  if (!is.null(parameter_sample)) {
    stopifnot(isa(parameter_sample, "stanfit"))
    n_params_sample <- parameter_sample@sim$iter - parameter_sample@sim$warmup
  } else {
    # p, scale, and shape must be given
    if (is.null(p) || is.null(shape) || is.null(scale)) {
      stop("if no parameter sample is given all of p, scale, shape must be given") # nolint
    }
    if (length(p) != n_groups) stop()
    if (all(dim(shape) != c(n_groups, 3))) stop()
    if (all(dim(scale) != c(n_groups, 3))) stop()
  }
  # extract subject and group id levels for conversion to and back from integer
  subject_id_levels <- unique(as.character(data$subject_id))
  group_id_levels <- attr(model, "group_id") # important to maintain ordering
  if (is.null(p)) {
    # extract parameter arrays from stanfit object
    # p[i,j] is probability for ith sample for jth group
    p <- rstan::extract(parameter_sample, "p")[[1]]
  } else {
    # expand fixed parameters to same format as rstan parameters
    p <- t(array(p, dim = c(n_groups, n_params_sample)))
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
  if (is.null(shape)) {
    shape <- rstan::extract(parameter_sample, "shape")[[1]]
  } else {
    # need to add one dimension (iterations to fit format)
    shape <- aperm(
      array(shape, dim = c(n_groups, 3, n_params_sample)),
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
    response_probabilities <- p[idx[iter], , drop = FALSE]
    shapes <- matrix(shape[idx[iter], , ], ncol = 3)
    scales <- matrix(scale[idx[iter], , ], ncol = 3)
    # sample
    res <- f(data, response_probabilities, shapes, scales,
        visit_spacing = attr(model, "visit_spacing"),
        max_time = attr(model, "max_time")
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



# helper to create empty standata for model
.nodata.srp_model <- function(model) { # nolint
  tibble(
    subject_id = integer(),
    group_id = integer(),
    from = integer(),
    to = integer(),
    t_min = numeric(),
    t_max = numeric(),
    t_sot = numeric()
  )
}



# helper to create all-missing standata for model
.emptydata.srp_model <- function(model, n_per_group, seed = NULL) { # nolint
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n <- sum(n_per_group)
  group_ids <- attr(model, "group_id")
  rr <- attr(model, "recruitment_rate")
  res <- tibble()
  for (i in seq_along(group_ids)) {
    if (n_per_group[i] < 1) {
      next
    }
    subject_ids <- get_identifier(n = n_per_group[i])
    recruitment_times <- cumsum(rexp(n_per_group[i], rate = rr[i]))
    res <- bind_rows(res, tibble(
      subject_id = subject_ids,
      group_id = group_ids[i],
      t = recruitment_times,
      state = "stable" # first visits are always stable
    ))
  }
  return(res)
}



# convert time to event data to stan data list
data2standata.srp_model <- function(data, model) { # nolint
  lst_stan_data <- data %>%
    mutate(
      group_id = as.integer(factor(.data$group_id,
                                   levels = attr(model, "group_id"))),
      subject_id = as.integer(factor(.data$subject_id)),
      from = as.integer(factor(.data$from, levels = attr(model, "states"))),
      to = if_else(
          is.na(.data$to),
          "unknown",
          as.character(.data$to)
        ) %>%
        factor(levels = c(attr(model, "states"), "unknown")) %>%
        as.integer(),
      t_min = .data$t_min - .data$t_sot,
      t_max = .data$t_max - .data$t_sot
    ) %>%
    arrange(.data$subject_id, .data$from) %>%
    as.list()

  # make sure everything is an array
  for (i in seq_along(lst_stan_data)) {
    lst_stan_data[[i]] <- as.array(lst_stan_data[[i]])
  }

  lst_stan_data$M_groups <- length(attr(model, "group_id"))
  lst_stan_data$N <- nrow(data)
  lst_stan_data$N_subjects <- length(unique(data$subject_id))

  lst_stan_data <- c(lst_stan_data, as.list(model))

  return(lst_stan_data)
}

#' @importFrom stringr str_extract
#' @export
parameter_sample_to_tibble.srp_model <- function(model, sample, ...) { # nolint
  stopifnot(isa(sample, "stanfit"))
  as.matrix(sample) %>%
    as_tibble() %>%
    mutate(
      iter = row_number()
    ) %>%
    tidyr::pivot_longer(-"iter") %>%
    filter(.data$name != "lp__") %>%
    tidyr::separate("name", into = c("parameter", "group_id"),
                    sep = "\\[", fill = "right") %>%
    tidyr::separate("group_id", into = c("group_id", "transition"),
                    sep = "[\\]|,]", fill = "right", extra = "drop") %>%
    mutate(
      group_id = attr(model, "group_id")[as.integer(stringr::str_extract(
        .data$group_id, "[0-9]+"))],
      transition = as.integer(.data$transition)
    )
}



#' @inheritParams plot_mstate
#' @name srp_model
#' @export
plot_mstate.srp_model <- function(data, model, now = max(tbl_mstate$t_max), # nolint
                                  relative_to_sot = TRUE, ...) {

  starting_state <- attr(model, "states")[1]

  tbl_mstate <- data %>%
    rename(`Group ID` = "group_id")

  if (relative_to_sot) {
    tbl_mstate <- tbl_mstate %>%
      mutate(
        t_min = .data$t_min - .data$t_sot,
        t_max = .data$t_max - .data$t_sot,
        t_sot = 0
      )
  }

  tbl_points <- tbl_mstate %>%
    mutate(
      tmp = purrr::pmap(
        list(.data$from, .data$to, .data$t_min, .data$t_max, .data$t_sot),
        ~if (is.na(..2)) {
          tibble(t = c(..3, ..5), state = c(..1, starting_state))
        } else {
          tibble(t = c(..3, ..4, ..5), state = c(..1, ..2, starting_state))
        }
      )
    ) %>%
    select(all_of(c("subject_id", "Group ID", "tmp"))) %>%
    tidyr::unnest("tmp") %>%
    filter(is.finite(.data$t), .data$t < now) %>%
    distinct() %>%
    arrange(.data$subject_id, .data$t)

  tbl_intervals <- tbl_mstate %>%
    bind_rows(
      select(tbl_mstate, all_of(c("subject_id", "Group ID", "t_sot"))) %>%
        distinct() %>%
        mutate(
          from = starting_state,
          to = starting_state,
          t_min = .data$t_sot,
          t_max = .data$t_sot
        )
    ) %>%
    arrange(.data$subject_id, .data$t_min, .data$t_max) %>%
    distinct() %>%
    group_by(.data$subject_id) %>%
    transmute(
      .data$subject_id,
      .data$`Group ID`,
      state = if_else(.data$to == lead(.data$from), lead(.data$from),
                      NA_character_),
      tmp1 = .data$t_max,
      tmp2 = lead(.data$t_min)
    ) %>%
    ungroup() %>%
    filter(!is.na(.data$state), is.finite(.data$tmp1), is.finite(.data$tmp2),
           .data$tmp2 > .data$tmp1)

  tbl_at_risk <- tbl_mstate %>%
    filter(.data$t_max == Inf) %>%
    transmute(
      .data$subject_id,
      .data$`Group ID`,
      t = .data$t_min,
      state = .data$from
    )

  tbl_censored <- tbl_mstate %>%
    filter(.data$t_max == -Inf) %>%
    transmute(
      .data$subject_id,
      .data$`Group ID`,
      t = .data$t_min,
      state = .data$from
    )

  scale <- max(tbl_points$t)

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$tmp1, xend = .data$tmp2,
                   y = .data$subject_id, yend = .data$subject_id,
                   color = .data$state),
      data = tbl_intervals
    ) +
    ggplot2::geom_point(ggplot2::aes(.data$t, .data$subject_id,
                                     color = .data$state), data = tbl_points) +
    ggplot2::geom_segment(
      ggplot2::aes(.data$t, .data$subject_id, xend = .data$t + scale / 33,
                   yend = .data$subject_id, color = .data$state),
      arrow = ggplot2::arrow(type = "closed", angle = 10,
                             length = ggplot2::unit(0.05, "npc")),
      data = tbl_at_risk
    ) +
    ggplot2::geom_point(ggplot2::aes(.data$t, .data$subject_id,
                                     color = .data$state), shape = "x",
                        size = 5, data = tbl_censored) +
    ggplot2::geom_vline(xintercept = now) +
    ggplot2::labs(x = if (relative_to_sot) "Time since SoT" else
      "Time since first SoT", y = "Subject ID") +
    ggplot2::scale_color_discrete("") +
    ggplot2::facet_wrap(~.data$`Group ID`, ncol = 1,
                        labeller = ggplot2::label_both,
                        strip.position = "right", scales = "free_y") +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 5),
      legend.position = "right"
    )

}



#' @inheritParams visits_to_mstate
#' @name srp_model
#' @export
visits_to_mstate.srp_model <- function(tbl_visits, model,
                                       now = max(tbl_visits$t),
                                       eof_indicator = "EOF") {
  # make sure everything is sorted
  tbl_visits <- arrange(tbl_visits, .data$subject_id, .data$t)

  tbl_mstate <- list()

  subject_id_lagged <- 0L
  state_lagged <- 0L
  t_sot <- 0 # start of treatment
  for (i in seq_len(nrow(tbl_visits))) {
    if (tbl_visits$subject_id[i] != subject_id_lagged || i == 1) {
      # switch to new subject
      subject_id_lagged <- tbl_visits$subject_id[i]
      state_lagged <- "stable"
      t_sot <- tbl_visits$t[i]
      if (tbl_visits$state[i] != state_lagged) {
        # record jump
        stop(sprintf(
          "first visit must be in starting state; subject_id=%s, state=%s",
          tbl_visits$subject_id[i],
          tbl_visits$state[i]
        ))
      }
    }
    # handle jumps
    if (tbl_visits$state[i] != state_lagged) {
      if (tbl_visits$state[i] == eof_indicator) { # record eof
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
          subject_id = tbl_visits$subject_id[i],
          group_id = tbl_visits$group_id[i],
          from = state_lagged,
          to = NA,
          t_min = tbl_visits$t[i],
          # - Inf indicates censoring and end of follow up
          # (event can no longer be observed)
          t_max = -Inf,
          t_sot = t_sot
        ))
        state_lagged <- tbl_visits$state[i]
      } else { # record jump
        tbl_mstate <- bind_rows(tbl_mstate, tibble(
          subject_id = tbl_visits$subject_id[i],
          group_id = tbl_visits$group_id[i],
          from = state_lagged,
          to = tbl_visits$state[i],
          t_min = tbl_visits$t[i - 1],
          t_max = tbl_visits$t[i],
          t_sot = t_sot
        ))
      }
      # updated current state
      state_lagged <- tbl_visits$state[i]
    }
    # handle non-eof censoring
    censored <- FALSE
    if (!(tbl_visits$state[i] %in% c("progression", eof_indicator))) {
      if (i < nrow(tbl_visits)) {
        if (tbl_visits$subject_id[i] != tbl_visits$subject_id[i + 1]) {
          censored <- TRUE
        }
      } else {
        censored <- TRUE
      }
    }
    if (censored) {
      tbl_mstate <- bind_rows(tbl_mstate, tibble(
        subject_id = tbl_visits$subject_id[i],
        group_id = tbl_visits$group_id[i],
        from = tbl_visits$state[i],
        to = NA,
        t_min = now,
        t_max = Inf,
        t_sot = t_sot
      ))
    }
  }
  return(tbl_mstate)
}



#' @export
plot.srp_model <- function(x, dt, sample = NULL, seed = NULL,
                           n_grid = 50, ...) {
  if (is.null(sample)) {
    sample <- sample_prior(x, seed = seed, ...)
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("the patchwork package is required to plot SRP models")
  }
  tbl_sample <- parameter_sample_to_tibble(x, sample)
  # plot transition times
  p1 <- tbl_sample %>%
    filter(.data$parameter %in% c("shape", "scale")) %>%
    tidyr::pivot_wider(names_from = "parameter",
                       values_from = "value") %>%
    tidyr::expand_grid(dt = seq(dt[1], dt[2], length.out = n_grid)) %>%
    mutate(
      survival = 1 - stats::pweibull(.data$dt,
                                     shape = .data$shape, scale = .data$scale)
    ) %>%
    group_by(.data$group_id, .data$transition, .data$dt) %>%
    summarize(
      survival = mean(.data$survival), .groups = "drop"
    ) %>%
    filter(
      is.finite(.data$survival)
    ) %>%
    mutate(
      transition = case_when(
        .data$transition == 1 ~ "stable to response",
        .data$transition == 2 ~ "stable to progression",
        .data$transition == 3 ~ "response to progression",
      ) %>%
      factor(levels = c("stable to response", "stable to progression",
                        "response to progression"))
    ) %>%
    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(.data$dt, .data$survival,
                                      color = .data$group_id)) +
      ggplot2::labs(x = "time to next event", y = "'Survival' fraction") +
      ggplot2::scale_color_discrete("") +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, by = .1)) +
      ggplot2::facet_wrap(~.data$transition, nrow = 1) +
      ggplot2::theme(
        legend.position = "top",
        panel.grid.minor = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(1.5, "lines")
      )
  # plot response probability
  p2 <- tbl_sample %>%
    filter(.data$parameter == "p") %>%
    ggplot2::ggplot() +
      ggplot2::stat_ecdf(ggplot2::aes(.data$value,
                                      color = .data$group_id), geom = "line") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(x = "response probability", y = "CDF") +
      ggplot2::scale_x_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, by = .1),
                                  expand = c(0, 0)) +
      ggplot2::scale_y_continuous(limits = c(0, 1),
                                  breaks = seq(0, 1, by = .1)) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      )
  # plot pfs
  tbl_pfs_survival <- sample_pfs_rate(
      x,
      # ToDo: make sure this works from 0
      t = seq(0.5, dt[2], length.out = n_grid),
      sample = sample
    ) %>%
    # integrate over prior sample
    group_by(.data$group_id, .data$t) %>%
    summarize(pfs = mean(.data$pfs), .groups = "drop")
  p3 <- ggplot2::ggplot(tbl_pfs_survival) +
    ggplot2::geom_line(ggplot2::aes(x = .data$t,
                                    y = .data$pfs, color = .data$group_id)) +
    ggplot2::labs(x = "time", y = "PFS") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = .1)) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  design <- "
  111
  234
  "
  p1 + p2 + p3 + patchwork::guide_area() +
    patchwork::plot_layout(design = design, guides = "collect")
}



#' @export
sample_pfs_rate.srp_model <- function( # nolint
  model,
  t, # PFS_r is 1 - Pr[progression or death before time t]
  sample = NULL,
  warmup = 500L,
  nsim = 2000L,
  seed = NULL,
  ...
) {
  if (is.null(sample)) {
    sample <- sample_prior(model, warmup = warmup, nsim = nsim, seed = seed,
      rstan_output = TRUE, pars = attr(model, "parameter_names"), ...
    )
  }
  sample <- parameter_sample_to_tibble(model, sample)
  pr_direct_progression <- function(shape_2, scale_2, t) {
    return(stats::pweibull(t, shape_2, scale_2))
  }
  pr_indirect_progression <- function(shape_1, shape_3, scale_1, scale_3, t) {
    # need to integrate over response time
    integrand <- function(t_response) {
      stats::dweibull(t_response, shape_1, scale_1) *
        stats::pweibull(t - t_response, shape_3, scale_3)
    }
    # can reduce absolute tolerance substantially ok for probabilities
    res <- stats::integrate(
      integrand, lower = 0, upper = t, rel.tol = 1e-5, abs.tol = 1e-5
    )
    return(res$value)
  }
  tbl_pfs <- sample %>%
    # pivot parameters
    tidyr::pivot_wider(
      names_from = all_of(c("parameter", "transition")),
      values_from = "value"
    ) %>%
    rename(pr_response = "p_NA") %>%
    # cross with time points
    tidyr::expand_grid(t = t) %>%
    # compute PFS before t
    mutate(
      pfs = purrr::pmap_dbl(
        list(
          .data$pr_response,
          .data$scale_1, .data$scale_2, .data$scale_3,
          .data$shape_1, .data$shape_2, .data$shape_3,
          t
        ),
        function(pr_response, scale_1, scale_2, scale_3, shape_1,
                 shape_2, shape_3, t) {
          pr_progression_t <- pr_response *
            pr_indirect_progression(shape_1, shape_3, scale_1, scale_3, t) +
            (1 - pr_response) * pr_direct_progression(shape_2, scale_2, t)
          return(1 - pr_progression_t)
        }
      )
    ) %>%
    # only keep interesting stuff
    select(all_of(c("iter", "group_id", "t", "pfs")))
  return(tbl_pfs)
}
