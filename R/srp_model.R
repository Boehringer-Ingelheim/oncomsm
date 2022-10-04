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
  logodds_sd,
  median_time_to_next_event_mean,
  median_time_to_next_event_sd,
  visit_spacing,
  logodds_min = rep(logodds(.001), length(group_id)),
  logodds_max = rep(logodds(.999), length(group_id)),
  shape_min = matrix(.99, nrow = length(group_id), ncol = 3),
  shape_max = matrix(1.01, nrow = length(group_id), ncol = 3)
) {
  mdl <- as.list(environment()) # store all input parameters
  mdl$group_id <- NULL
  mdl$visit_spacing <- NULL
  mdl <- lapply(mdl, base::as.array)
  attr(mdl, "group_id") <- group_id
  attr(mdl, "states") <- c("stable", "response", "progression")
  attr(mdl, "visit_spacing") <- as.array(visit_spacing)
  attr(mdl, "stanmodel") <- stanmodels[["srp_model"]]
  attr(mdl, "parameter_names") <- c("p", "shape", "scale",
                                    "median_time_to_next_event")
  class(mdl) <- c("srp_model", "Model", class(mdl))
  is_valid(mdl)
  return(mdl)
}


#' To check the validity of input parameters to create_srp
#'
#' @param mdl an object of type srp_model
#'
#' @return True statements or corresponding parameter specific errors
#' @export
#'
#' @importFrom checkmate assert_vector
is_valid.srp_model <- function(mdl) { # nolint
  with(mdl, {
    stopifnot(logodds_mean[1] < logodds_max[1])
    stopifnot(logodds_mean[2] < logodds_max[2])
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
.impute.srp_model <- function(model, data, nsim, parameter_sample, # nolint
                              seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # TODO: convert data to matrix and process in c++
  stopifnot(isa(parameter_sample, "stanfit"))
  # extract subject and group id levels for conversion to and back from integer
  subject_id_levels <- unique(as.character(data$subject_id))
  group_id_levels <- attr(model, "group_id") # important to maintain ordering
  # extract parameter matrices
  p <- rstan::extract(parameter_sample, "p")[[1]]
  scale <- rstan::extract(parameter_sample, "scale")[[1]]
  shape <- rstan::extract(parameter_sample, "shape")[[1]]

  data <- data %>%
    arrange(.data$t_sot, .data$subject_id, (.data$t_min + .data$t_max) / 2) %>%
    mutate(
      subject_id = as.integer(factor(as.character(.data$subject_id),
                                     levels = subject_id_levels)),
      group_id = as.integer(factor(.data$group_id, levels = group_id_levels))
    )

  res <- tribble(
      ~subject_id, ~group_id, ~from, ~to, ~t_min, ~t_max, ~t_sot, ~iter
    )

  visit_spacing <- attr(model, "visit_spacing")
  idx <- sample(seq_len(dim(p)[1]), size = nsim, replace = TRUE)

  for (i in seq_len(nrow(data))) {
  if (!is.na(data$to[i])) { # observed transition, nothing to sample
    res <- bind_rows(res, tidyr::expand_grid(data[i, ], iter = 1:nsim))
  } else { # a censored transition
    for (j in 1:nsim) {
      g <- data$group_id[i]
      k <- idx[j]
      sshape <- shape[k, g, ]
      sscale <- scale[k, g, ]
      if (data$from[i] == "stable") {
        # first sample response/progression
        pr_response_raw <- p[k, g] # use survival information!
        pr_survival_response <- 1 - stats::pweibull(
          data$t_min[i] -
            data$t_sot[i],
          sshape[1], sscale[1]
        )
        pr_survival_progression <- 1 - stats::pweibull(
          data$t_min[i] -
            data$t_sot[i],
          sshape[2], sscale[2]
        )
        pr_response <- pr_response_raw * pr_survival_response / (
          pr_response_raw * pr_survival_response +
            (1 - pr_response_raw) * pr_survival_progression
        )
        response <- stats::rbinom(1, 1, pr_response)
        if (response) {
          # sample exact response time
          t_response <- rtruncweibull(
            sshape[1],
            # t_min since time of SoT is known
            scale = sscale[1], data$t_min[i], Inf
          )
          # apply visit scheme
          n_visits_response <- t_response %/% visit_spacing[g]
          tmin_response <- data$t_min[i] +
            visit_spacing[g] * n_visits_response
          tmax_response <- data$t_min[i] +
            visit_spacing[g] * (n_visits_response + 1)
          # sample subsequent progression,
          dt_progression <- stats::rweibull(1, sshape[3], sscale[3])
          # apply visit scheme
          n_visits_progression <- (dt_progression +
            t_response) %/% visit_spacing[g]
          tmin_progression <- data$t_min[i] +
            visit_spacing[g] * n_visits_progression
          tmax_progression <- data$t_min[i] +
            visit_spacing[g] * (n_visits_progression + 1)
          res <- bind_rows(
            res, tribble(
              ~subject_id, ~group_id, ~from, ~to, ~t_min, ~t_max, ~t_sot,
                ~iter,
              data$subject_id[i], g, "stable", "response", tmin_response,
                tmax_response, data$t_sot[i], j,
              data$subject_id[i], g, "response", "progression",
                tmin_progression, tmax_progression, data$t_sot[i], j
            )
          )
        } else { # sample progression directly
          dt_progression <- rtruncweibull(
            sshape[2],
            scale = sscale[2],
            data$t_min[i], Inf # t_min since time of SoT is known
          )
          # apply visit scheme
          n_visits_progression <- dt_progression %/% visit_spacing[g]
          tmin_progression <- data$t_min[i] +
            visit_spacing[g] * n_visits_progression
          tmax_progression <- data$t_min[i] +
            visit_spacing[g] * (n_visits_progression + 1)
          res <- bind_rows(
            res, tribble(
              ~subject_id, ~group_id, ~from, ~to, ~t_min, ~t_max, ~t_sot, ~iter,
              data$subject_id[i], g, "stable", "progression", tmin_progression,
                tmax_progression, data$t_sot[i], j
            )
          )
        } # end stable -> progression
      } # end from == 1
      if (data$from[i] == "response") {
        if (data$from[i - 1] != "stable" ||
          data$subject_id[i - 1] != data$subject_id[i]) {
          stop()
        }
        # sample exact response time
        t_response <- rtruncweibull(
          sshape[1],
          scale = sscale[1], data$t_min[i - 1], data$t_max[i - 1]
        )
        # sample progression time
        dt_progression <- rtruncweibull(
          sshape[3],
          scale = sscale[3], data$t_min[i] - t_response, Inf
        )
        # apply visit scheme
        n_visits_progression <- dt_progression %/% visit_spacing[g]
        tmin_progression <- data$t_min[i] +
          visit_spacing[g] * n_visits_progression
        tmax_progression <- data$t_min[i] +
          visit_spacing[g] * (n_visits_progression + 1)
        res <- bind_rows(
          res, tribble(
            ~subject_id, ~group_id, ~from, ~to, ~t_min, ~t_max, ~t_sot, ~iter,
            data$subject_id[i], g, "response", "progression", tmin_progression,
              tmax_progression, data$t_sot[i], j
          )
        )
      } # end from == 2
    } # end iterate of j
  } # end if/else
} # end iteration over i
  # convert subject and group id back
  res <- res %>%
    mutate(
      subject_id = as.character(
          factor(.data$subject_id, levels = seq_along(subject_id_levels),
                 labels = subject_id_levels)
        ),
      group_id = as.character(
          factor(.data$group_id, levels = seq_along(group_id_levels),
                 labels = group_id_levels)
        )
    )
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
.emptydata.srp_model <- function(model, n_per_group) { # nolint
  n <- sum(n_per_group)
  tibble(
    subject_id = 1:n,
    group_id = attr(model, "group_id")[rep(seq_len(length(n_per_group)),
                                           times = n_per_group)],
    from = rep("stable", n),
    to = rep(NA, n),
    # assume everyone is recruited at zero
    # only t_min/max - t_sot is relevant anyhow
    t_min = rep(0, n),
    t_max = rep(Inf, n),
    t_sot = rep(0, n)
  )
}



# convert time to event data to stan data list
data2standata.srp_model <- function(model, data) { # nolint
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
    tidyr::pivot_longer(-.data$iter) %>%
    filter(.data$name != "lp__") %>%
    tidyr::separate(.data$name, into = c("parameter", "group_id"),
                    sep = "\\[", fill = "right") %>%
    tidyr::separate(.data$group_id, into = c("group_id", "transition"),
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
plot_mstate.srp_model <- function(model, data, now = max(tbl_mstate$t_max), # nolint
                                  relative_to_sot = TRUE, ...) {

  starting_state <- attr(model, "states")[1]

  tbl_mstate <- data %>%
    rename(`Group ID` = .data$group_id)

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
    select(.data$subject_id, .data$`Group ID`, .data$tmp) %>%
    tidyr::unnest(.data$tmp) %>%
    filter(is.finite(.data$t), .data$t < now) %>%
    distinct() %>%
    arrange(.data$subject_id, .data$t)

  tbl_intervals <- tbl_mstate %>%
    bind_rows(
      select(tbl_mstate, .data$subject_id, .data$`Group ID`, .data$t_sot) %>%
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
    tidyr::pivot_wider(names_from = .data$parameter,
                       values_from = .data$value) %>%
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
      names_from = c(.data$parameter, .data$transition),
      values_from = .data$value
    ) %>%
    rename(pr_response = .data$p_NA) %>%
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
    select(.data$iter, .data$group_id, .data$t, .data$pfs)
  return(tbl_pfs)
}
