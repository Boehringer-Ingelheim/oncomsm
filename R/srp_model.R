#' todo
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
  attr(mdl, "parameter_names") <- c("p", "shape", "scale", "median_time_to_next_event")
  class(mdl) <- c("srp_model", "Model", class(mdl))
  return(mdl)
}




# see Model.R
.impute.srp_model <- function(model, data, nsim, parameter_sample, seed = NULL, ...) {
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
    arrange(t_sot, subject_id, (t_min + t_max)/2) %>%
    mutate(
      subject_id = as.integer(factor(as.character(subject_id), levels = subject_id_levels)),
      group_id = as.integer(factor(group_id, levels = group_id_levels))
    )
  data %>%
    group_by(subject_id) %>%
    filter(from != lag(to)) %>%
    nrow() %>%
    {assertthat::assert_that(. == 0)}

  res <- tribble(
      ~subject_id, ~group_id, ~from, ~to, ~t_min, ~t_max, ~t_sot, ~iter
    )

  visit_spacing <- attr(model, "visit_spacing")
  idx <- sample(1:dim(p)[1], size = nsim, replace = TRUE)

  for (i in 1:nrow(data)) {
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
          pr_survival_response <- 1 - stats::pweibull(data$t_min[i] - data$t_sot[i], sshape[1], sscale[1])
          pr_survival_progression <- 1 - stats::pweibull(data$t_min[i] - data$t_sot[i], sshape[2], sscale[2])
          pr_response <- pr_response_raw * pr_survival_response / (
            pr_response_raw * pr_survival_response +
              (1 - pr_response_raw) * pr_survival_progression
          )
          response <- rbinom(1, 1, pr_response)
          if (response) {
            # sample exact response time
            t_response <- bhmbasket.predict:::rtruncweibull(
              sshape[1], scale = sscale[1], data$t_min[i], Inf # t_min since time of SoT is known
            )
            # apply visit scheme
            n_visits_response <- t_response %/% visit_spacing[g]
            tmin_response <- data$t_min[i] + visit_spacing[g] * n_visits_response
            tmax_response <- data$t_min[i] + visit_spacing[g] * (n_visits_response + 1)
            # sample subsequent progression,
            dt_progression <- stats::rweibull(1, sshape[3], sscale[3])
            # apply visit scheme
            n_visits_progression <- (dt_progression + t_response) %/% visit_spacing[g]
            tmin_progression <- data$t_min[i] + visit_spacing[g] * n_visits_progression
            tmax_progression <- data$t_min[i] + visit_spacing[g] * (n_visits_progression + 1)
            res <- bind_rows(
              res, tribble(
                       ~subject_id, ~group_id,      ~from,           ~to,           ~t_min,           ~t_max,        ~t_sot, ~iter,
                data$subject_id[i],         g,   "stable",    "response",    tmin_response,    tmax_response, data$t_sot[i],      j,
                data$subject_id[i],         g, "response", "progression", tmin_progression, tmax_progression, data$t_sot[i],      j
              )
            )
          } else { # sample progression directly
            dt_progression <- bhmbasket.predict:::rtruncweibull(
              sshape[2], scale = sscale[2], data$t_min[i], Inf # t_min since time of SoT is known
            )
            # apply visit scheme
            n_visits_progression <- dt_progression %/% visit_spacing[g]
            tmin_progression <- data$t_min[i] + visit_spacing[g] * n_visits_progression
            tmax_progression <- data$t_min[i] + visit_spacing[g] * (n_visits_progression + 1)
            res <- bind_rows(
              res, tribble(
                       ~subject_id, ~group_id,    ~from,           ~to,           ~t_min,           ~t_max,        ~t_sot, ~iter,
                data$subject_id[i],         g, "stable", "progression", tmin_progression, tmax_progression, data$t_sot[i], j
              )
            )
          } # end stable -> progression
        } # end from == 1
        if (data$from[i] == "response") {
          if (data$from[i - 1] != "stable" || data$subject_id[i - 1] != data$subject_id[i]) {
            stop()
          }
          # sample exact response time
          t_response <- bhmbasket.predict:::rtruncweibull(
            sshape[1], scale = sscale[1], data$t_min[i - 1], data$t_max[i - 1]
          )
          # sample progression time
          dt_progression <- bhmbasket.predict:::rtruncweibull(
            sshape[3], scale = sscale[3], data$t_min[i] - t_response, Inf
          )
          # apply visit scheme
          n_visits_progression <- dt_progression %/% visit_spacing[g]
          tmin_progression <- data$t_min[i] + visit_spacing[g] * n_visits_progression
          tmax_progression <- data$t_min[i] + visit_spacing[g] * (n_visits_progression + 1)
          res <- bind_rows(
            res, tribble(
                     ~subject_id, ~group_id,      ~from,           ~to,           ~t_min,           ~t_max,        ~t_sot, ~iter,
              data$subject_id[i],         g, "response", "progression", tmin_progression, tmax_progression, data$t_sot[i],     j
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
          factor(subject_id, levels = seq_along(subject_id_levels), labels = subject_id_levels)
        ),
      group_id = as.character(
          factor(group_id, levels = seq_along(group_id_levels), labels = group_id_levels)
        )
    )
  return(res)
}



# helper to create empty standata for model
.nodata.srp_model <- function(model) {
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
.emptydata.srp_model <- function(model, n_per_arm) {
  n <- sum(n_per_arm)
  tibble(
    subject_id = 1:n,
    group_id = attr(model, "group_id")[rep(1:length(n_per_arm), times = n_per_arm)],
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
data2standata.srp_model <- function(model, data) {
  lst_stan_data <- data %>%
    mutate(
      group_id = as.integer(factor(.data$group_id, levels = attr(model, "group_id"))),
      subject_id = as.integer(factor(.data$subject_id)),
      from = as.integer(factor(.data$from, levels = attr(model, "states"))),
      to = .data$to %>%
        {if_else(is.na(.), "unknown", as.character(.))} %>%
        factor(levels = c(attr(model, "states"), "unknown")) %>%
        as.integer(),
      t_min = t_min - t_sot,
      t_max = t_max - t_sot
    ) %>%
    arrange(subject_id, from) %>% # make sure verything is sorted properly
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
parameter_sample_to_tibble.srp_model <- function(model, sample) {
  stopifnot(isa(sample, "stanfit"))
  as.matrix(sample) %>%
    as_tibble() %>%
    mutate(
      iter = row_number()
    ) %>%
    tidyr::pivot_longer(-.data$iter) %>%
    filter(.data$name != "lp__") %>%
    tidyr::separate(.data$name, into = c("parameter", "group_id"), sep = "\\[", fill = "right") %>%
    tidyr::separate(.data$group_id, into = c("group_id", "transition"), sep = "[\\]|,]", fill = "right", extra = "drop") %>%
    mutate(
      group_id = attr(model, "group_id")[as.integer(stringr::str_extract(.data$group_id, "[0-9]+"))],
      transition = as.integer(transition)
    )
}



#' @export
plot_mstate.srp_model <- function(model, tbl_mstate, now = max(tbl_mstate$t_max), relative_to_sot = TRUE, ...) {

  starting_state <- attr(mdl, "states")[1]

  tbl_mstate <- tbl_mstate %>%
    rename(`Group ID` = group_id)

  if (relative_to_sot) {
    tbl_mstate <- tbl_mstate %>%
      mutate(
        t_min = t_min - t_sot,
        t_max = t_max - t_sot,
        t_sot = 0
      )
  }

  tbl_points <- tbl_mstate %>%
    # filter(t_max != -Inf) %>%
    mutate(
      tmp = purrr::pmap(
        list(from, to, t_min, t_max, t_sot),
        ~if (is.na(..2)) {
          tibble(t = c(..3, ..5), state = c(..1, starting_state))
        } else {
          tibble(t = c(..3, ..4, ..5), state = c(..1, ..2, starting_state))
        }
      )
    ) %>%
    select(subject_id, `Group ID`, tmp) %>%
    tidyr::unnest(tmp) %>%
    filter(is.finite(t), t < now) %>%
    distinct() %>%
    arrange(subject_id, t)

  tbl_intervals <- tbl_mstate %>%
    bind_rows(
      select(tbl_mstate, subject_id, `Group ID`, t_sot) %>%
        distinct() %>%
        mutate(from = starting_state, to = starting_state, t_min = t_sot, t_max = t_sot)
    ) %>%
    arrange(subject_id, t_min, t_max) %>%
    distinct() %>%
    group_by(subject_id) %>%
    transmute(
      subject_id,
      `Group ID`,
      state = if_else(to == lead(from), lead(from), NA_character_),
      tmp1 = t_max,
      tmp2 = lead(t_min)
    ) %>%
    ungroup() %>%
    filter(!is.na(state), is.finite(tmp1), is.finite(tmp2), tmp2 > tmp1)

  tbl_at_risk <- tbl_mstate %>%
    filter(t_max == Inf) %>%
    transmute(
      subject_id,
      `Group ID`,
      t = t_min,
      state = from
    )

  tbl_censored <- tbl_mstate %>%
    filter(t_max == -Inf) %>%
    transmute(
      subject_id,
      `Group ID`,
      t = t_min,
      state = from
    )

  scale <- max(tbl_points$t)

  ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = tmp1, xend = tmp2, y = subject_id, yend = subject_id, color = state), data = tbl_intervals) +
    ggplot2::geom_point(ggplot2::aes(t, subject_id, color = state), data = tbl_points) +
    ggplot2::geom_segment(
      ggplot2::aes(t, subject_id, xend = t + scale/33, yend = subject_id, color = state),
      arrow = ggplot2::arrow(type = "closed", angle = 10, length = ggplot2::unit(0.05, "npc")),
      data = tbl_at_risk
    ) +
    ggplot2::geom_point(ggplot2::aes(t, subject_id, color = state), shape = "x", size = 5, data = tbl_censored) +
    ggplot2::geom_vline(xintercept = now) +
    ggplot2::labs(x = if (relative_to_sot) "Time since SoT" else "Time since first SoT", y = "Subject ID") +
    ggplot2::scale_color_discrete("") +
    ggplot2::facet_wrap(~`Group ID`, ncol = 1, labeller = ggplot2::label_both, strip.position = "right", scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 5),
      legend.position = "right"
    )

}



#' @export
plot.srp_model <- function(model, dt, sample = NULL, seed = NULL, n_grid = 50, ...) {
  if(is.null(sample)) {
    sample <- sample_prior(model, seed = seed, ...)
  }
  sample <- parameter_sample_to_tibble(model, sample)
  # plot transition times
  p1 <- sample %>%
    filter(parameter %in% c("shape", "scale")) %>%
    tidyr::pivot_wider(names_from = parameter, values_from = value) %>%
    tidyr::expand_grid(dt = seq(dt[1], dt[2], length.out = n_grid)) %>%
    mutate(
      pdf = pweibull(dt, shape = shape, scale = scale)
    ) %>%
    group_by(group_id, transition, dt) %>%
    summarize(pdf = mean(pdf), .groups = "drop") %>%
    filter(is.finite(pdf)) %>%
    mutate(
      transition = case_when(
        transition == 1 ~ "stable to response",
        transition == 2 ~ "stable to progression",
        transition == 3 ~ "response to progression",
      ) %>%
      factor(levels = c("stable to response", "stable to progression", "response to progression"))
    ) %>%
    ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(dt, pdf, color = group_id)) +
      ggplot2::labs(x = "time to next event", y = "CDF") +
      ggplot2::scale_color_discrete("") +
      ggplot2::facet_wrap(~transition, nrow = 1) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "top",
        panel.grid.minor = ggplot2::element_blank()
      )
  # plot response probability
  p2 <- sample %>%
    filter(parameter == "p") %>%
    ggplot2::ggplot() +
      ggplot2::stat_ecdf(aes(value, color = group_id), geom = "line") +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      ggplot2::labs(x = "response probability", y = "ECDF") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      )
  design <- "
  111
  223
  "
  p1 + p2 + patchwork::guide_area() + patchwork::plot_layout(design = design, guides = "collect")
}



#' @export
sample_pfs_rate.srp_model <- function(
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
  pr_direct_progression <- function(shape_3, scale_3, t) {
    return(pweibull(t, shape_3, scale_3))
  }
  pr_indirect_progression <- function(shape_1, shape_2, scale_1, scale_2, t) {
    # need to integrate over response time
    integrand <- function(t_response) {
      dweibull(t_response, shape_1, scale_1) * pweibull(t - t_response, shape_2, scale_2)
    }
    # can reduce absolute tolerance substantially - doesn't matter for probabilities
    res <- integrate(
      integrand, lower = 0, upper = t, rel.tol = 1e-5, abs.tol = 1e-5
    )
    return(res$value)
  }
  tbl_pfs <- sample %>%
    # pivot parameters
    tidyr::pivot_wider(names_from = c(parameter, transition), values_from = value) %>%
    rename(pr_response = p_NA) %>%
    # cross with time points
    tidyr::expand_grid(t = t) %>%
    # compute PFS before t
    mutate(
      pfs = purrr::pmap_dbl(
        list(pr_response, scale_1, scale_2, scale_3, shape_1, shape_2, shape_3, t),
        function(pr_response, scale_1, scale_2, scale_3, shape_1, shape_2, shape_3, t) {
          pr_progression_t <- pr_response * pr_indirect_progression(shape_1, shape_2, scale_1, scale_2, t) +
            (1 - pr_response) * pr_direct_progression(shape_3, scale_3, t)
          return(1 - pr_progression_t)
        }
      )
    ) %>%
    select(iter, group_id, t, pfs) # only keep interesting stuff
  return(tbl_pfs)
}
