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
.impute.srp_model <- function(model, data, nsim, now, parameter_sample, seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # TODO: convert data to matrix and process in c++
  stopifnot(isa(parameter_sample, "stanfit"))
  # extract parameter matrices
  p <- rstan::extract(parameter_sample, "p")[[1]]
  scale <- rstan::extract(parameter_sample, "scale")[[1]]
  shape <- rstan::extract(parameter_sample, "shape")[[1]]

  data <- data %>%
    arrange(subject_id, (dt_min + dt_max)/2)
  data %>%
    group_by(subject_id) %>%
    filter(from != lag(to)) %>%
    nrow() %>%
    {assertthat::assert_that(. == 0)}

  res <- tribble(
      ~subject_id, ~group_id, ~from, ~to, ~dt_min, ~dt_max, ~iter
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
          pr_survival_response <- 1 - stats::pweibull(data$dt_min[i],  sshape[1], sscale[1])
          pr_survival_progression <- 1 - stats::pweibull(data$dt_min[i],  sshape[2], sscale[2])
          pr_response <- pr_response_raw * pr_survival_response / (
            pr_response_raw * pr_survival_response +
              (1 - pr_response_raw) * pr_survival_progression
          )
          response <- rbinom(1, 1, pr_response)
          if (response) {
            # sample exact response time
            dt_response <- bhmbasket.predict:::rtruncweibull(
              sshape[1], scale = sscale[1], data$dt_min[i], Inf
            )
            # apply visit scheme
            n_visits_response <- (dt_response - data$dt_min[i]) %/% visit_spacing[g]
            dtmin_response <- data$dt_min[i] + visit_spacing[g] * n_visits_response
            dtmax_response <- data$dt_min[i] + visit_spacing[g] * (n_visits_response + 1)
            # sample subsequent progression,
            dt_progression <- stats::rweibull(1, sshape[3], sscale[3])
            # apply visit scheme
            n_visits_progression <- (dt_progression + dt_response - data$dt_min[i]) %/% visit_spacing[g]
            dtmin_progression <- data$dt_min[i] + visit_spacing[g] * n_visits_progression
            dtmax_progression <- data$dt_min[i] + visit_spacing[g] * (n_visits_progression + 1)
            res <- bind_rows(
              res, tribble(
                       ~subject_id, ~group_id,      ~from,           ~to,           ~dt_min,           ~dt_max, ~iter,
                data$subject_id[i],         g,   "stable",    "response",    dtmin_response,    dtmax_response,     j,
                data$subject_id[i],         g, "response", "progression", dtmin_progression, dtmax_progression,     j
              )
            )
          } else { # sample progression/cure directly
            dt_progression <- bhmbasket.predict:::rtruncweibull(
              sshape[2], scale = sscale[2], data$dt_min[i], Inf
            )
            # apply visit scheme
            n_visits_progression <- (dt_progression - data$dt_min[i]) %/% visit_spacing[g]
            dtmin_progression <- data$dt_min[i] + visit_spacing[g] * n_visits_progression
            dtmax_progression <- data$dt_min[i] + visit_spacing[g] * (n_visits_progression + 1)
            res <- bind_rows(
              res, tribble(
                       ~subject_id, ~group_id,    ~from,           ~to,           ~dt_min,           ~dt_max, ~iter,
                data$subject_id[i],         g, "stable", "progression", dtmin_progression, dtmax_progression,     j
              )
            )
          } # end stable -> progression
        } # end from == 1
        if (data$from[i] == "resp") {
          dt_progression <- bhmbasket.predict:::rtruncweibull(
            sshape[3], scale = sscale[3], data$dt_min[i], Inf
          )
          # apply visit scheme
          n_visits_progression <- (dt_progression - data$dt_min[i]) %/% visit_spacing[g]
          dtmin_progression <- data$dt_min[i] + visit_spacing[g] * n_visits_progression
          dtmax_progression <- data$dt_min[i] + visit_spacing[g] * (n_visits_progression + 1)
          res <- bind_rows(
            res, tribble(
                     ~subject_id, ~group_id,      ~from,           ~to,           ~dt_min,           ~dt_max, ~iter,
              data$subject_id[i],         g, "response", "progression", dtmin_progression, dtmax_progression,     j
            )
          )
        } # end from == 2
      } # end iterate of j
    } # end if/else
  } # end iteration over i
  return(res)
}



# helper to create empty standata for model
.nodata.srp_model <- function(model) {
  tibble(
    subject_id = integer(),
    group_id = integer(),
    from = integer(),
    to = integer(),
    dt_min = numeric(),
    dt_max = numeric(),
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
    dt_min = rep(0, n),
    dt_max = rep(Inf, n),
    t_sot = rep(NA, n)
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
      dt_min = pmax(1/30, dt_min), # Weibull CDF does not like exactly 0
      dt_max = pmax(dt_min + 1/30, dt_max)
    ) %>%
    as.list()

  lst_stan_data$M_groups <- length(attr(model, "group_id"))
  lst_stan_data$N <- nrow(data)

  lst_stan_data <- c(lst_stan_data, as.list(model))

  return(lst_stan_data)
}

#' @importFrom stringr str_extract
.parameter_sample_to_tibble.srp_model <- function(model, sample) {
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

