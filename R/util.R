#' Log-odds function
#'
#' computed the log odds of a probability.
#'
#' @param p numeric of probabilities
#'
#' @return log(p/(1-p))
#'
#' @export
logodds <- function(p) log(p / (1 - p))



# create_example_data
.generate_example_data <- function() {
  mdl <- create_srp_model(
    # names of the arms/groups
    group_id = c("A", "B", "C"),
    # per-group logodds of response|stable
    logodds_mean = c(logodds(.20), logodds(.3), logodds(.3)),
    logodds_sd = c(.5, .5, .5),
    # m[i,j] is the median time to next event for group i and transition j
    median_time_to_next_event = matrix(c(
      3, 2, 6,
      2, 8, 9,
      2, 6, 24
    ), byrow = TRUE, nrow = 3, ncol = 3),
    # fixed standard deviation of the prior for all median times
    median_time_to_next_event_sd = matrix(
      1,
      byrow = TRUE, nrow = 3, ncol = 3
    ),
    # uniform prior over the shape parameter, difficult to identify,
    # better keep it tight to avoid issues with the sampler
    shape_min = matrix(
      .75,
      byrow = TRUE, nrow = 3, ncol = 3
    ),
    shape_max = matrix(
      2,
      byrow = TRUE, nrow = 3, ncol = 3
    ),
    # the visit interval
    visit_spacing = c(1.2, 1.2, 1.2)
  )
  recruitment_rate_overall <- 3
  sample_predictive(
      mdl,
      sample = smpl_prior,
      n_per_group = c(30L, 30L, 30L),
      nsim = 1,
      seed = 3423423
    ) %>%
    select(-t_sot) %>%
    {
      left_join(
        .,
        select(., subject_id, group_id) %>%
          distinct() %>%
          arrange(runif(n())) %>% # permute groups
          mutate(
            # poisson recruitment process
            t_sot = cumsum(rexp(n = n(), rate = recruitment_rate_overall)),
          ),
        by = c("subject_id", "group_id")
      )
    } %>%
    mutate(
      t_min = t_min + t_sot,
      t_max = t_max + t_sot
    ) %>%
    mstate_to_visits(mdl, .) %>%
    filter(t <= 10) %>%
    visits_to_mstate(
      start_state = "stable",
      absorbing_states = "progression",
      now = 10
    )
}
