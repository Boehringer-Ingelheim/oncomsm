test_that("Testing that fixed seed works", {
  mdl <- create_srp_model(
    group_id = 1:3,
    logodds_mean = rep(0, 3),
    logodds_sd = rep(0.5, 3),
    visit_spacing = rep(1.2, 3),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE, nrow = 3, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 3, ncol = 3)
  )
  smpl_prior1 <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)
  tbl1 <- sample_predictive(mdl,
    sample = smpl_prior1,
    n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423
  )
  smpl_prior2 <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)
  tbl2 <- sample_predictive(mdl,
    sample = smpl_prior1,
    n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423
  )
  # check that both tables are exactly the same
  expect_true(all(tbl1 == tbl2))
  # check n per group
  expect_true(
    tbl1 %>%
      group_by(group_id) %>%
      summarize(n = length(unique(subject_id))) %>%
      {
        all(.$n == c(1, 1) & .$group_id == c("1", "3"))
      }
  )
  # sample directly
  tbl1 <- sample_predictive(mdl,
                            n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423
  )
  tbl2 <- sample_predictive(mdl,
                            n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423
  )
  expect_true(all(tbl1 == tbl2))
})



test_that("Testing marginal calibration of sampling from the prior", {
  # to test that the response probability is properly calibrated, all factors
  # need to be taken into account.
  # The discrete visit spacing can lead to missing a few responses if the
  # transition to progression happens within a visit interval.
  # To minimize this effect, a shape > 1, small visit interval, and long
  # time from response->progression is used
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .01, # low sd means we basically sample from a fixed parameter
    visit_spacing = 0.1,
    median_time_to_next_event = matrix(c(
      3, 3, 24
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01,
      byrow = TRUE, nrow = 1,
      ncol = 3
    ),
    max_time = 12 * 50 # large t to avoid limiting effects
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  tbl_prior_predictive <- sample_predictive(mdl,
    sample = smpl_prior,
    n_per_group = 1, nsim = 1e3,
    p = 0.5, # sample under fixed value
    seed = 6563
  )
  # test that observed response rate is close enough to 0.5 (logodds(0.5) = 0)
  res_test <- tbl_prior_predictive %>%
    group_by(iter, subject_id) %>%
    summarize(
      responder = any(state == "response"),
      .groups = "drop"
    ) %>%
    {
      binom.test(sum(.$responder), nrow(.), p = 0.5)
    } # nolint
  expect_true(res_test$p.value >= 0.01)
  # check calibration of times to next event, use midpoints of intervals as
  # approximation
  # work out the theoretical mean given scale = 1 and specified median = 3
  # (see mdl definition and https://en.wikipedia.org/wiki/Weibull_distribution)
  theoretical_means <- mdl$median_time_to_next_event_mean / log(2) * gamma(2)
  # check that stable to response timings are roughly calibrated
  tbl_means <- tbl_prior_predictive %>%
    distinct(subject_id, iter, state, .keep_all = TRUE) %>%
    group_by(iter, group_id, subject_id) %>%
    summarize(
      dt = t - lag(t),
      from = lag(state),
      to = state,
      .groups = "drop"
    ) %>%
    filter(to != "stable") %>%
    group_by(group_id, from, to) %>%
    summarize(
      estimated_mean = mean(dt),
      se = sd(dt) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      theoretical_mean = case_when(
        from == "stable" & to == "response" ~ theoretical_means[1],
        from == "stable" & to == "progression" ~ theoretical_means[2],
        from == "response" & to == "progression" ~ theoretical_means[3],
      )
    )
  # testing for comparison with theoretical mean, allowing for estimation error
  expect_true(all(with(
    tbl_means,
    abs(estimated_mean - theoretical_mean) <= 2 * se
  )))
  # check for two groups ------------------------------------------------------
  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean = logodds(c(.33, .8)),
    logodds_sd = rep(0.1, 2),
    visit_spacing = rep(0.2, 2),
    max_time = 12 * 50,
    median_time_to_next_event = matrix(c(
      3, 3, 12,
      3, 3, 12
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1,
      byrow = TRUE,
      nrow = 2, ncol = 3
    )
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  tbl <- sample_predictive(mdl,
    sample = smpl_prior,
    n_per_group = c(1L, 1L), nsim = 1e3, seed = 423
  )
  # test that observed response rate is close enough to true rate
  res_test <- tbl %>%
    filter(group_id == "1") %>%
    group_by(iter, subject_id) %>%
    summarize(
      responder = any(state == "response"),
      .groups = "drop"
    ) %>%
    {
      binom.test(sum(.$responder), nrow(.), p = 0.33)
    } # nolint
  expect_true(res_test$p.value >= 0.01)
  res_test <- tbl %>%
    filter(group_id == "2") %>%
    group_by(iter, subject_id) %>%
    summarize(
      responder = any(state == "response"),
      .groups = "drop"
    ) %>%
    {
      binom.test(sum(.$responder), nrow(.), p = .8)
    } # nolint
  expect_true(res_test$p.value >= 0.01)
  # Testing marginal calibration of sampling from the fixed scale -------------
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .1,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01,
      byrow = TRUE, nrow = 1,
      ncol = 3
    )
  )
  # shape is implicitly set to 1
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  scale_true <- matrix(c(10, 10, 10), nrow = 1)
  # sample setting to different response probabilities
  tbl_prior_predictive <- sample_predictive(mdl,
    sample = smpl_prior,
    scale = scale_true,
    n_per_group = 1L, nsim = 1e3,
    seed = 342
  )
  # check calibration of times to next event, use midpoints of intervals as
  # approximation
  # work out the theoretical mean given scale = 1 and true scale
  # (see https://en.wikipedia.org/wiki/Weibull_distribution)
  theoretical_means <- scale_true * gamma(1 + 1 / 1) # shape is one
  # check that stable to response timings are roughly calibrated, use midpoints
  # of intervals
  tbl_means <- tbl_prior_predictive %>%
    distinct(subject_id, iter, state, .keep_all = TRUE) %>%
    group_by(iter, group_id, subject_id) %>%
    summarize(
      dt = t - lag(t),
      from = lag(state),
      to = state,
      .groups = "drop"
    ) %>%
    filter(to != "stable") %>%
    group_by(group_id, from, to) %>%
    summarize(
      estimated_mean = mean(dt),
      se = sd(dt) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      theoretical_mean = case_when(
        from == "stable" & to == "response" ~ theoretical_means[1],
        from == "stable" & to == "progression" ~ theoretical_means[2],
        from == "response" & to == "progression" ~ theoretical_means[3],
      )
    )
  # testing for comparison with theoretical mean, allowing for estimation error
  expect_true(all(with(
    tbl_means,
    abs(estimated_mean - theoretical_mean) <= 2 * se
  )))
})



test_that("Testing if a single indivdual can be sampled only once in a group", { # nolint
  mdl <- create_srp_model(
    group_id = 1,
    # very low chance of response -> direct progression
    logodds_mean = logodds(0.001),
    logodds_sd = 0.5,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 6, 9
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 1, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)
  tbl <- sample_predictive(mdl,
    sample = smpl_prior, n_per_group = 1L,
    nsim = 1, seed = 42
  )
  # there should only be a single individual, single iteration
  expect_true(length(unique(tbl$subject_id)) == 1)
})


test_that("Testing sampling multiple individuals", {
  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean = logodds(c(.33, .8)),
    logodds_sd = rep(0.1, 2),
    visit_spacing = rep(0.2, 2),
    max_time = 12 * 50,
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(0.1,
      byrow = TRUE,
      nrow = 2, ncol = 3
    )
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)
  tbl <- sample_predictive(mdl,
    sample = smpl_prior, n_per_group = c(20L, 40L),
    nsim = 5, seed = 42
  )
  tmp <- tbl %>%
    distinct(group_id, subject_id, iter, .keep_all = TRUE) %>%
    group_by(group_id, subject_id) %>%
    summarize(niter = n(), .groups = "drop")
  expect_true(nrow(tmp) == 60L)
  expect_true(all(tmp$niter == 5L))
})



test_that("impute remainder of trial from interim data", {
  mdl <- create_srp_model(
    group_id = c("A", "B"),
    logodds_mean = c(0, 0),
    logodds_sd = c(.1, .1),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 2, ncol = 3)
  )
  # sample some data
  tbl_data1 <- sample_predictive(mdl, c(20, 20),
    nsim = 1, seed = 43L,
    nsim_parameters = 1500L
  )
  # impute another 20/group conditional on observed data
  tbl_data2 <- impute(mdl, tbl_data1, c(40, 40), nsim = 25, seed = 4453L)
  expect_true(
    tbl_data2 %>%
      group_by(subject_id, group_id, iter) %>%
      n_groups() == 40 * 25 * 2
  )
})
