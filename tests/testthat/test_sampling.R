test_that("Testing that fixed seed works", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    B = define_srp_prior(),
    C = define_srp_prior()
  )
  smpl_prior1 <- sample_prior(mdl, seed = 36L, nsim = 500L)
  tbl1 <- sample_predictive(mdl,
    sample = smpl_prior1,
    n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423
  )
  smpl_prior2 <- sample_prior(mdl, seed = 36L, nsim = 500L)
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
        all(.$n == c(1, 1) & .$group_id == c("A", "C"))
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
  eps <- 0.01
  mdl <- create_srpmodel(
    A = define_srp_prior(
        p_mean = 0.5,
        p_n = 1e2,
        median_t_q05 = c(4, 6, 12) - eps,
        median_t_q95 = c(4, 6, 12) + eps,
        shape_q05 = c(1, 1, 1) - eps,
        shape_q95 = c(1, 1, 1) + eps,
        visit_spacing = 0.1
      ),
    maximal_time = 12 * 100
  )
  smpl_prior <- sample_prior(mdl, seed = 36L, nsim = 500L)
  tbl_prior_predictive <- sample_predictive(mdl,
    sample = smpl_prior,
    n_per_group = 1, nsim = 1e3,
    p = 0.5, # sample under fixed value
    seed = 948935435
  )
  # test that observed response rate is close enough to 0.5
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
  theoretical_means <- c(4, 6, 12) / log(2) * gamma(2)
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
    abs(estimated_mean - theoretical_mean) <= 3 * se
  )))
  # check for two groups ------------------------------------------------------
  mdl <- create_srpmodel(
    A = define_srp_prior(
      p_mean = 0.33, p_n = 1e2,
      visit_spacing = 0.1
    ),
    B = define_srp_prior(
      p_mean = 0.8, p_n = 1e2,
      visit_spacing = 0.1
    ),
    maximal_time = 12 * 100
  )
  smpl_prior <- sample_prior(mdl, seed = 36L, nsim = 500L)
  tbl <- sample_predictive(mdl,
    sample = smpl_prior,
    n_per_group = c(1L, 1L), nsim = 1e3, seed = 423
  )
  # test that observed response rate is close enough to true rate
  res_test <- tbl %>%
    filter(group_id == "A") %>%
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
    filter(group_id == "B") %>%
    group_by(iter, subject_id) %>%
    summarize(
      responder = any(state == "response"),
      .groups = "drop"
    ) %>%
    {
      binom.test(sum(.$responder), nrow(.), p = .8)
    } # nolint
  expect_true(res_test$p.value >= 0.01)
  # Testing marginal calibration of sampling from the fixed shape & scale -----
  mdl <- create_srpmodel(
    A = define_srp_prior(
        p_mean = 0.5, p_n = 1e2,
        median_t_q05 = c(4, 6, 12) - eps,
        median_t_q95 = c(4, 6, 12) + eps,
        shape_q05 = c(1, 1, 1) - eps,
        shape_q95 = c(1, 1, 1) + eps,
        visit_spacing = 0.1
      ),
    maximal_time = 100 * 12
  )
  smpl_prior <- sample_prior(mdl, seed = 36L, nsim = 500L)
  # define test function to estimate scale/shape from sampled data
  scale_default <- matrix(c(4, 6, 12), ncol = 1, nrow = 3)
  test_calibration <- function(scale_factor, shape) {
    # sample setting to different response probabilities
    tbl_prior_predictive <- sample_predictive(
      mdl,
      sample = smpl_prior,
      scale = scale_default * scale_factor,
      shape = matrix(shape, ncol = 1, nrow = 3),
      n_per_group = 1L, nsim = 1e3,
      seed = 342
    )
    # check calibration of times to next event, use midpoints of intervals as
    # approximation
    # work out the theoretical mean given scale = 1 and true scale
    # (see https://en.wikipedia.org/wiki/Weibull_distribution)
    theoretical_means <- scale_default * scale_factor * gamma(1 + 1 / shape)
    # check that stable to response timings are roughly calibrated,
    # use midpoints of intervals
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
    # testing for comparison with theoretical mean, allowing for
    # estimation error
    expect_true(all(with(
      tbl_means,
      abs(estimated_mean - theoretical_mean) <= 3 * se
    )))
  }
  # test over a grid of shape/scale combinations
  for (scale_factor in c(1, 3)) {
    for (shape in c(0.5, 1, 5)) {
      test_calibration(scale_factor, shape)
    }
  }
})



test_that("Testing if a single indivdual can be sampled only once in a group", { # nolint
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  smpl_prior <- sample_prior(mdl, seed = 36L)
  tbl <- sample_predictive(mdl,
    sample = smpl_prior, n_per_group = 1L,
    nsim = 1, seed = 42
  )
  # there should only be a single individual, single iteration
  expect_true(length(unique(tbl$subject_id)) == 1)
})


test_that("Testing sampling multiple individuals", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    B = define_srp_prior()
  )
  smpl_prior <- sample_prior(mdl, seed = 36L)
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



test_that("posterior shifts as expected", {
  # use very tight visitig spacing to avoid overlooking responses
  mdl <- create_srpmodel(
    A = define_srp_prior(
      p_n = 5, p_eta = .1,
      visit_spacing = 0.01, recruitment_rate = 100),
    B = define_srp_prior(
      p_n = 5, p_eta = .1,
      visit_spacing = 0.01, recruitment_rate = 100)
  )
  tbl_data <- sample_predictive(
    mdl,
    n_per_group = c(20, 20),
    nsim = 1,
    seed = 42L,
    p = c(0, 1),
    # use very small scale to avoid censoring issues
    scale = matrix(5, nrow = 2, ncol = 3)
  )
  p_obs <- tbl_data %>%
    group_by(group_id, iter, subject_id) %>%
    summarize(
      responder = any(state == "response"),
      .groups = "drop"
    ) %>%
    group_by(group_id) %>%
    summarize(
      p_response = mean(responder)
    ) %>%
    pull(p_response)
  smpl_posterior <- sample_posterior(mdl, tbl_data)
  p <- parameter_sample_to_tibble(mdl, smpl_posterior) %>%
    filter(parameter == "p") %>%
    group_by(group_id) %>%
    summarize(mean = mean(value)) %>%
    pull(mean)
  expect_true(all(c(
    p[1] < 0.1,
    p[2] > 0.9
  )))
})
