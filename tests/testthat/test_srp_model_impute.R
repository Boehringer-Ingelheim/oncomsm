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
      ), byrow = TRUE,  nrow = 3, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 3, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)

  tbl1 <- sample_predictive(mdl, sample = smpl_prior,
                            n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423)
  tbl2 <- sample_predictive(mdl, sample = smpl_prior,
                            n_per_group = c(1L, 0L, 1L), nsim = 20, seed = 423)
  # check that both tables are exactly the same
  expect_true(all(tbl1 == tbl2))
  # check n per group
  expect_true(
    tbl1 %>%
      group_by(group_id) %>%
      summarize(n = length(unique(subject_id))) %>%
      {all(.$n == c(1, 1) & .$group_id == c("1", "3"))} # nolint
  )
})



test_that("Testing marginal calibration of sampling from the prior", {
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .01, # low sd means we basically sample from a fixed parameter
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
        3, 3, 6
      ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01, byrow = TRUE, nrow = 1,
                                          ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  tbl_prior_predictive <- sample_predictive(mdl, sample = smpl_prior,
                                            n_per_group = 1L, nsim = 1000,
                                            seed = 342)
  # test that observed response rate is close enough to 0.5 (logodds(0.5) = 0)
  tbl_observed_rr <- tbl_prior_predictive %>%
    group_by(subject_id, iter) %>%
    summarise(
      responder = any(to == "response"),
      .groups = "drop"
    ) %>%
    summarise(
      p_hat = mean(responder),
      p_se = sd(responder) / sqrt(n())
    )
  # allow for 2 standard errors:
  expect_true(with(tbl_observed_rr, abs(p_hat - 0.5) < 5 * p_se))
  # check calibration of times to next event, use midpoints of intervals as
  # approximation
  # work out the theoretical mean given scale = 1 and specified median = 3
  # (see mdl definition and https://en.wikipedia.org/wiki/Weibull_distribution)
  theoretical_means <- mdl$median_time_to_next_event_mean / log(2) * gamma(2)
  # Theoretical mean time to event for response -> progression computed from
  # the initial stable state
  theoretical_means[3] <- theoretical_means[3] + theoretical_means[2]
  # check that stable to response timings are roughly calibrated, use midpoints
  # of intervals
  tbl_means <- tbl_prior_predictive %>%
    group_by(from, to) %>%
    summarise(
      mean_hat = mean((t_min + t_max) / 2),
      mean_se = sd((t_min + t_max) / 2) / sqrt(n()),
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
  expect_true(all(with(tbl_means,
    abs(mean_hat - theoretical_mean) <= 3 * mean_se
  )))
})



test_that("Testing marginal calibration of sampling from the fixed p", {
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .01, # low sd means we basically sample from a fixed parameter
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01, byrow = TRUE, nrow = 1,
                                          ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  # sample setting to different response probabilities
  for (p_true in c(0.1, 0.25, 0.75, 0.9)) {
    tbl_prior_predictive <- sample_predictive(mdl, sample = smpl_prior,
                                              p = p_true, n_per_group = 1L,
                                              nsim = 1e2, seed = 342)
    # test that observed response rate is close enough to true rate
    res <- tbl_prior_predictive %>%
      group_by(subject_id, iter) %>%
      summarise(
        responder = any(to == "response"),
        .groups = "drop"
      ) %>%
      {binom.test(x = sum(.$responder), n = nrow(.), p = p_true)} # nolint
    expect_true(res$p.value > 0.01)
  }
  # check for two groups
  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean = rep(0, 2),
    logodds_sd = rep(0.01, 2),
    visit_spacing = rep(1.2, 2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE,  nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01, byrow = TRUE,
                                          nrow = 2, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  tbl <- sample_predictive(mdl, sample = smpl_prior,
                           n_per_group = c(1L, 1L), nsim = 1e2, seed = 423,
                           p = c(0.1, 0.9))
  # test that observed response rate is close enough to true rate
  res <- tbl %>%
    group_by(group_id, subject_id, iter) %>%
    summarise(
      responder = any(to == "response"),
      .groups = "drop"
    ) %>%
    group_by(group_id) %>%
    tidyr::nest() %>%
    mutate(
      p_true = case_when(group_id == "1" ~ 0.1, group_id == "2" ~ 0.9)
    ) %>%
    mutate(
      p_value = purrr::map2_dbl(data, p_true,
                         ~binom.test(x = sum(..1$responder), n = nrow(..1),
                                     p = ..2)$p.value)
    )
  expect_true(all(res$p_value > 0.01))
})



test_that("Testing marginal calibration of sampling from the fixed scale", {
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .01,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01, byrow = TRUE, nrow = 1,
                                          ncol = 3)
  )
  # shape is implicitly set to 1
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 1000, seed = 36L)
  scale_true <- matrix(c(10, 10, 10), nrow = 1)
  # sample setting to different response probabilities
  tbl_prior_predictive <- sample_predictive(mdl, sample = smpl_prior,
                                            scale = scale_true,
                                            n_per_group = 1L, nsim = 1e3,
                                            seed = 342)
  # check calibration of times to next event, use midpoints of intervals as
  # approximation
  # work out the theoretical mean given scale = 1 and true scale
  # (see https://en.wikipedia.org/wiki/Weibull_distribution)
  theoretical_means <- scale_true * gamma(1 + 1 / 1) # shape is one
  # Theoretical mean time to event for response -> progression computed from
  # the initial stable state
  theoretical_means[3] <- theoretical_means[3] + theoretical_means[2]
  # check that stable to response timings are roughly calibrated, use midpoints
  # of intervals
  tbl_means <- tbl_prior_predictive %>%
    group_by(from, to) %>%
    summarise(
      mean_hat = mean((t_min + t_max) / 2),
      mean_se = sd((t_min + t_max) / 2) / sqrt(n()),
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
  expect_true(all(with(tbl_means,
                       abs(mean_hat - theoretical_mean) <= 3 * mean_se
  )))
})



test_that("Testing short time response->progression (within same visit interval)", { # nolint
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = 0.5,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 0.1, 6
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 1, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)
  tbl <- sample_predictive(mdl, sample = smpl_prior, n_per_group = 5L,
                           nsim = 100, seed = 3423423)
  # count how often the two jumps occur in the same visit interval
  # should not happen, since one would then combine them into a direct
  # transition between the two endpoints
  n_double_jumps <- tbl %>%
    group_by(iter, subject_id) %>%
    filter(t_min == lag(t_min) & t_max == lag(t_max)) %>%
    nrow()
  expect_true(n_double_jumps == 0)
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
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 1, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 250, nsim = 500, seed = 36L)
  tbl <- sample_predictive(mdl, sample = smpl_prior, n_per_group = 1L,
                           nsim = 1, seed = 42)
  # there should only be a single individual, single iteration
  expect_true(nrow(tbl) == 1)
})



test_that("can generate visit data from SRP model", {
  mdl <- create_srp_model(
    group_id = c("A", "B"),
    logodds_mean =  c(logodds(.5), logodds(.5)),
    logodds_sd = c(.1, .1),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE,  nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 2, ncol = 3)
  )
  generate_visit_data(mdl, n_per_group = c(20, 20),
                      recruitment_rate = c(.1, 5), seed = 112341)
  expect_true(TRUE) # TODO: implement check
})



test_that("impute remainder of trial from interim data", {
  mdl <- create_srp_model(
    group_id = c("A", "B"),
    logodds_mean =  c(logodds(.5), logodds(.5)),
    logodds_sd = c(.1, .1),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE,  nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 2, ncol = 3)
  )
  # sample some data
  tbl_data1 <- sample_predictive(mdl, c(20, 20), nsim = 1) %>%
    select(-iter)
  # impute another 20/group conditional on observed data
  tbl_data2 <- impute(mdl, tbl_data1, c(40, 40),
                      recruitment_rates = c(1, 1), nsim = 25)
  expect_true(
    tbl_data2 %>% group_by(subject_id, group_id, iter) %>% n_groups() ==
      40 * 25 * 2
  )
})
