test_that("Testing prior predictive impute function", {
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean =  c(logodds(.5)),
    logodds_sd = c(.01),
    visit_spacing = c(0.1),
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(0.01,
                                          byrow = TRUE,  nrow = 1, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 500, nsim = 2000, seed = 36L)

  tbl_cpp <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L),
                            nsim = 10000,
                            seed = 3423423)
  tbl_cpp %>%
    group_by(subject_id, iter) %>%
    summarise(
      responder = any(to == "response"),
      .groups = "drop"
    ) %>%
    summarise(
      p_hat = mean(responder),
      p_se = sd(responder) / sqrt(n())
    )

  means <- tbl_cpp %>%
    filter(from == "stable", to == "response") %>%
    summarise(
      mean_hat = mean((t_min + t_max) / 2),
      mean_se = sd((t_min + t_max) / 2) / sqrt(n())
    )
  expect_true(all((means$mean_hat - 4.33 * gamma(2)) < 0.1))

})
