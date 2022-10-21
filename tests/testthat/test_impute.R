test_that("Testing .impute function", {

  mdl <- create_srp_model(
    group_id = 1:3,
    logodds_mean =  c(logodds(.5), logodds(.5), logodds(.5)),
    logodds_sd = c(.5, .5, .5),
    visit_spacing = c(1.2, 1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE,  nrow = 3, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 3, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 500, nsim = 2000, seed = 36L)

  # sample only from one group
  tbl1 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 1L, 0L),
                            nsim = 20,
                            seed = 3423423)
  tbl2 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 1L, 0L),
                            nsim = 20,
                            seed = 3423423)
  expect_true(all(tbl1 == tbl2))
  # TODO: Take care of the corner case where progression and response intervals
  # are the same
  # TODO: Write a test case with low median time, so that many corner cases are
  # generated
  mdl2 <- create_srp_model(
    group_id = 1,
    logodds_mean =  c(logodds(.5)),
    logodds_sd = c(.5),
    visit_spacing = c(1.2),
    median_time_to_next_event = matrix(c(
      0.5, 0.7, 0.6
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 1, ncol = 3)
  )
  smpl_prior2 <- sample_prior(mdl2, warmup = 500, nsim = 2000, seed = 36L)

  tbl3 <- sample_predictive(mdl2,
                            sample = smpl_prior2,
                            n_per_group = c(1L),
                            nsim = 3,
                            seed = 3423423)
  # TODO: check that only one individual from group one is sampled once
  tbl4 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 2L, 1L),
                            nsim = 2,
                            seed = 3423423,
                            debug = TRUE)
  tbl_data <- tbl4 %>%
    filter(group_id == 1) %>%
    count(subject_id)
  expect_true(tbl_data["subject_id"] == 1)

})
