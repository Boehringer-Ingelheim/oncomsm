test_that("Testing .impute function", {

  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean =  c(logodds(.5), logodds(.5)),
    logodds_sd = c(.5, .5),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE,  nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 2, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 500, nsim = 2000, seed = 36L)

  # sample only from one group
  tbl1 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 0L),
                            nsim = 1,
                            seed = 3423423,
                            debug = FALSE)
  tbl2 <- sample_predictive(mdl,
                    sample = smpl_prior,
                    n_per_group = c(1L, 0L),
                    nsim = 1,
                    seed = 3423423,
                    debug = FALSE)
  # TODO: check that fixed seed works and results are the same
  # TODO: check that only one individual from group one is sampled once
  # TODO: temporary tests during refactoring from R to cpp implementation
  tbl3 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 0L),
                            nsim = 1,
                            seed = 3423423,
                            debug = TRUE)
  # TODO: check that c++ implementation returns the same as r implementation
  expect_true(all(dim(tbl3) == dim(tbl2)))
  expect_true(all(tbl3 == tbl2))

})
