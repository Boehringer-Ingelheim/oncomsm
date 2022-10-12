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
  sample_predictive(mdl,
                    sample = smpl_prior,
                    n_per_group = c(30L, 30L),
                    nsim = 1,
                    seed = 3423423,
                    DEBUG = FALSE)
  expect_true(isa(mdl, c("srp_model", "Model", "list")))

})
