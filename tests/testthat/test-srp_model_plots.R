test_that("default plotting works as intended", {

  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean =  c(logodds(.5), logodds(.75)),
    logodds_sd = c(.5, .75),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      4, 6, 10
    ), byrow = TRUE,  nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 2, ncol = 3)
  )
  # use stored sample to ensure cross-platform consistency
  # smpl <- sample_prior(mdl, nsim = 500, seed = 1414322)
  # saveRDS(smpl, file = "inst/testthat/prior_smpl_plotting_1.rds", version = 3)
  smpl <- readRDS(
    system.file("testthat/prior_smpl_plotting_1.rds", package = "oncomsm")
  )
  plt <- plot(mdl, dt = c(0, 36), sample = smpl, n_grid = 10)
  vdiffr::expect_doppelganger("plot.srp_model_1", plt)

  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .5,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 1, ncol = 3)
  )
  # smpl <- sample_prior(mdl, nsim = 500, seed = 1414322)
  # saveRDS(smpl, file = "inst/testthat/prior_smpl_plotting_2.rds", version = 3)
  smpl <- readRDS(
    system.file("testthat/prior_smpl_plotting_2.rds", package = "oncomsm")
  )
  plt <- plot(mdl, dt = c(0, 36), sample = smpl, n_grid = 10)
  vdiffr::expect_doppelganger("plot.srp_model_2", plt)

})
