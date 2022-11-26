# plot() =======================================================================
test_that("default plotting works as intended", {
  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean = c(logodds(.5), logodds(.75)),
    logodds_sd = c(.5, .75),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      4, 6, 10
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 2, ncol = 3)
  )
  # use stored sample to ensure cross-platform consistency
  # smpl <- sample_prior(mdl, nsim = 500, seed = 1414322) # nolint
  # saveRDS(smpl, file = "inst/testthat/prior_smpl_plotting_1.rds", version = 3) # nolint
  smpl <- readRDS(
    system.file("testthat/prior_smpl_plotting_1.rds", package = "oncomsm")
  )
  plt <- plot(mdl, parameter_sample = smpl, n_grid = 10)
  vdiffr::expect_doppelganger("plot.srp_model_1", plt)

  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = 0,
    logodds_sd = .5,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 1, ncol = 3)
  )
  # smpl <- sample_prior(mdl, nsim = 500, seed = 1414322) # nolint
  # saveRDS(smpl, file = "inst/testthat/prior_smpl_plotting_2.rds", version = 3) # nolint
  smpl <- readRDS(
    system.file("testthat/prior_smpl_plotting_2.rds", package = "oncomsm")
  )
  plt <- plot(mdl, parameter_sample = smpl, n_grid = 10)
  vdiffr::expect_doppelganger("plot.srp_model_2", plt)
  # check that plotting without sample works
  plot(mdl)
})



# plot_mstate() ================================================================
test_that("can plot mstate data for SRP model", {

  # one group
  mdl <- create_srp_model(
    group_id = 1,
    logodds_mean = logodds(.5),
    logodds_sd = 0.5,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 1, ncol = 3)
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
    "1", "1", 0, "stable",
    "1", "1", 1.2, "stable",
    "1", "1", 2.4, "response",
    "1", "1", 3.0, "response",
    "1", "1", 3.6, "progression",
    "1", "2", 1.5, "stable",
    "1", "2", 2, "stable",
    "1", "2", 3, "response",
    "1", "3", 1.5, "stable",
    "1", "3", 3, "response",
    "1", "3", 4, "response",
    "1", "3", 4.25, "EOF"
  )

  tbl_mstate <- visits_to_mstate(tbl_visits, mdl)
  plt <- plot_mstate(tbl_mstate, mdl, relative_to_sot = FALSE)
  vdiffr::expect_doppelganger("plot_mstate.srp_model_1", plt)

  # two groups
  mdl <- create_srp_model(
    group_id = 1:2,
    logodds_mean = rep(logodds(.5), 2),
    logodds_sd = rep(0.5, 2),
    visit_spacing = rep(1.2, 2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 2, ncol = 3)
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
    "1", "1", 0, "stable",
    "1", "1", 1.2, "stable",
    "1", "1", 2.4, "response",
    "1", "1", 3.0, "response",
    "1", "1", 3.6, "progression",
    "1", "2", 1.5, "stable",
    "1", "2", 2, "stable",
    "1", "2", 3, "response",
    "1", "3", 1.5, "stable",
    "2", "3", 3, "response",
    "2", "3", 4, "response",
    "2", "3", 4.25, "EOF"
  )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl)
  plt <- plot_mstate(tbl_mstate, mdl, relative_to_sot = FALSE)
  vdiffr::expect_doppelganger("plot_mstate.srp_model_2", plt)
  # check relative to SOT
  plt <- plot_mstate(tbl_mstate, mdl, relative_to_sot = TRUE)
  vdiffr::expect_doppelganger("plot_mstate.srp_model_3", plt)
})
