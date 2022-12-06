# plot() =======================================================================
test_that("default plotting works as intended", {
  mdl <- create_srpmodel(
    A = define_srp_prior(
      p_n = 10,
      median_t_q05 = c(1, 1, 3), median_t_q95 = c(3, 6, 12),
      shape_q05 = c(0.5, 3, 3), shape_q95 = c(1, 5, 5)
    ),
    B = define_srp_prior(
      p_mean = 0.75, p_n = 25, p_eta = 0.3
    )
  )
  # use stored sample to ensure cross-platform consistency
  # smpl <- sample_prior(mdl, seed = 1414322, nsim = 500L) # nolint
  # saveRDS(smpl, file = "inst/testthat/prior_smpl_plotting_1.rds", version = 3) # nolint
  smpl <- readRDS(
    system.file("testthat/prior_smpl_plotting_1.rds", package = "oncomsm")
  )
  plt <- plot(mdl, parameter_sample = smpl, confidence = 0.9, dt_n_grid = 5)
  vdiffr::expect_doppelganger("plot.srp_model_1", plt)

  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # smpl <- sample_prior(mdl, seed = 1414322, nsim = 500L) # nolint
  # saveRDS(smpl, file = "inst/testthat/prior_smpl_plotting_2.rds", version = 3) # nolint
  smpl <- readRDS(
    system.file("testthat/prior_smpl_plotting_2.rds", package = "oncomsm")
  )
  plt <- plot(mdl, parameter_sample = smpl, n_grid = 5)
  vdiffr::expect_doppelganger("plot.srp_model_2", plt)
  # check that plotting without sample works
  plot(mdl, dt_n_grid = 5, nsim = 250L)
  # check that plotting sub functions work (output is subsumed in plot())
  plot_pfs(mdl, dt_n_grid = 5, nsim = 250L)
  plot_transition_times(mdl, dt_n_grid = 5, nsim = 250L)
  plot_response_probability(mdl, nsim = 250L)
})



# plot_mstate() ================================================================
test_that("can plot mstate data for SRP model", {
  # one group
  mdl <- create_srpmodel(
    `1` = define_srp_prior()
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
  mdl <- create_srpmodel(
    `1` = define_srp_prior(),
    `2` = define_srp_prior()
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
