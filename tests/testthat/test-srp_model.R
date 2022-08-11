test_that("can create SRP model", {

  mdl <<- create_srp_model(
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
  expect_true(isa(mdl, c("srp_model", "Model", "list")))

}) # "can create SRP model"



test_that("can create empty standata for SRP model", {

  lst_standata <- oncomsm:::data2standata.srp_model(
    mdl, oncomsm:::.nodata.srp_model(mdl)
  )
  # TODO: implement check

})



test_that("can convert data to standata for SRP model", {

  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,    ~t,        ~state,
          "1",         "1",     0,      "stable",
          "1",         "1",   1.2,      "stable",
          "1",         "1",   2.4,    "response",
          "1",         "1",   3.0,    "response",
          "1",         "1",   3.6, "progression",
          "1",         "2",   1.5,      "stable",
          "1",         "2",     2,      "stable",
          "1",         "2",     3,    "response",
          "1",         "3",   1.5,      "stable",
          "1",         "3",     3,    "response",
          "1",         "3",     4,    "response",
          "1",         "3",   4.25,         "EOF"
  )
  tbl_mstate <<- visits_to_mstate(
    tbl_visits,
    start_state = "stable",
    absorbing_states = c("progression")
  )
  lst_standata <- oncomsm:::data2standata.srp_model(mdl, tbl_mstate)
  # TODO: implement check

})



test_that("can plot mstate data for SRP model", {

  create_plot <- function() plot_mstate(mdl, tbl_mstate, relative_to_sot = FALSE)

  vdiffr::expect_doppelganger("plot_mstate.srp_model", create_plot)

})



test_that("can sample from prior", {

  smpl_prior <- sample_prior(mdl, seed = 1414322)

  p <- parameter_sample_to_tibble(mdl, smpl_prior) %>%
    filter(parameter == "p") %>%
    group_by(group_id) %>%
    summarize(mean = mean(value)) %>%
    pull(mean)

  expect_true(all(abs(p - 0.5) < .01))

})



test_that("can generate data from SRP model", {

  tbl_prior_predictive1 <<- sample_predictive(
    mdl,
    n_per_group = c(20, 20),
    nsim = 25,
    seed = 42L
  )

  expect_true(
    tbl_prior_predictive1 %>%
      group_by(group_id) %>%
      summarize(n = length(unique(subject_id))) %>%
      pull(n) %>%
      {. == c(20, 20)} %>%
      all()
  )

  p <- tbl_prior_predictive1 %>% # we throw together iterations to get better estimate
    filter(from == "stable") %>%
    group_by(group_id, to) %>%
    summarize(n = length(unique(subject_id)), .groups = "drop_last") %>%
    mutate(p = n / sum(n)) %>%
    ungroup() %>%
    arrange(group_id, to) %>%
    pull(p)

  expect_true(all(abs(p - 0.5) < .01))

})



test_that("prior predictive seed works", {

  tbl_prior_predictive2 <-sample_predictive(
    mdl,
    n_per_group = c(20, 20),
    nsim = 25,
    seed = 42L
  )
  expect_true(all(tbl_prior_predictive1 == tbl_prior_predictive2))

})



test_that("can generate visit data from SRP model", {

  generate_visit_data(mdl, n_per_group = c(20, 20), seed = 112341)

  # TODO: plausi checks on data

})



test_that("can sample from posterior", {

  tbl_data <- tbl_prior_predictive1 %>%
    filter(iter == 1) %>%
    select(-iter)

  p_obs <- tbl_data %>%
    filter(from == "stable") %>%
    group_by(group_id) %>%
    summarize(p = mean(to == "response")) %>%
    pull(p)

  smpl_posterior <- sample_posterior(mdl, tbl_data)

  p <- parameter_sample_to_tibble(mdl, smpl_posterior) %>%
    filter(parameter == "p") %>%
    group_by(group_id) %>%
    summarize(mean = mean(value)) %>%
    pull(mean)

  expect_true(all(case_when(
    p_obs > 0.5 & p > 0.5 ~ TRUE,
    p_obs < 0.5 & p < 0.5 ~ TRUE,
    TRUE ~ FALSE
  )))

})



test_that("can sample from posterior predictive", {

  tbl_data <- tbl_prior_predictive1 %>%
    filter(iter == 1) %>%
    select(-iter)

  impute_predictive(mdl, data = tbl_data, nsim = 10)
  expect_true(TRUE) # TODO: implement check

})
