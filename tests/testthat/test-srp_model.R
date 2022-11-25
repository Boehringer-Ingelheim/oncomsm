test_that("can create SRP model", {
  mdl <<- create_srp_model(
    group_id = c("A", "B"),
    logodds_mean = c(logodds(.5), logodds(.5)),
    logodds_sd = c(.1, .1),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 2, ncol = 3)
  )
  # check class
  expect_true(isa(mdl, c("srp_model", "Model", "list")))
  # check print method
  expect_true(capture.output(print(mdl)) == "srp_model<A,B> ")

  mdl <- create_srp_model(
    group_id = "A",
    logodds_mean = 0,
    logodds_sd = 0.5,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 1, ncol = 3)
  )
  # check print method
  expect_true(capture.output(print(mdl)) == "srp_model<A> ")
}) # "can create SRP model"



test_that("private function is_valid throws correct errors", {
  expect_error(
    create_srp_model(
      group_id = 1:2,
      logodds_mean = c(logodds(0.7), logodds(.5), logodds(.6)),
      logodds_sd = c(.5, .5),
      visit_spacing = c(1.2, 1.2),
      median_time_to_next_event = matrix(c(
        3, 3, 6,
        3, 3, 6
      ), byrow = TRUE, nrow = 2, ncol = 3),
      median_time_to_next_event_sd = matrix(
        1,
        byrow = TRUE, nrow = 2, ncol = 3
      )
    ),
    regexp = "Assertion on 'logodds_mean' failed: Must have length 2, but has length 3" # nolint
  )

  expect_error(
    create_srp_model(
      group_id = 1,
      logodds_mean = c(logodds(.9)),
      logodds_max = c(logodds(.75)),
      logodds_sd = c(.5),
      visit_spacing = c(1.2),
      median_time_to_next_event = matrix(c(
        3, 3, 6
      ), byrow = TRUE, nrow = 1, ncol = 3),
      median_time_to_next_event_sd = matrix(
        1,
        byrow = TRUE, nrow = 1, ncol = 3
      )
    ),
    regexp = "Assertion on 'logodds_mean < logodds_max' failed: Must be TRUE"
  )

  expect_error(
    create_srp_model(
      group_id = 1:2,
      logodds_mean = c(logodds(0.25), logodds(.5)),
      logodds_sd = c(.5, .5),
      visit_spacing = c(1.2, 1.2),
      median_time_to_next_event = matrix(c(
        3, -3, 6,
        3, 3, 6
      ), byrow = TRUE, nrow = 2, ncol = 3),
      median_time_to_next_event_sd = matrix(
        1,
        byrow = TRUE, nrow = 2, ncol = 3
      )
    ),
    regexp = "Assertion on 'median_time_to_next_event_mean' failed: Element 1 is not >= 0" # nolint
  )
})



test_that("can convert data to standata for SRP model", {
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
  lst_standata <- oncomsm:::data2standata.srp_model(tbl_mstate, mdl)
  expect_true(TRUE) # TODO: implement check
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



test_that("posterior shifts as expected", {
  tbl_data <- sample_predictive(
    mdl,
    n_per_group = c(20, 20),
    nsim = 1,
    seed = 42L,
    p = c(0.1, 0.9)
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
  expect_true(all(case_when(
    p_obs > 0.5 & p > 0.5 ~ TRUE,
    p_obs < 0.5 & p < 0.5 ~ TRUE,
    TRUE ~ FALSE
  )))
})
