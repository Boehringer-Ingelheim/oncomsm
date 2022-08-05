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

  lst_standata <- bhmbasket.predict:::data2standata.srp_model(
    mdl, bhmbasket.predict:::.nodata.srp_model(mdl)
  )
  expect_true(TRUE) # TODO: implement check

})



test_that("can convert data to standata for SRP model", {

  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,    ~t,        ~state,
            1,           1,     0,      "stable",
            1,           1,   1.2,      "stable",
            1,           1,   2.4,    "response",
            1,           1,   3.6, "progression",
            1,           2,     2,      "stable",
            1,           2,     3,    "response",
            1,           3,   1.5,      "stable",
            1,           3,     3,    "response",
            1,           3,   3.5,         "EOF"
  )
  tbl_mstate <- visits_to_mstate(
    tbl_visits,
    start_state = "stable",
    absorbing_states = c("progression")
  )
  lst_standata <- bhmbasket.predict:::data2standata.srp_model(mdl, tbl_mstate)
  expect_true(TRUE) # TODO: implement check

})



test_that("can sample from prior", {

  tbl_prior_sample <- sample_prior(mdl, seed = 1414322)

  p <- tbl_prior_sample %>%
    filter(parameter == "p") %>%
    group_by(group_id) %>%
    summarize(mean = mean(value)) %>%
    pull(mean)

  expect_true(all(abs(p - 0.5) < .01))

})



test_that("can generate data from SRP model and seed works", {

  tbl_prior_predictive1 <<- sample_prior_predictive(mdl, c(10, 10), 10, seed = 42L)
  tbl_prior_predictive2 <- sample_prior_predictive(mdl, c(10, 10), 10, seed = 42L)

  expect_true(
    tbl_prior_predictive1 %>%
      group_by(group_id) %>%
      summarize(n = length(unique(subject_id))) %>%
      pull(n) %>%
      {. == c(10, 10)} %>%
      all()
  )
  expect_true(all(tbl_prior_predictive1 == tbl_prior_predictive2))

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

  tbl_posterior_sample <- sample_posterior(mdl, tbl_data)

  p <- tbl_posterior_sample %>%
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

  impute_posterior_predictive(mdl, tbl_data, nsim = 10)
  expect_true(TRUE) # TODO: implement check

})
