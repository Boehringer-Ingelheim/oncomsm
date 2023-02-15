test_that("impute remainder of trial from interim data", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    B = define_srp_prior()
  )
  # sample some data
  tbl_data1 <- sample_predictive(
    mdl,
    c(20, 20),
    nsim = 1, seed = 43L,
    nsim_parameters = 1500L
  )
  # impute another 20/group conditional on observed data
  tbl_data2 <- impute(mdl, tbl_data1, c(40, 40), nsim = 25, seed = 4453L)
  expect_true(
    tbl_data2 %>%
      group_by(subject_id, group_id, iter) %>%
      n_groups() == 40 * 25 * 2
  )
})



test_that("impute remainder of trial from interim data, no new individuals", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # sample some data and reduce to first visits
  tbl_data1 <- sample_predictive(mdl, 10, nsim = 1, seed = 43L,
                                 nsim_parameters = 1500L
    ) %>%
    group_by(subject_id) %>%
    filter(row_number() == 1) %>%
    ungroup()
  # impute conditional on observed data,
  # need to adapt sampler a bit here to avoid diverging transitions!
  tbl_data2 <- impute(mdl,
                      tbl_data1, 10, nsim = 1, seed = 4453L)
  expect_true(
    n_groups(group_by(tbl_data2, subject_id, group_id)) == 10L
  )
})



test_that("impute remainder, without adding subjects", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # sample some data and reduce to first visits
  tbl_data1 <- sample_predictive(mdl, 10, nsim = 1, seed = 43L,
                                 nsim_parameters = 1500L
    ) %>%
    group_by(subject_id) %>%
    filter(row_number() == 1) %>%
    ungroup()
  # impute conditional on observed data,
  # need to adapt sampler a bit here to avoid diverging transitions!
  tbl_data2 <- impute(mdl,
                      tbl_data1, nsim = 1, seed = 4453L)
  expect_true(
    n_groups(group_by(tbl_data2, subject_id, group_id)) == 10L
  )
})



test_that("impute to mstate works", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # sample some data and reduce to first visits
  tbl_data <- sample_predictive(mdl, 10, nsim = 1, seed = 55L,
                                as_mstate = TRUE)
  expect_true(
    n_groups(group_by(tbl_data, subject_id, group_id)) == 10L
  )
})



test_that("impute throws correct errors", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # create some data in wrong format
  tbl_data <- tribble(
    ~subject_id, ~group_id, ~t,     ~state,
          " s1",       "A",  0, "response",
  )
  expect_error(
    impute(mdl,
           data = tbl_data, n_per_group = 1, nsim = 1, seed = 32,
           nuts_control = list(adapt_delta = 0.99))
  )
  # create some data in wrong format
  tbl_data <- tribble(
    ~subject_id, ~group_id, ~t,     ~state,
    " s1",       "A",  2, "stable",
    " s1",       "B",  1, "stable"
  )
  expect_error(
    impute(mdl,
           data = tbl_data, n_per_group = 1, nsim = 1, seed = 32,
           nuts_control = list(adapt_delta = 0.99))
  )
  # create some data in wrong format
  tbl_data <- tribble(
    ~subject_id, ~group_id, ~t,     ~state,
    " s1",       "A",  2, "stable",
    " s2",       "B",  1, "response"
  )
  expect_error(
    impute(mdl,
           data = tbl_data, n_per_group = 1, nsim = 1, seed = 32,
           nuts_control = list(adapt_delta = 0.99))
  )
})



test_that("impute from 'response' works", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  tbl_data <- tribble(
    ~subject_id, ~group_id, ~t,     ~state,
           "s1",       "A",  0,   "stable",
           "s1",       "A",  0, "response"
  )
  tbl_test <- impute(mdl,
       data = tbl_data, n_per_group = 1, nsim = 1, seed = 32,
     )
  expect_true(
    all(tbl_test$subject_id == "s1"),
    all(tbl_test$state[2:(nrow(tbl_test) - 1)] == "response"),
    tbl_test$state[nrow(tbl_test)] == "progression"
  )
})
