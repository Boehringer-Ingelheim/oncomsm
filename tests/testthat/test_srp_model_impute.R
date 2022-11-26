test_that("impute remainder of trial from interim data", {
  mdl <- create_srp_model(
    group_id = c("A", "B"),
    logodds_mean = c(0, 0),
    logodds_sd = c(.1, .1),
    visit_spacing = c(1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE, nrow = 2, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE, nrow = 2, ncol = 3)
  )
  # sample some data
  tbl_data1 <- sample_predictive(mdl, c(20, 20),
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
  mdl <- create_srp_model(
    group_id = "A",
    logodds_mean = 0,
    logodds_sd = 10,
    visit_spacing = 1.2,
    median_time_to_next_event = matrix(c(
      3, 3, 6
    ), byrow = TRUE, nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(10, byrow = TRUE, nrow = 1, ncol = 3)
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
