test_that("can generate visit data", {

  set.seed(42L)

  tbl_data <- generate_visit_data(
    group_id = c("A", "B", "c"),
    n = c(40, 30, 30),
    response_rate = c(0.33, .4, .5),
    recruitment_rate = c(2, 1, 1),
    visit_spacing = c(1.25, 1.25, 1.25),
    max_duration = 48,
    responder_weibull_scale = rep(4, 3),
    responder_weibull_shape = rep(5, 3),
    nonresponder_weibull_scale = rep(4, 3),
    nonresponder_weibull_shape = rep(2, 3)
  )

  expect_true(all(
    names(tbl_data) == c("group_id", "subject_id", "t", "status", "eof")
  ))

  first_statuses <- tbl_data %>%
    group_by(subject_id) %>%
    summarize(status = first(status)) %>%
    pull(status)
  expect_true(
      all(first_statuses == "S")
  )

  expect_true(max(tbl_data$t) <= 48)

})
