test_that("visits to mstate conversion works", {

    # base case
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,    ~t,   ~state,
            1,           1,     0, "stable"
  )
  tbl_mstate <- visits_to_mstate(
    tbl_visits,
    start_state = "stable",
    absorbing_states = c("progression"),
    now = 0
  )
  expect_true(with(tbl_mstate, all(c(
     subject_id == 1,
     group_id == 1,
     from == "stable",
     is.na(to),
     t_min == 0,
     t_max == Inf,
     t_sot == 0
  ))))

  # make sure "now" parameter works as intended
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,    ~t,   ~state,
            1,           1,     0, "stable"
  )
  tbl_mstate <- visits_to_mstate(
    tbl_visits,
    start_state = "stable",
    absorbing_states = c("progression"),
    now = 1
  )
  expect_true(tbl_mstate$t_min == 1)

  # check longer trajectory
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,    ~t,        ~state,
            1,           1,     0,      "stable",
            1,           1,   1.2,      "stable",
            1,           1,   2.4,    "response",
            1,           1,   3.6, "progression"
  )
  tbl_mstate <- visits_to_mstate(
    tbl_visits,
    start_state = "stable",
    absorbing_states = c("progression"),
    now = max(tbl_visits$t) + 1
  )
  expect_true(all(tbl_mstate$from == c("stable", "response")))
  expect_true(all(tbl_mstate$to == c("response", "progression")))
  expect_equal(tbl_mstate$t_min, c(1.2, 2.4))
  expect_equal(tbl_mstate$t_max, c(2.4, 3.6))

  # check EOF trajectory
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,    ~t,        ~state,
            1,           1,     0,      "stable",
            1,           1,   1.2,      "stable",
            1,           1,   2.4,    "response",
            1,           1,   3.6,         "EOF"
  )
  tbl_mstate <- visits_to_mstate(
    tbl_visits,
    start_state = "stable",
    absorbing_states = c("progression"),
    now = max(tbl_visits$t) + 1
  )
  expect_true(all(tbl_mstate$from == c("stable", "response")))
  expect_true(tbl_mstate$to[1] == "response")
  expect_true(is.na(tbl_mstate$to[2]))
  expect_equal(tbl_mstate$t_min, c(1.2, 3.6))
  expect_equal(tbl_mstate$t_max, c(2.4, -Inf))

  # check multiple individuals
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
  expect_true(all(tbl_mstate$from == c("stable", "response", "stable",
                                       "response", "stable", "response")))
  expect_equal(tbl_mstate$t_min, c(1.2, 2.4, 2, 3.6, 1.5, 3.5))
  expect_equal(tbl_mstate$t_max, c(2.4, 3.6, 3.0, Inf, 3.0, -Inf))
  expect_equal(tbl_mstate$t_sot, c(0, 0, 2, 2, 1.5, 1.5))

  # throw error
  tbl_visits <- tibble::tribble(
      ~group_id, ~subject_id,    ~t,        ~state,
              1,           1,     0,      "stable",
              1,           1,   1.2,      "stable",
              1,           1,   2.4,    "response",
              1,           1,   3.6, "progression",
              1,           2,     3,    "response"
    )
  expect_error(
    tbl_mstate <- visits_to_mstate(
      tbl_visits,
      start_state = "stable",
      absorbing_states = c("progression")
    ),
    regexp = "first visit must be in starting state; subject_id=2, state=response" # nolint
  )

})
