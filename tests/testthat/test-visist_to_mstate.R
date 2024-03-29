test_that("visits to mstate conversion works", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # base case
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
    "1", "1", 0, "stable"
  )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl, now = 0)
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
    ~group_id, ~subject_id, ~t, ~state,
    "1", "1", 0, "stable"
  )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl, now = 1)
  expect_true(tbl_mstate$t_min == 1)

  # check longer trajectory
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
      1, 1, 0, "stable",
      1, 1, 1.2, "stable",
      1, 1, 2.4, "response",
      1, 1, 3.6, "progression"
    ) %>%
    mutate(
      group_id = as.character(group_id),
      subject_id = as.character(subject_id)
    )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl, now = max(tbl_visits$t) + 1)
  expect_true(all(tbl_mstate$from == c("stable", "response")))
  expect_true(all(tbl_mstate$to == c("response", "progression")))
  expect_equal(tbl_mstate$t_min, c(1.2, 2.4))
  expect_equal(tbl_mstate$t_max, c(2.4, 3.6))

  # check EOF trajectory
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
      1, 1, 0, "stable",
      1, 1, 1.2, "stable",
      1, 1, 2.4, "response",
      1, 1, 3.6, "EOF"
    ) %>%
    mutate(
      group_id = as.character(group_id),
      subject_id = as.character(subject_id)
    )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl, now = max(tbl_visits$t) + 1)
  expect_true(all(tbl_mstate$from == c("stable", "response")))
  expect_true(tbl_mstate$to[1] == "response")
  expect_true(is.na(tbl_mstate$to[2]))
  expect_equal(tbl_mstate$t_min, c(1.2, 3.6))
  expect_equal(tbl_mstate$t_max, c(2.4, -Inf))

  # check different labels
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C"),
    censored = "blub"
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
            1,           1, 0.0,    "A",
            1,           1, 1.2,    "A",
            1,           1, 2.4,    "B",
            1,           1, 3.6, "blub"
    ) %>%
    mutate(
      group_id = as.character(group_id),
      subject_id = as.character(subject_id)
    )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl, now = max(tbl_visits$t) + 1)
  expect_true(all(tbl_mstate$from == c("A", "B")))
  expect_true(tbl_mstate$to[1] == "B")
  expect_equal(tbl_mstate$t_max, c(2.4, -Inf))

  # check multiple individuals
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
      1, 1, 0, "stable",
      1, 1, 1.2, "stable",
      1, 1, 2.4, "response",
      1, 1, 3.6, "progression",
      1, 2, 2, "stable",
      1, 2, 3, "response",
      1, 3, 1.5, "stable",
      1, 3, 3, "response",
      1, 3, 3.5, "EOF"
    ) %>%
    mutate(
      group_id = as.character(group_id),
      subject_id = as.character(subject_id)
    )
  tbl_mstate <- visits_to_mstate(tbl_visits, mdl)
  expect_true(all(tbl_mstate$from == c(
    "stable", "response", "stable",
    "response", "stable", "response"
  )))
  expect_equal(tbl_mstate$t_min, c(1.2, 2.4, 2, 3.6, 1.5, 3.5))
  expect_equal(tbl_mstate$t_max, c(2.4, 3.6, 3.0, Inf, 3.0, -Inf))
  expect_equal(tbl_mstate$t_sot, c(0, 0, 2, 2, 1.5, 1.5))

  # throw error
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id, ~t, ~state,
      1, 1, 0, "stable",
      1, 1, 1.2, "stable",
      1, 1, 2.4, "response",
      1, 1, 3.6, "progression",
      1, 2, 3, "response"
    ) %>%
    mutate(
      group_id = as.character(group_id),
      subject_id = as.character(subject_id)
    )
  expect_error(
    tbl_mstate <- visits_to_mstate(tbl_visits, mdl),
    regexp = "first visit must be in starting state; subject_id=2, state=response" # nolint
  )
})
