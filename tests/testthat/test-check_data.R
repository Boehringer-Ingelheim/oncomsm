test_that("missing columns in data are flagged correctly", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C")
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,
          "A",         "1",
          "A",         "1"
  )
  expect_error(
    check_data(tbl_visits, mdl),
    regexp = "columns: t, state; missing from data"
  )
  tbl_visits <- tibble::tribble(
    ~group_id,
    "A",
    "A"
  )
  expect_error(
    check_data(tbl_visits, mdl),
    regexp = "columns: subject_id, t, state; missing from data"
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
          "A",         "1", 0.0,    "A",
          "A",         "1", 1.2,    "A",
          "A",         "1", 2.4,    "B",
          "A",         "1", 3.6,  "EOF"
  )
  expect_no_error(
    check_data(tbl_visits, mdl)
  )
})



test_that("missing data is flagged correctly", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C")
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t,        ~state,
          "A",         "1", 0.0, NA_character_,
          "A",         "1", 1.2,           "A"
  )
  expect_error(
    check_data(tbl_visits, mdl),
    regexp = "no missing data in visits allowed"
  )
})



test_that("wrong state labels are flagged", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C")
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
    "A",         "1", 0.0,    "A",
    "A",         "2", 1.2,    "D"
  )
  expect_error(
    check_data(tbl_visits, mdl),
    regexp = "state must be consistent with model specification"
  )
})


test_that("unsorted data are flagged", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C")
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
          "A",         "1", 0.0,    "A",
          "A",         "2", 1.2,    "A",
          "A",         "1", 2.4,    "B",
          "A",         "2", 3.6,  "EOS"
  )
  expect_error(
    check_data(tbl_visits, mdl),
    regexp = "data needs to be sorted by 'subject_id' and 't'"
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
    "A",         "1", 0.0,    "A",
    "A",         "1", 1.2,    "B",
    "A",         "1", 2.4,    "A",
    "A",         "1", 3.6,  "EOF"
  )
  expect_error(
    check_data(tbl_visits, mdl),
    regexp = "no reverse jumps in state are allowed"
  )
})



test_that("EOF after terminal", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C")
  )
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
      "A",         "1", 0.0,    "A",
      "A",         "1", 2.4,    "C",
      "A",         "1", 3.6,  "EOF"
  )
  expect_message(
    check_data(tbl_visits, mdl),
    regexp = "Removed censoring events after terminal state, make sure this is intended" # nolint
  )
})
