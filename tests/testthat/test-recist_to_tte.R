test_that("RECIST -> TTE example 1", {
  tbl_tmp <- tibble::tibble(
      subject_id = factor(c("A001")),
      group_id = factor(c("G001")),
      visit_date = as.Date(c("2042-01-01")),
      visit_status = factor("SOT", levels = c("SOT", "PD", "SD", "PR", "CR")),
      end_of_treatment = c(FALSE)
    ) %>%
    recist_to_tte()
  expect_true(tbl_tmp$subject_id == "A001")
  expect_equal(tbl_tmp$t_sot, 0)
  expect_equal(tbl_tmp$dt1, 0)
  expect_true(is.na(tbl_tmp$dt2))
})

test_that("RECIST -> TTE example 2", {
  tbl_tmp <- tibble::tibble(
      subject_id = factor(c("A001")),
      group_id = factor(c("G001")),
      visit_date = as.Date(c("2042-01-01")),
      visit_status = factor("SOT", levels = c("SOT", "PD", "SD", "PR", "CR")),
      end_of_treatment = c(TRUE)
    ) %>%
    recist_to_tte()
  expect_equal(tbl_tmp$t_sot, 0)
  expect_equal(tbl_tmp$dt1, Inf)
  expect_equal(tbl_tmp$dt2, Inf)
})

test_that("RECIST -> TTE example 3", {
  expect_error(
    tibble::tibble(
      subject_id = factor(c("A001")),
      group_id = factor(c("G001")),
      visit_date = as.Date(c("2042-01-01")),
      visit_status = factor("PD", levels = c("SOT", "PD", "SD", "PR", "CR")),
      end_of_treatment = c(FALSE)
    ) %>%
    recist_to_tte(),
    regexp = "first visit must be SOT"
  )
})

# todo: test cutoff + multiple individuals
