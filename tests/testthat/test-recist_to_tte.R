test_that("convert visit data to time-to-event", {

  tbl_tte_1 <- tibble::tribble(
      ~group_id, ~subject_id, ~t, ~status, ~eof,
      "A",         "1",  0,    "S", TRUE
    ) %>%
    visits_to_tte()
  expect_true(all(
    tbl_tte_1$group_id == "A",
    tbl_tte_1$subject_id == "1",
    tbl_tte_1$t_recruitment == 0,
    tbl_tte_1$dt1 == Inf,
    tbl_tte_1$dt2 == Inf
  ))

  tbl_tte_2 <- tibble::tribble(
      ~group_id, ~subject_id, ~t, ~status, ~eof,
      "A",         "1",  0,    "S", FALSE,
      "A",         "1",  1,    "S", FALSE
    ) %>%
    visits_to_tte()
  expect_true(all(
    tbl_tte_2$group_id == "A",
    tbl_tte_2$subject_id == "1",
    tbl_tte_2$t_recruitment == 0,
    tbl_tte_2$dt1 == 1,
    is.na(tbl_tte_2$dt2)
  ))

  tbl_tte_3 <- tibble::tribble(
      ~group_id, ~subject_id, ~t, ~status, ~eof,
      "A",         "1",  0,    "S", FALSE,
      "A",         "1",  1,    "P", TRUE
    ) %>%
    visits_to_tte()
  expect_true(all(
    tbl_tte_3$group_id == "A",
    tbl_tte_3$subject_id == "1",
    tbl_tte_3$t_recruitment == 0,
    tbl_tte_3$dt1 == Inf,
    tbl_tte_3$dt2 == Inf
  ))

  tbl_tte_4 <- tibble::tribble(
      ~group_id, ~subject_id, ~t, ~status, ~eof,
      "A",         "1",  0,    "S", FALSE,
      "A",         "1",  1,    "R", FALSE
    ) %>%
    visits_to_tte()
  expect_true(all(
    tbl_tte_4$group_id == "A",
    tbl_tte_4$subject_id == "1",
    tbl_tte_4$t_recruitment == 0,
    tbl_tte_4$dt1 == 0,
    tbl_tte_4$dt2 == 1
  ))

})

test_that("convert visit data to time-to-event preserves individuals", {

  tbl_tte <- tibble::tribble(
      ~group_id, ~subject_id, ~t, ~status, ~eof,
            "A",          "1", 2,     "S", TRUE
    ) %>%
    visits_to_tte(1)
  expect_true(all(
    tbl_tte$group_id == "A",
    tbl_tte$subject_id == "1",
    is.na(tbl_tte$t_recruitment),
    is.na(tbl_tte$dt1),
    is.na(tbl_tte$dt2)
  ))

})
