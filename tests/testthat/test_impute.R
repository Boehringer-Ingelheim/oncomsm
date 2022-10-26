test_that("Testing consistency of .impute function to return same values", {

  mdl <<- create_srp_model(
    group_id = 1:3,
    logodds_mean =  c(logodds(.5), logodds(.5), logodds(.5)),
    logodds_sd = c(.5, .5, .5),
    visit_spacing = c(1.2, 1.2, 1.2),
    median_time_to_next_event = matrix(c(
      3, 3, 6,
      3, 3, 6,
      3, 3, 6
    ), byrow = TRUE,  nrow = 3, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 3, ncol = 3)
  )
  smpl_prior <- sample_prior(mdl, warmup = 500, nsim = 2000, seed = 36L)

  tbl1 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 1L, 0L),
                            nsim = 20,
                            seed = 3423423)
  tbl2 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 1L, 0L),
                            nsim = 20,
                            seed = 3423423)
  expect_true(all(tbl1 == tbl2))
})

test_that("Testing corner case values when response and progression occur within
          the same visit interval", {

  mdl2 <- create_srp_model(
    group_id = 1,
    logodds_mean =  c(logodds(.5)),
    logodds_sd = c(.5),
    visit_spacing = c(1.2),
    median_time_to_next_event = matrix(c(
      0.5, 0.7, 0.6
    ), byrow = TRUE,  nrow = 1, ncol = 3),
    median_time_to_next_event_sd = matrix(1, byrow = TRUE,  nrow = 1, ncol = 3)
  )
  smpl_prior2 <- sample_prior(mdl2, warmup = 500, nsim = 2000, seed = 36L)

  tbl3 <- sample_predictive(mdl2,
                            sample = smpl_prior2,
                            n_per_group = c(5L),
                            nsim = 100,
                            seed = 3423423)
  corner_check <- tbl3 %>%
         group_by(iter, subject_id) %>%
         filter(t_min == lag(t_min) & t_max == lag(t_max)) %>%
    nrow()
  expect_true(corner_check == 0)
  })
test_that("Testing if a unique indivudual can is sampled only once in a group",
  {
    smpl_prior <- sample_prior(mdl, warmup = 500, nsim = 2000, seed = 36L)
    tbl4 <- sample_predictive(mdl,
                            sample = smpl_prior,
                            n_per_group = c(1L, 2L, 1L),
                            nsim = 2,
                            seed = 3423423,
                            debug = TRUE)
    tbl_data <- tbl4 %>%
    filter(group_id == 1) %>%
    count(subject_id)
    expect_true(tbl_data["subject_id"] == 1)
  })
