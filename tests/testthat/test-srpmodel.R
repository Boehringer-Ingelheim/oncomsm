test_that("can create model", {
  mdl <- create_srpmodel(
    B = define_srp_prior(),
    A = define_srp_prior()
  )
  # check class
  expect_true(
    isa(mdl, c("srpmodel", "list"))
  )
  # check print method
  expect_true(
    format(mdl) == "srpmodel<B,A>",
    capture.output(print(mdl)) == "srpmodel<B,A> ",
    all(mdl$group_id == c("B", "A"))
  )
  # single-group model
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  # check print method
  expect_true(
    format(mdl) == "srpmodel<A>",
    capture.output(print(mdl)) == "srpmodel<A> ",
    all(mdl$group_id == "A")
  )
})



test_that("can calculate PFS rate", {
  mdl <- create_srpmodel(
    A = define_srp_prior()
  )
  smpl <- sample_prior(mdl, 500L, seed = 42132L)
  res1 <- compute_pfs(mdl, t = 12, parameter_sample = smpl)
  # should be same as resampling under same seed
  res2 <- compute_pfs(mdl, t = 12, nsim = 500L, seed = 42132L)
  expect_true(
    all(res1 == res2)
  )
  # check monotonicity of average survival
  expect_true(all(
    compute_pfs(mdl, t = 0:25, parameter_sample = smpl) %>%
      group_by(t) %>%
      summarize(
        pfs = mean(pfs)
      ) %>%
      pull(pfs) %>%
      diff() < 0
  ))
})



test_that("can modify state labels", {
  mdl <- create_srpmodel(
    A = define_srp_prior(),
    states = c("A", "B", "C"),
    censored = "blub"
  )
  # can we sample predictive?
  tbl_smpl <- sample_predictive(mdl, nsim = 10L, n_per_group = 10L, seed = 42L)
  expect_true({
    all(tbl_smpl$state %in% mdl$states)
  })
  # can we sample posterior?
  tbl_visits <- tibble::tribble(
    ~group_id, ~subject_id,  ~t, ~state,
          "A",         "1", 0.0,    "A",
          "A",         "1", 1.2,    "A",
          "A",         "1", 2.4,    "B",
          "A",         "1", 3.6, "blub"
    )
  testthat::expect_no_error(
    sample_posterior(mdl, tbl_visits, nsim = 500, seed = 42L)
  )
})
