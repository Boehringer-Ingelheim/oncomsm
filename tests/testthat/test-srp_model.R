test_that("can create SRP model", {
  mdl <- create_srp_model(
    B = srp_group_prior(),
    A = srp_group_prior()
  )
  # check class
  expect_true(
    isa(mdl, c("srp_model", "Model", "list"))
  )
  # check print method
  expect_true(
    format(mdl) == "srp_model<B,A>",
    capture.output(print(mdl)) == "srp_model<B,A> ",
    all(mdl$group_id == c("B", "A"))
  )
  # single-group model
  mdl <- create_srp_model(
    A = srp_group_prior()
  )
  # check print method
  expect_true(
    format(mdl) == "srp_model<A>",
    capture.output(print(mdl)) == "srp_model<A> ",
    all(mdl$group_id == "A")
  )
})



test_that("can calculate PFS rate", {
  mdl <- create_srp_model(
    A = srp_group_prior()
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
