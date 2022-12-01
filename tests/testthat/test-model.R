test_that("can create model", {
  mdl <- create_model(
    B = group_prior(),
    A = group_prior()
  )
  # check class
  expect_true(
    isa(mdl, c("model", "list"))
  )
  # check print method
  expect_true(
    format(mdl) == "model<B,A>",
    capture.output(print(mdl)) == "model<B,A> ",
    all(mdl$group_id == c("B", "A"))
  )
  # single-group model
  mdl <- create_model(
    A = group_prior()
  )
  # check print method
  expect_true(
    format(mdl) == "model<A>",
    capture.output(print(mdl)) == "model<A> ",
    all(mdl$group_id == "A")
  )
})



test_that("can calculate PFS rate", {
  mdl <- create_model(
    A = group_prior()
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
