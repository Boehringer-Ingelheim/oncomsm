test_that("decision rules can be simulated", {
  mdl <- create_srpmodel(
    A = define_srp_prior(recruitment_rate = 2)
  )
  # create some interim data with only successes
  tbl_interim <- sample_predictive(mdl, 10, nsim = 1, seed = 45968, p = 1) %>%
    filter(t <= 10)
  # presample interim posterior
  interim_post <- sample_posterior(mdl, tbl_interim, nsim = 500L, seed = 243756)
  # define decision rule, simply check whether posterior of response probability
  # being above a threshold is larger than confidence
  rule <- function(model, data, threshold = 0.7, confidence = 0.75, nsim = 500L) {
    n_groups <- length(model$group_id)
    smpl <- sample_posterior(model, data = data, warmup = 250L, nsim = nsim)
    p_post <- parameter_sample_to_tibble(mdl, smpl) %>%
        filter(parameter == "p") %>%
        transmute(
          iter = 1:(n_groups * nsim),
          group_id,
          p = value
        )
    res <- p_post %>%
      group_by(group_id) %>%
      summarize(
        go = mean(p_post >= threshold) >= confidence
      )
    return(res)
  }
  # test posterior predictive
  tbl_decisions <- simulate_decision_rule(
    mdl, 40L, rule, data = tbl_interim, parameter_sample = interim_post,
    seed = 234, nsim = 10L
  )
  expect_true(
    mean(tbl_decisions$go) > 0.5
  )
  # same without using pre sampled
  tbl_decisions <- simulate_decision_rule(
    mdl, 40L, rule, data = tbl_interim, seed = 56750789, nsim = 10L
  )
  expect_true(
    mean(tbl_decisions$go) > 0.5
  )
  # same without any interim data
  tbl_decisions <- simulate_decision_rule(
    mdl, 40L, rule, seed = 5675, nsim = 10L
  )
  expect_true(
    mean(tbl_decisions$go) < 0.5
  )
})

