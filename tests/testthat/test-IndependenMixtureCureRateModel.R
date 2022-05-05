test_that("can create IndependentMixtureCureRateModel for single arm", {

  set.seed(125L)

  mdl <- model(
    tte_model = independent_mixture_cure_rate_model(
      group_id = "A",
      logodds_mean = log(.3/.7),
      logodds_sd = 1,
      shape_mean = 5,
      shape_sd = 1,
      scale_mean = 3,
      scale_sd = 2,
      visit_spacing = 1.4
    ),
    recruitement_model = independent_poisson_recruitment_model(
      group_id = "A",
      log_monthly_rate_mean = log(2),
      log_monthly_rate_sd = .5,
      maximal_recruitment_interval = 2
    )
  )

  # plot(mdl, n = 30)

})
