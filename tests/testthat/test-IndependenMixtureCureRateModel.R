test_that("can create IndependentMixtureCureRateModel for single arm", {

  set.seed(125L)

  mdl <- model(
    tte_model = independent_mixture_cure_rate_model(
      group_id = "A",
      logodds_mean = log(.3/.7),
      logodds_sd = 1,
      shape_mean = 5,
      shape_sd = 1,
      median_time_to_response_mean = 3,
      median_time_to_response_sd = 2,
      visit_spacing = 1.4
    ),
    recruitement_model = independent_poisson_recruitment_model(
      group_id = "A",
      log_monthly_rate_mean = log(2),
      log_monthly_rate_sd = .5,
      maximal_recruitment_interval = 2
    )
  )

  tbl_data <- generate_visit_data(
    group_id = "A",
    n = 40,
    response_rate = .33,
    recruitment_rate = 2,
    visit_spacing = 1.2,
    max_duration = 48,
    responder_weibull_scale = 4,
    responder_weibull_shape = 5,
    nonresponder_weibull_scale = 4,
    nonresponder_weibull_shape = 5
  )

  draw_samples(mdl, visits_to_tte(tbl_data))

})
