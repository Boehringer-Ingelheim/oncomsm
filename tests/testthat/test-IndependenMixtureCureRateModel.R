# test_that("can create IndependentMixtureCureRateModel for single arm", {
#
#   set.seed(125L)
#
#   mdl <- model(
#     tte_model = independent_mixture_cure_rate_model(
#       group_id = "A",
#       logodds_mean = log(.3/.7),
#       logodds_sd = 1,
#       shape_mean = 5,
#       shape_sd = 1,
#       median_time_to_response_mean = 3,
#       median_time_to_response_sd = 2,
#       max_time_to_response = 10,
#       visit_spacing = 1.4
#     ),
#     recruitement_model = independent_poisson_recruitment_model(
#       group_id = "A",
#       log_monthly_rate_mean = log(2),
#       log_monthly_rate_sd = .5,
#       maximal_recruitment_interval = 2
#     )
#   )
#
#   tbl_data <- generate_visit_data(
#     group_id = "A",
#     n = 40,
#     response_rate = .33,
#     recruitment_rate = 2,
#     visit_spacing = 1.2,
#     max_duration = 48,
#     responder_weibull_scale = 4,
#     responder_weibull_shape = 5,
#     nonresponder_weibull_scale = 4,
#     nonresponder_weibull_shape = 5
#   )
#
#   plot(mdl, n = 30, nsim = 1e3)
#
#   tmp <- draw_samples(mdl$tte_model, n = 30, nsim = 1e3)
#
#   tmp <- draw_samples(mdl$tte_model, n = 30, nsim = 1e3, return_raw_stan_output = TRUE)
#
#   tmpp <- rstan::extract(tmp, c("median_time_to_response"))
#
#   summary(tmpp$median_time_to_response)
#
#   tmp2 <- draw_samples(
#     mdl,
#     tibble::tibble(
#       group_id = "A",
#       subject_id = "1",
#       t_recruitment = NA_real_,
#       dt1 = NA_real_,
#       dt2 = NA_real_
#     )
#   )
#
# })
#
#
#
# test_that("can create IndependentMixtureCureRateModel for single arm", {
#
#   set.seed(125L)
#
#   mdl <- model(
#     tte_model = independent_mixture_cure_rate_model(
#       group_id = c("A", "B", "C"),# c("C", "A", "B"),
#       logodds_mean = rep(log(.5/.5), 3),
#       logodds_sd = rep(1, 3),
#       shape_mean = rep(5, 3),
#       shape_sd = rep(1, 3),
#       median_time_to_response_mean = rep(3, 3),
#       median_time_to_response_sd = rep(2, 3),
#       max_time_to_response = rep(10, 3),
#       visit_spacing = rep(1.4, 3)
#     ),
#     recruitement_model = independent_poisson_recruitment_model(
#       group_id = c("A", "B", "C"),
#       log_monthly_rate_mean = rep(log(2), 3),
#       log_monthly_rate_sd = rep(.5, 3),
#       maximal_recruitment_interval = rep(2, 3)
#     )
#   )
#
#   tbl_data <- generate_visit_data(
#     #group_id = c("B", "A", "C"), delay = c(15, 00, 30),
#     group_id = c("A", "B", "C"), delay = c(00, 15, 30),
#     n = rep(40, 3),
#     response_rate = c(.5, .0, 1),
#     recruitment_rate = rep(2, 3),
#     visit_spacing = rep(1.2, 3),
#     max_duration = 100,
#     responder_weibull_scale = rep(4, 3),
#     responder_weibull_shape = rep(5, 3),
#     nonresponder_weibull_scale = rep(4, 3),
#     nonresponder_weibull_shape = rep(5, 3)
#   )
#
#   tmp <- draw_samples(mdl$tte_model, data = visits_to_tte(tbl_data), nsim = 1e3, return_raw_stan_output = TRUE)
#
#   plot(mdl, n = 30, nsim = 1e3)
#
#   tmp <- draw_samples(mdl$tte_model, n = 30, nsim = 1e3)
#
#   tmp <- draw_samples(mdl$tte_model, n = 30, nsim = 1e3, return_raw_stan_output = TRUE)
#
#   tmpp <- rstan::extract(tmp, c("median_time_to_response"))
#
#   summary(tmpp$median_time_to_response)
#
#   tmp2 <- draw_samples(
#     mdl,
#     tibble::tibble(
#       group_id = "A",
#       subject_id = "1",
#       t_recruitment = NA_real_,
#       dt1 = NA_real_,
#       dt2 = NA_real_
#     )
#   )
#
# })
