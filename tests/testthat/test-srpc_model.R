test_that("can fit SRPC model", {
  library(tidyverse)
  tbl_data <- tbl_example_data %>%
    mutate(
      group = sample(1:2, size = n(), replace = TRUE)
    )

  tbl_data

  fit <- rstan::sampling(
    bhmbasket.predict:::stanmodels$srpc_model,
    data = list(
      M_groups = 2,
      N = nrow(tbl_data),
      group_id = as.array(tbl_data$group),
      from = as.array(tbl_data$from),
      to = as.array(tbl_data$to),
      tstart = as.array(tbl_data$tmin + sqrt(.Machine$double.eps)),
      tstop = as.array(tbl_data$tmax + 2* sqrt(.Machine$double.eps)),

      logodds_mean = matrix(c(
          logodds(.2), logodds(.5),
          logodds(.2), logodds(.5)
        ), byrow = TRUE, nrow = 2, ncol = 2),
      logodds_sd = matrix(1, nrow = 2, ncol = 2),
      logodds_min = matrix(logodds(.001), nrow = 2, ncol = 2),
      logodds_max = matrix(logodds(.999), nrow = 2, ncol = 2),
      shape = matrix(.5, nrow = 2, ncol = 3),
      median_time_to_event_mean = matrix(c(
          100, 200, 200,
          100, 200, 200
        ), byrow = TRUE, nrow = 2, ncol = 3),
      median_time_to_event_sd = matrix(300, nrow = 2, ncol = 3)
    ),
    chains = 1L, iter = 2000L, warmup = 1000L
  )

  dim(rstan::extract(fit, "p")$p)
  rstan::extract(fit, "p")$p[1, , ]

})
