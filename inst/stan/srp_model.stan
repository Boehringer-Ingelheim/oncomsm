// stable/response/progression(or death) model
data {

  int<lower=1> M_groups;
  int<lower=0> N;

  int<lower=1> group_id[N];
  int<lower=1,upper=2> from[N];
  int<lower=2,upper=4> to[N]; // 4 := censored

  real<lower=0> dt_min[N];
  real<lower=0> dt_max[N];

  // event probability (via logit)
  // logodds[i] is logodds of group i for response from stable
  real logodds_mean[M_groups];
  real<lower=machine_precision()> logodds_sd[M_groups];
  real logodds_min[M_groups];
  real logodds_max[M_groups];
  // alpha (weibull shape)
  real<lower=machine_precision()> shape_min[M_groups, 3];
  real<lower=machine_precision()> shape_max[M_groups, 3];
  // weibull median
  real median_time_to_next_event_mean[M_groups, 3];
  real<lower=machine_precision()> median_time_to_next_event_sd[M_groups, 3];

}



parameters {

  real<lower=0,upper=1> logodds_raw[M_groups]; // the boundaries here are for scaling!
  real<lower=sqrt(machine_precision())> median_time_to_next_event[M_groups, 3];
  real<lower=0,upper=1> shape_raw[M_groups, 3];

}



transformed parameters {

  real logodds[M_groups];
  real p[M_groups];
  real scale[M_groups, 3];
  real shape[M_groups, 3];

  for (g in 1:M_groups) {
    logodds[g] = logodds_min[g] + (logodds_max[g] - logodds_min[g]) * logodds_raw[g]; // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html
    p[g] = 1/(1 + exp(-logodds[g]));
    for (j in 1:3) {
      shape[g, j] = shape_min[g, j] + (shape_max[g, j] - shape_min[g, j]) * shape_raw[g, j];
      scale[g, j] = median_time_to_next_event[g, j]/(log(2)^(1/shape[g, j])); // solve for scale given shape and median
    }
  }

}



model {

  int group; // buffer for current group
  real eps = 1e-6; // small value used for numerical stability of log(), see below

  // prior
  for (g in 1:M_groups) {
    logodds[g] ~ normal(logodds_mean[g], logodds_sd[g]) T[logodds_min[g],logodds_max[g]];
    for (j in 1:3) {
      shape[g, j] ~ uniform(shape_min[g, j], shape_max[g, j]);
      // this is a linear transformation and does not require a jacobian
      median_time_to_next_event[g, j] ~ normal(
          median_time_to_next_event_mean[g, j],
          median_time_to_next_event_sd[g, j]
        ) T[sqrt(machine_precision()),];
    }
  }

  // likelihood
  for (i in 1:N) {
    group = group_id[i];
    if (from[i] == to[i]) {
      reject("from[i] == to[i] for i = ", i);
    }
    if (from[i] == 1) {
      if (to[i] == 2) { // stable -> response
        target += log(
          p[group] * (
            weibull_cdf(dt_max[i], shape[group, 1], scale[group, 1]) -
            weibull_cdf(dt_min[i], shape[group, 1], scale[group, 1])
          ) +
          eps // numerical stability
        );
      }
      if (to[i] == 3) { // stable -> progression
        target += log(
          (1 - p[group]) * (
            weibull_cdf(dt_max[i], shape[group, 2], scale[group, 2]) -
            weibull_cdf(dt_min[i], shape[group, 2], scale[group, 2])
          ) +
          eps // numerical stability
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at riks, right censored)
        target += log( // mixture of the two options
          p[group] * ( 1 - weibull_cdf(dt_min[i], shape[group, 1], scale[group, 1]) ) +
          (1 - p[group]) * ( 1 - weibull_cdf(dt_min[i], shape[group, 2], scale[group, 2]) ) +
          eps // numerical stability
        );
      }
    } // end from[i] == 1
    if (from[i] == 2) {
      if (to[i] == 3) { // response -> progression
        target += log(
            weibull_cdf(dt_max[i], shape[group, 3], scale[group, 3]) -
            weibull_cdf(dt_min[i], shape[group, 3], scale[group, 3]) +
            eps // numerical stability
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at riks, right censored)
        target += log(
          1 - weibull_cdf(dt_min[i], shape[group, 3], scale[group, 3]) +
          eps // numerical stability
        );
      }
    } // end from[i] == 2
  }

}
