// stable/response/progression(or death)/cured model
// cured is unobserved
data {

  int<lower=1> M_groups;
  int<lower=0> N;

  int<lower=1> group_id[N];
  int<lower=1,upper=2> from[N];
  int<lower=2,upper=4> to[N]; // 4 := censored

  real<lower=0> tstart[N];
  real<lower=0> tstop[N];

  // event probability (via logit)
  // logodds[i, 1] is logodds of group i for response from stable
  // logodds[i, 2] is logodds of group i for progression from response
  real logodds_mean[M_groups, 2];
  real<lower=machine_precision()> logodds_sd[M_groups, 2];
  real logodds_min[M_groups, 2];
  real logodds_max[M_groups, 2];
  // alpha (weibull shape)
  real<lower=machine_precision()> shape[M_groups, 3];
  // weibull median
  real median_time_to_event_mean[M_groups, 3];
  real<lower=machine_precision()> median_time_to_event_sd[M_groups, 3];

}



parameters {

  real<lower=0,upper=1> logodds_raw[M_groups, 2]; // the boundaries here are for scaling!
  real<lower=sqrt(machine_precision())> median_time_to_event[M_groups, 3];

}



transformed parameters {

  real p[M_groups, 2];
  real scale[M_groups, 3];
  real logodds[M_groups, 2];

  for (g in 1:M_groups) {
    for (i in 1:2) {
      logodds[g, i] = logodds_min[g, i] + (logodds_max[g, i] - logodds_min[g, i]) * logodds_raw[g, i]; // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html
      p[g, i] = 1/(1 + exp(-logodds[g, i]));
    }
    for (j in 1:3) {
      scale[g, j] = median_time_to_event[g, j]/(log(2)^(1/shape[g, j])); // solve for scale given shape and median
    }
  }

}



model {

  int group; // buffer for current group

  // prior
  for (g in 1:M_groups) {
    for (i in 1:2) {
      logodds[g, i] ~ normal(logodds_mean[g, i], logodds_sd[g, i]) T[logodds_min[g, i],logodds_max[g, i]];
    }
    for (j in 1:3) {
      // this is a linear transformation and does not require a jacobian
      median_time_to_event[g, j] ~ normal(
          median_time_to_event_mean[g, j],
          median_time_to_event_sd[g, j]
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
          p[group, 1] * (
            weibull_cdf(tstop[i], shape[group, 1], scale[group, 1]) -
            weibull_cdf(tstart[i], shape[group, 1], scale[group, 1])
          )
        );
      }
      if (to[i] == 3) { // stable -> progression
        target += log(
          (1 - p[group, 1]) * (
            weibull_cdf(tstop[i], shape[group, 2], scale[group, 2]) -
            weibull_cdf(tstart[i], shape[group, 2], scale[group, 2])
          )
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at riks, right censored)
        target += log( // mixture of the two options
          p[group, 1] * ( 1 - weibull_cdf(tstart[i], shape[group, 1], scale[group, 1]) ) +
          (1 - p[group, 1]) * ( 1 - weibull_cdf(tstart[i], shape[group, 2], scale[group, 2]) )
        );
      }
    } // end from[i] == 1
    if (from[i] == 2) {
      if (to[i] == 3) { // response -> progression
        target += log(
          p[group, 2] * (
            weibull_cdf(tstop[i], shape[group, 3], scale[group, 3]) -
            weibull_cdf(tstart[i], shape[group, 3], scale[group, 3])
          )
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at riks, right censored)
        target += log( // mixture of the two options
          p[group, 2] * ( 1 - weibull_cdf(tstart[i], shape[group, 3], scale[group, 3]) ) +
          (1 - p[group, 2]) * 1 // cured fraction
        );
      }
    } // end from[i] == 2
  }

}
