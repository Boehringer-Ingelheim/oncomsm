// stable/response/progression(or death) model

functions {

  real weibull_trunc_lpdf (real t, real shape, real scale, real a, real b) {
    real eps = 1e-6; // numerical stability, minimum interval width
    real delta_cdf = weibull_cdf(b, shape, scale) - weibull_cdf(a, shape, scale);
    return weibull_lpdf(t | shape, scale) - log(delta_cdf + eps);
  }

}



data {

  int<lower=1> M_groups;
  int<lower=0> N;
  int<lower=0> N_subjects;

  int<lower=1> subject_id[N];
  int<lower=1> group_id[N];
  int<lower=1,upper=2> from[N];
  int<lower=2,upper=4> to[N]; // 4 := censored

  real<lower=0> t_min[N];
  real<lower=0> t_max[N];

  // event probability (via logit)
  // logodds[i] is logodds of group i for response from stable
  real logodds_mean[M_groups];
  real<lower=machine_precision()> logodds_sd[M_groups];
  real logodds_min[M_groups];
  real logodds_max[M_groups];
  // alpha (weibull shape)
  real<lower=machine_precision()> shape_min[M_groups, 3];
  real<lower=machine_precision()> shape_max[M_groups, 3];
  real shape_mu[M_groups, 3];
  real<lower=machine_precision()> shape_sigma[M_groups, 3];
  // weibull median
  real median_time_to_next_event_mean[M_groups, 3];
  real<lower=machine_precision()> median_time_to_next_event_sd[M_groups, 3];
  real median_t_mu[M_groups, 3];
  real<lower=machine_precision()> median_t_sigma[M_groups, 3];
  real<lower=machine_precision()> median_t_max[M_groups, 3];

}



parameters {

  real<lower=0,upper=1> logodds_raw[M_groups]; // the boundaries here are for scaling!
  real<lower=0,upper=1> t_jump_from_stable_raw[N_subjects]; // the boundaries here are for scaling!
  real<lower=sqrt(machine_precision())> median_time_to_next_event[M_groups, 3];
  real<lower=0,upper=1> shape_raw[M_groups, 3];

}



transformed parameters {

  real logodds[M_groups];
  real p[M_groups];
  real scale[M_groups, 3];
  real shape[M_groups, 3];
  real<lower=0> t_jump_from_stable[N_subjects];

  for (g in 1:M_groups) {
    logodds[g] = logodds_min[g] + (logodds_max[g] - logodds_min[g]) * logodds_raw[g]; // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html
    p[g] = 1/(1 + exp(-logodds[g]));
    for (j in 1:3) {
      shape[g, j] = shape_min[g, j] + (shape_max[g, j] - shape_min[g, j]) * shape_raw[g, j];
      scale[g, j] = median_time_to_next_event[g, j]/(log(2)^(1/shape[g, j])); // solve for scale given shape and median
    }
  }

  for (i in 1:N) {
    if (from[i] == 1) {
      if (to[i] != 4) {
        t_jump_from_stable[subject_id[i]] = t_min[i] + (t_max[i] - t_min[i]) * t_jump_from_stable_raw[subject_id[i]]; // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html
      } else {
        t_jump_from_stable[subject_id[i]] = t_min[i] + (9999 - t_min[i]) * t_jump_from_stable_raw[subject_id[i]]; // upper limit is infinity
      }
    }
  }

}



model {

  // buffer for current group_id and subject_id in loop
  int g;
  int s;
  // buffer for jump times in loop, dt since it is reset at previous junmp time
  // hence a difference (clock reset at state switching)
  real dt_jump = 0;
  real dt_jump_min = 0;
  real dt_jump_max = 0;
  // small constant used for numerical stability
  real eps = 1e-6;

  // handle prior contribution to log likelihood target
  for (gg in 1:M_groups) {
    logodds[gg] ~ normal(logodds_mean[gg], logodds_sd[gg]) T[logodds_min[gg],logodds_max[gg]];
    for (j in 1:3) {
      shape[gg, j] ~ lognormal(
        shape_mu[gg, j], shape_sigma[gg, j]
      ) T[shape_min[gg, j], shape_max[gg, j]];
      // this is a linear transformation and does not require a jacobian
      median_time_to_next_event[gg, j] ~ lognormal(
          median_t_mu[gg, j], median_t_sigma[gg, j]
        ) T[0, median_t_max[gg, j]];
    }
  }

  // handle data contribution to log likelihood target
  for (i in 1:N) {
    s = subject_id[i]; // convenience only
    g = group_id[i]; // convenience only
    if (from[i] == 1) {
      // Weibull can become unstable around 0, hence we shift everything away
      // from 0 by a tiny ammount; if/else would hinder sampler, this is smoother
      // and has minimal effect on the log likelihood
      dt_jump = t_jump_from_stable[s] + eps;
      dt_jump_min = t_min[i] + eps; // jump needs to have occured between this lower and ...
      dt_jump_max = t_max[i] + eps; // .. this upper bound
      // pick jump-specific weibull parameters
      if (to[i] == 2) { // stable -> response
        target += log(p[g] + eps) + // prevent log(0)
          weibull_trunc_lpdf(dt_jump | shape[g, 1], scale[g, 1], dt_jump_min, dt_jump_max);
      }
      if (to[i] == 3) { // stable -> progression
        target += log(1 - p[g] + eps) + // prevent log(0)
          weibull_trunc_lpdf(dt_jump | shape[g, 2], scale[g, 2], dt_jump_min, dt_jump_max);
      }
      if (to[i] == 4) { // stable -> ??? (still at risk or right censored)
        // here we need to use the mixture distribution since the next state
        // is unknown
        target += log(
                 p[g]  * exp(weibull_trunc_lpdf(dt_jump | shape[g, 1], scale[g, 1], dt_jump_min, 9999))
          + (1 - p[g]) * exp(weibull_trunc_lpdf(dt_jump | shape[g, 2], scale[g, 2], dt_jump_min, 9999))
          + eps // numerical stability, prevent log(0)
        );
      }
    } // end from[i] == "stable"
    if (from[i] == 2) {
      // this is a terminal jump. we can integrate out the exact jump time
      // directly to avoid sampling it since it is not needed to calculate
      // the likelihood contribution of a subsequent jump
      //
      // We need to substract the previous (unobserved) jump time
      // to get the boundaries for the new jump time 2->3
      // again adding small eps to avoid 0
      //
      // TODO: consider replacing fmax() with a smooth maximum for better gradients
      dt_jump_min = fmax(eps, t_min[i] - t_jump_from_stable[s] + eps);
      dt_jump_max = fmax(dt_jump_min + eps, t_max[i] - t_jump_from_stable[s] + eps);
      if (to[i] == 3) { // response -> progression
        target += log(
          weibull_cdf(dt_jump_max, shape[g, 3], scale[g, 3]) -
          weibull_cdf(dt_jump_min, shape[g, 3], scale[g, 3]) +
          eps // numerical stability
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at riks, right censored)
        target += log(
          1 - weibull_cdf(dt_jump_min, shape[g, 3], scale[g, 3]) +
          eps // numerical stability
        );
      }
    } // end from[i] == "response"
  }
}
