// stable/response/progression(or death) model

functions {

  real beta_mix_trunc_lpdf(real p, real mean, real n, real eta, real lower, real upper) {
    real eps = 1e-6;
    real a = mean * n;
    real b = (1 - mean) * n;
    // normalizing factor (CDF upper - lower)
    real norm = (1 - eta) * (beta_cdf(upper, a, b) - beta_cdf(lower, a, b)) +
      eta * (upper - lower);
    real pdf = ((1 - eta) * exp(beta_lpdf(p | a, b)) + eta) / norm;
    if (p < lower || p > upper) {
      reject(p);
    }
    return log(pdf + eps);
  }

}



data {

  int<lower=1> M_groups;
  int<lower=0> N;
  int<lower=0> N_subjects;
  real<lower=0> maximal_time;

  int<lower=1> subject_id[N];
  int<lower=1> group_id[N];
  int<lower=1,upper=2> from[N];
  int<lower=2,upper=4> to[N]; // 4 := censored

  real<lower=0> t_min[N];
  real<lower=0> t_max[N];

  // event probability
  real<lower=machine_precision(),upper=1-machine_precision()> p_mean[M_groups];
  real<lower=machine_precision()> p_n[M_groups];
  real<lower=0.0,upper=1.0> p_eta[M_groups];
  real<lower=0.0,upper=1.0> p_min[M_groups];
  real<lower=0.0,upper=1.0> p_max[M_groups];
  // weibull shape
  real shape_mu[M_groups, 3];
  real<lower=machine_precision()> shape_sigma[M_groups, 3];
  // weibull median
  real median_t_mu[M_groups, 3];
  real<lower=machine_precision()> median_t_sigma[M_groups, 3];

}



parameters {

  real<lower=0.5> shape[M_groups, 3];
  real<lower=0,upper=1> p_raw[M_groups]; // the boundaries here are for scaling!
  real<lower=sqrt(machine_precision())> median_t[M_groups, 3];

}



transformed parameters {

  real p[M_groups];
  real scale[M_groups, 3];

  for (g in 1:M_groups) {
    // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html
    p[g] = p_min[g] + (p_max[g] - p_min[g]) * p_raw[g];
    for (j in 1:3) {
      // solve for scale given shape and median
      scale[g, j] = median_t[g, j] / (log(2)^(1/shape[g, j]));
    }
  }

}



model {

  // buffer for current group_id and subject_id in loop
  real eps = 1e-6;
  int g;
  int s;
  // instead of sampling the actual jumping times, we just approximate them as
  // midpoints
  real t_jump_from_stable[N_subjects];

  // handle prior contribution to log likelihood target
  for (gg in 1:M_groups) {
    // use mixture prior with weight eta on uniform
    target += beta_mix_trunc_lpdf(p[gg] | p_mean[gg], p_n[gg], p_eta[gg],
      p_min[gg], p_max[gg]);
    // transition time Weibull parameter priors
    for (j in 1:3) {
      shape[gg, j] ~ lognormal(shape_mu[gg, j], shape_sigma[gg, j]);
      // this is a linear transformation and does not require a jacobian
      median_t[gg, j] ~ lognormal(median_t_mu[gg, j], median_t_sigma[gg, j]);
    }
  }

  // handle data contribution to log likelihood target
  for (i in 1:N) {
    s = subject_id[i]; // convenience only
    g = group_id[i]; // convenience only
    if (from[i] == 1) {
      // approximate jump time with midpoint instead of treating as random
      // should reduce number of paramters and avoid too much correlation
      // between paramters.
      t_jump_from_stable[s] = (t_min[i] + t_max[i]) / 2;
      if (to[i] == 2) { // stable -> response
        target += log(
          p[g] * (weibull_cdf(t_max[i] + eps, shape[g, 1], scale[g, 1]) -
            weibull_cdf(t_min[i] + eps, shape[g, 1], scale[g, 1])) + eps
        );
      }
      if (to[i] == 3) { // stable -> progression
        target += log(
          (1 - p[g]) * (weibull_cdf(t_max[i]+ eps, shape[g, 2], scale[g, 2]) -
            weibull_cdf(t_min[i]+ eps, shape[g, 2], scale[g, 2])) + eps
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at risk or right censored)
        // here we need to use the mixture distribution since the next state
        // is unknown
        target += log(
          p[g] * (1 - weibull_cdf(t_min[i]+ eps, shape[g, 1], scale[g, 1])) +
          (1 - p[g]) * (1 - weibull_cdf(t_min[i]+ eps, shape[g, 2], scale[g, 2])) + eps
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
      if (to[i] == 3) { // response -> progression
        target += log(
          weibull_cdf(t_max[i] - t_jump_from_stable[s] + eps, shape[g, 3], scale[g, 3]) -
          weibull_cdf(t_min[i] - t_jump_from_stable[s] + eps, shape[g, 3], scale[g, 3]) + eps
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at risk, right censored)
        target += log(
          1 - weibull_cdf(t_min[i] - t_jump_from_stable[s] + eps, shape[g, 3], scale[g, 3]) + eps
        );
      }
    } // end from[i] == "response"
  }
}
