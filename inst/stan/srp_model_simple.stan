// stable/response/progression(or death) model

functions {

  real beta_mix_trunc_lpdf(real p, real mean, real n, real eta, real lower_inp, real upper_inp) {
    real log_eps = log(1e-6);
    real a = mean * n;
    real b = (1 - mean) * n;
    real log1m_eta = log1m(eta);
    real log_eta = log(eta);
    // normalizing factor (CDF upper - lower)
    real norm = log_sum_exp(log1m_eta + log_diff_exp(beta_lcdf(upper_inp | a, b), beta_lcdf(lower_inp | a, b)),
                            log_eta + log(upper_inp - lower_inp));
    real lpdf = log_sum_exp(log1m_eta + beta_lpdf(p | a, b), log_eta) - norm;
    if (p < lower_inp || p > upper_inp) {
      reject(p);
    }
    return log_sum_exp(lpdf, log_eps);
  }

}



data {

  int<lower=1> M_groups;
  int<lower=0> N;
  int<lower=0> N_subjects;
  real<lower=0> maximal_time;

  array[N] int<lower=1> subject_id;
  array[N] int<lower=1> group_id;
  array[N] int<lower=1,upper=2> from;
  array[N] int<lower=2,upper=4> to; // 4 := censored

  array[N] real<lower=0> t_min;
  array[N] real<lower=0> t_max;

  // event probability
  array[M_groups] real<lower=machine_precision(),upper=1-machine_precision()> p_mean;
  array[M_groups] real<lower=machine_precision()> p_n;
  array[M_groups] real<lower=0.0,upper=1.0> p_eta;
  array[M_groups] real<lower=0.0,upper=1.0> p_min;
  array[M_groups] real<lower=0.0,upper=1.0> p_max;
  // weibull shape
  array[M_groups, 3] real shape_mu;
  array[M_groups, 3] real<lower=machine_precision()> shape_sigma;
  // weibull median
  array[M_groups, 3] real median_t_mu;
  array[M_groups, 3] real<lower=machine_precision()> median_t_sigma;

}



parameters {

  array[M_groups, 3] real<lower=0.5> shape;
  array[M_groups] real<lower=0,upper=1> p_raw; // the boundaries here are for scaling!
  array[M_groups, 3] real<lower=sqrt(machine_precision())> median_t;

}



transformed parameters {

  array[M_groups] real p;
  array[M_groups, 3] real scale;

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
  real log_eps = log(eps);
  int g;
  int s;
  // instead of sampling the actual jumping times, we just approximate them as
  // midpoints
  array[N_subjects] real t_jump_from_stable;

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
        target += log_sum_exp(
          log(p[g]) + log_diff_exp(weibull_lcdf(t_max[i] + eps | shape[g, 1], scale[g, 1]),
                                    weibull_lcdf(t_min[i] + eps | shape[g, 1], scale[g, 1])),
          log_eps
        );
      }
      if (to[i] == 3) { // stable -> progression
        target += log_sum_exp(
          log1m(p[g]) + log_diff_exp(weibull_lcdf(t_max[i]+ eps | shape[g, 2], scale[g, 2]),
                                      weibull_lcdf(t_min[i]+ eps | shape[g, 2], scale[g, 2])),
          log_eps
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at risk or right censored)
        // here we need to use the mixture distribution since the next state
        // is unknown
        target += log_sum_exp(
          log_sum_exp(
            log(p[g]) + log1m_exp(weibull_lcdf(t_min[i]+ eps | shape[g, 1], scale[g, 1])),
            log1m(p[g]) + log1m_exp(weibull_lcdf(t_min[i]+ eps | shape[g, 2], scale[g, 2]))
          ),
          log_eps
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
        target += log_sum_exp(
          log_diff_exp(weibull_lcdf(t_max[i] - t_jump_from_stable[s] + eps | shape[g, 3], scale[g, 3]),
                        weibull_lcdf(t_min[i] - t_jump_from_stable[s] + eps | shape[g, 3], scale[g, 3])),
          log_eps
        );
      }
      if (to[i] == 4) { // stable -> ??? (still at risk, right censored)
        target += log_sum_exp(
          log1m_exp(weibull_lcdf(t_min[i] - t_jump_from_stable[s] + eps | shape[g, 3], scale[g, 3])),
           log_eps
        );
      }
    } // end from[i] == "response"
  }
}
