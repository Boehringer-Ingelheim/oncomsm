functions {

  real ttweibull_rng(real alpha, real sigma, real t1, real t2) {
    real p1;
    real p2;
    real u;
    real x;
    p1 = weibull_cdf(t1, alpha, sigma);
    p2 = weibull_cdf(t2, alpha, sigma);
    if (p2 - p1 < 1e-4){
      u = p1;
    } else {
      u = uniform_rng(p1, p2);
    }
    u = fmax(u, 1e-4);
    u = fmin(u, 1 - 1e-4);
    x = sigma * (-log1m(u))^(1/alpha);
    return x;
  }

}


// A = interval censored
// B = right censored
// C = non-responder
// D = to be recruited

data {

  int<lower=1> M_groups;

  int<lower=0> N_A;
  int<lower=0> N_B;
  int<lower=0> N_C;
  int<lower=0> N_D;

  int<lower=1> group_id_A[N_A];
  int<lower=1> group_id_B[N_B];
  int<lower=1> group_id_C[N_C];
  int<lower=1> group_id_D[N_D];

  int<lower=1> subject_id_A[N_A];
  int<lower=1> subject_id_B[N_B];
  int<lower=1> subject_id_C[N_C];
  int<lower=1> subject_id_D[N_D];

  real<lower=machine_precision()> dt1_A[N_A];
  real<lower=machine_precision()> dt2_A[N_A];
  real<lower=machine_precision()> dt1_B[N_B];

  // prior hyperparamters
  real visit_spacing[M_groups];
  // response probability (via logit)
  real logodds_mean[M_groups];
  real<lower=machine_precision()> logodds_sd[M_groups];
  real logodds_min[M_groups];
  real logodds_max[M_groups];
  // alpha (weibull shape)
  real shape_mean[M_groups];
  real<lower=machine_precision()> shape_sd[M_groups];
  // weibull median
  real median_time_to_response_mean[M_groups];
  real<lower=machine_precision()> median_time_to_response_sd[M_groups];
  // response time truncation
  real max_time_to_response[M_groups];

}


transformed data {

  int N_all = N_A + N_B + N_C + N_D;

}



parameters {

  real logodds[M_groups];
  real<lower=1-machine_precision(),upper=99> shape[M_groups]; // shape aka k, important to bound from 0
  real<lower=1.0/30.0,upper=99> median_time_to_response[M_groups];

}



transformed parameters {

  real<lower=0,upper=1> p[M_groups];
  real<lower=machine_precision()> scale[M_groups];

  for (g in 1:M_groups) {
    p[g] = 1/(1 + exp(-logodds[g]));
    scale[g] = median_time_to_response[g]/(log(2)^(1/shape[g]));
  }

}



model {

  int group_id;

  // prior
  for (g in 1:M_groups) {
    logodds[g] ~ normal(logodds_mean[g], logodds_sd[g]) T[logodds_min[g],logodds_max[g]];
    shape[g] ~ normal(shape_mean[g], shape_sd[g]) T[1 - machine_precision(),99];
    median_time_to_response[g] ~ normal(median_time_to_response_mean[g], median_time_to_response_sd[g]) T[machine_precision(),99];
  }

  // likelihood for the definite non-responders
  for (i in 1:N_C) {
    group_id = group_id_C[i];
    target += log(1 - p[group_id]);
  }

  // likelihood for the interval censored definite responders
  for (i in 1:N_A) {
    group_id = group_id_A[i];
    target += log(
      p[group_id] * (weibull_cdf(dt2_A[i], shape[group_id], scale[group_id]) -
        weibull_cdf(dt1_A[i], shape[group_id], scale[group_id])
      )
    );
  }

  // likelihood for the right censored individuals (still at risk)
  for (i in 1:N_B) { // vectorize?
    group_id = group_id_B[i];
    target += log(
      (1 - p[group_id]) + p[group_id] * (1 - weibull_cdf( dt1_B[i], shape[group_id], scale[group_id] ))
    );
  }

}



generated quantities {

  int<lower=1> group_id[N_all];
  int<lower=1> subject_id[N_all];
  int<lower=1> gg;
  real<lower=0> dt[N_all];
  real<lower=0> dt1[N_all];
  real<lower=0> dt2[N_all];

  real p_cond = 0.0; // buffer variable to make code more legible
  real S_t = 0.0; // buffer variable to make code more legible
  int offset = 0;
  int idx;

  // definite responder
  for (i in 1:N_A) { // definite responder, directly sample from truncated weibull
    gg = group_id_A[i];
    group_id[i] = group_id_A[i];
    subject_id[i] = subject_id_A[i];
    dt[i] = ttweibull_rng(shape[gg], scale[gg], dt1_A[i], dt2_A[i]);
    dt1[i] = dt1_A[i];
    dt2[i] = dt2_A[i];
  }
  offset = N_A;

  // right censored
  for (i in 1:N_B) {
    gg = group_id_B[i];
    idx = i + offset;
    group_id[idx] = group_id_B[i];
    subject_id[idx] = subject_id_B[i];
    // Pr[ responder | no response up to t ] is not constant, need to use Bayes Theorem
    S_t = 1 - weibull_cdf(dt1_B[i], shape[gg], scale[gg]);
    p_cond = S_t*p[gg] / ( (1 - p[gg]) + p[gg]*S_t );
    if (bernoulli_rng(p_cond) == 1) { // responder
      // upper truncation is just for numerical stability
      dt[idx] = ttweibull_rng(shape[gg], scale[gg], dt1_B[i], max_time_to_response[gg]);
      // assuming fixed visit intervals, increment lower boundary until next visit is after dt
      dt1[idx] = 0;
      while (dt1[idx] + visit_spacing[gg] < dt[idx]) {
        dt1[idx] += visit_spacing[gg];
      };
      dt2[idx] = dt1[idx] + visit_spacing[gg]; // corresponding upper bound
    } else { // non responder, set to infinity
      dt[idx] = positive_infinity();
      dt1[idx] = positive_infinity();
      dt2[idx] = positive_infinity();
    }
  }
  offset += N_B;

  // definite non-responder
  for (i in 1:N_C) {
    idx = i + offset;
    group_id[idx] = group_id_C[i];
    subject_id[idx] = subject_id_C[i];
    dt[idx] = positive_infinity();
    dt1[idx] = positive_infinity();
    dt2[idx] = positive_infinity();
  }
  offset += N_C;


  // new individuals
  for (i in 1:N_D) {
    gg = group_id_D[i];
    idx = i + offset;
    group_id[idx] = group_id_D[i];
    subject_id[idx] = subject_id_D[i];
    if (bernoulli_rng(p[gg]) == 1) { // responder
      // truncation is just for numerical stability
      dt[idx] = ttweibull_rng(shape[gg], scale[gg], 1.0/30.0, max_time_to_response[gg]);
      // assuming fixed visit intervals, increment lower boundary until next visit is after dt
      dt1[idx] = 0;
      while (dt1[idx] + visit_spacing[gg] < dt[idx]) {
        dt1[idx] += visit_spacing[gg];
      };
      dt2[idx] = dt1[idx] + visit_spacing[gg]; // corresponding upper bound
    } else { // non-responder
      dt[idx] = positive_infinity();
      dt1[idx] = positive_infinity();
      dt2[idx] = positive_infinity();
    }
  }

}

