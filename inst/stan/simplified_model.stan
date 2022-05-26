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

  real<lower=0> dt1_A[N_A];
  real<lower=0> dt2_A[N_A];
  real<lower=0> dt1_B[N_B];

  int<lower=0> n_recruited_per_group[M_groups];
  // wating times between recruitment as square matrix
  // dt_recruitment[i,j] is the waiting time from the i-1 th to the ith subject in group j
  real<lower=0> dt_recruitment[N_A+N_B+N_C, M_groups];

  // prior hyperparamters
  real visit_spacing[M_groups];
  // response probability (via logit)
  real logodds_mean[M_groups];
  real<lower=machine_precision()> logodds_sd[M_groups];
  vector[M_groups] logodds_min;
  vector[M_groups] logodds_max;
  // alpha (weibull shape)
  real log_shape_mean[M_groups];
  real<lower=machine_precision()> log_shape_sd[M_groups];
  // weibull median
  real median_time_to_response_mean[M_groups];
  real<lower=machine_precision()> median_time_to_response_sd[M_groups];
  // poisson recruitment
  real monthly_rate_mean[M_groups];
  real<lower=machine_precision()> monthly_rate_sd[M_groups];

}


transformed data {

  int N_all = N_A + N_B + N_C + N_D;

}



parameters {

  vector<lower=0,upper=1>[M_groups] logodds_raw;
  real log_shape[M_groups]; // shape aka k, important to bound from 0
  real<lower=sqrt(machine_precision())> median_time_to_response[M_groups];
  real<lower=sqrt(machine_precision())> monthly_rate[M_groups];

}



transformed parameters {

  real p[M_groups];
  real shape[M_groups] = exp(log_shape);
  real scale[M_groups];
  vector[M_groups] logodds = logodds_min + (logodds_max - logodds_min) .* logodds_raw; // https://mc-stan.org/docs/2_18/stan-users-guide/vectors-with-varying-bounds.html

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
    log_shape[g] ~ normal(log_shape_mean[g], log_shape_sd[g]);
    median_time_to_response[g] ~ normal(median_time_to_response_mean[g], median_time_to_response_sd[g]) T[sqrt(machine_precision()),];
    monthly_rate[g] ~ normal(monthly_rate_mean[g], monthly_rate_sd[g]) T[sqrt(machine_precision()),];
  }

  // recruitment
  for (i in 1:M_groups) {
    dt_recruitment[1:n_recruited_per_group[i], i] ~ exponential(monthly_rate[i]);
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
  for (i in 1:N_B) {
    group_id = group_id_B[i];
    target += log(
      (1 - p[group_id]) + p[group_id] * (1 - weibull_cdf( dt1_B[i], shape[group_id], scale[group_id] ))
    );
  }

}
