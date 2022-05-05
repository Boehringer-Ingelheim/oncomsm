functions {

  // truncated exponential random number generator
  // taken from https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/7
  real ttexponential_rng(real rate, real t1, real t2) {
    real p1;
    real p2;
    real u;
    real x;
    p1 = exponential_cdf(t1, rate);
    if (1 - p1 < 1e-4) {
      // unstable, simply use lower boundary
      return t1;
    }
    p2 = exponential_cdf(t2, rate);
    if (p2 - p1 < 1e-4) {
      // unstable, simply use midpoint
      return (t1 + t2)/2;
    }
    u = uniform_rng(p1, p2);
    x = -log(1 - u)/rate;
    return x;
  }

}



data {

  int<lower=1> M_groups;

  int<lower=0> N_old;
  int<lower=1> group_id_old[N_old];
  int<lower=1> subject_id_old[N_old];
  real<lower=0> t_old[N_old];
  int<lower=0> N_new;
  int<lower=1> group_id_new[N_new];
  int<lower=1> subject_id_new[N_new];

  real<lower=0> now;

  // prior hyperparamters
  real log_monthly_rate_mean[M_groups];
  real<lower=machine_precision()> log_monthly_rate_sd[M_groups];
  real<lower=machine_precision()> maximal_recruitment_interval[M_groups];

}



transformed data {

  int N_total = N_old + N_new;
  int idx_old_sorted[N_old] = sort_indices_asc(t_old);
  real<lower=0> dt_old_sorted[N_old];
  if (N_old > 0) {
    dt_old_sorted[1] = 0; // we ignore the waiting time to first recruitment
  }
  for (i in 2:N_old) {
    dt_old_sorted[i] = t_old[idx_old_sorted[i]] - t_old[idx_old_sorted[i - 1]];
  }

}



parameters {

  real log_monthly_rate[M_groups];

}



transformed parameters {

  real rate[M_groups] = exp(log_monthly_rate);

}



model {

  int idx;

  // prior
  for (g in 1:M_groups) {
    log_monthly_rate[g] ~ normal(log_monthly_rate_mean[g], log_monthly_rate_sd[g]);
  }

  // likelihood of recruitment times
  for (i in 1:N_old) {
    idx = idx_old_sorted[i]; // need to match group and sorted dts
    dt_old_sorted[i] ~ exponential(rate[group_id_old[idx]]);
  }

}



generated quantities {

  int<lower=1> group_id[N_total];
  int<lower=1> subject_id[N_total];
  real<lower=0> t_recruitment[N_total];
  int offset = 0;
  int idx;
  real max_wait;
  real t_last[M_groups];

  for (g in 1:M_groups) {
    t_last[g] = 0;
  }

  // store input data for observed individuals
  for (i in 1:N_old) {
    group_id[i] = group_id_old[i];
    subject_id[i] = subject_id_old[i];
    t_recruitment[i] = t_old[i];
    t_last[group_id[i]] = fmax(t_last[group_id[i]], t_recruitment[i]);
  }
  offset = N_old;

  // sample data for new individuals
  for (i in 1:N_new) {
    idx = i + offset;
    group_id[idx] = group_id_new[i];
    subject_id[idx] = subject_id_new[i];
    max_wait = maximal_recruitment_interval[group_id[idx]];
    t_recruitment[idx] = t_last[group_id[idx]] + ttexponential_rng(rate[group_id[idx]], fmax(0, now - t_last[group_id[idx]]), max_wait);
    t_last[group_id[idx]] = fmax(t_last[group_id[idx]], t_recruitment[idx]);
  }

}
