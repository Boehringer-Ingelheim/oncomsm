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

  int<lower=0> N_old;
  int<lower=1> subject_id_old[N_old];
  real<lower=0> t_old[N_old];
  int<lower=0> N_new;
  int<lower=1> subject_id_new[N_new];
  real<lower=0> now;

  // theta = log(rate) exponential waiting time for recruitment
  real prior_theta_loc;
  real<lower=machine_precision()> prior_theta_scale;
  real<lower=machine_precision()> dt_max;

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

  real theta;

}



transformed parameters {

  real rate  = exp(theta);

}



model {

  // prior
  theta ~ normal(prior_theta_loc, prior_theta_scale);

  // likelihood of recruitment times
  dt_old_sorted ~ exponential(rate);

}



generated quantities {

  real<lower=0> subject_id[N_total];
  real<lower=0> t[N_total];
  int offset = 0;
  int idx;

  // store input data for observed individuals
  for (i in 1:N_old) {
    subject_id[i] = subject_id_old[i];
    t[i] = t_old[i];
  }
  offset = N_old;

  // sample data for new individuals
  for (i in 1:N_new) {
    idx = i + offset;
    subject_id[idx] = subject_id_new[i];
    if (i == 1) {
      t[idx] = now + ttexponential_rng(rate, 0, dt_max);
    } else {
      t[idx] = t[idx - 1] + ttexponential_rng(rate, 0, dt_max);
    }
  }

}
