functions {

  // taken from https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/7
  real tweibull_rng(real alpha, real sigma, real t) {
    real p;
    real u;
    real x;
    p = weibull_cdf(t, alpha, sigma);
    if (1 - p < 1e-4) {
      // unstable, simply use lower boundary
      return t + 2*machine_precision();
    }
    u = uniform_rng(p, 1);
    x = sigma * (-log1m(u))^(1/alpha);
    return x;
  }

  real ttweibull_rng(real alpha, real sigma, real t1, real t2) {
    real p1;
    real p2;
    real u;
    real x;
    p1 = weibull_cdf(t1, alpha, sigma);
    if (1 - p1 < 1e-4) {
      // unstable, simply use lower boundary
      return t1;
    }
    p2 = weibull_cdf(t2, alpha, sigma);
    if (p2 - p1 < 1e-4) {
      // unstable, simply use midpoint
      return (t1 + t2)/2;
    }
    u = uniform_rng(p1, p2);
    x = sigma * (-log1m(u))^(1/alpha);
    return x;
  }

}


// A = interval censored
// B = right censored
// C = non-responder
// D = to be recruited

data {

  int<lower=0> N_A;
  int<lower=0> N_B;
  int<lower=0> N_C;
  int<lower=0> N_D;

  int<lower=1> subject_id_A[N_A];
  int<lower=1> subject_id_B[N_B];
  int<lower=1> subject_id_C[N_C];
  int<lower=1> subject_id_D[N_D];

  real<lower=0> dt1_A[N_A];
  real<lower=machine_precision()> dt2_A[N_A];
  real<lower=machine_precision()> dt1_B[N_B];

  real dvisit; // visit spacing

  // prior hyperparameters
  // response probability (via logit)
  real prior_logodds_loc;
  real<lower=machine_precision()> prior_logodds_scale;
  // alpha (weibull shape)
  real prior_alpha_loc;
  real<lower=machine_precision()> prior_alpha_scale;
  // sigma (weibull scale)
  real prior_sigma_loc;
  real<lower=machine_precision()> prior_sigma_scale;

}


transformed data {

  int N_all = N_A + N_B + N_C + N_D;

}



parameters {

  real theta;
  real<lower=1-machine_precision(),upper=99> alpha; // shape aka k
  real<lower=machine_precision(),upper=99> sigma; // scale aka lambda

}



transformed parameters {

  real<lower=0, upper=1> p = 1/(1 + exp(-theta));

}



model {

  // prior
  theta ~ normal(prior_logodds_loc, prior_logodds_scale);
  alpha ~ normal(prior_alpha_loc, prior_alpha_scale) T[1,99]; // the "risk" of response must be increasing
  sigma ~ normal(prior_sigma_loc, prior_sigma_scale) T[machine_precision(),99];

  // likelihood for the definite non-responders
  for (i in 1:N_C) {
    target += log(1 - p);
  }

  // likelihood for the interval censored individuals
  for (i in 1:N_A) {
    target += log(
      (1 - p) + p*(weibull_cdf(dt2_A[i], alpha, sigma) - weibull_cdf(dt1_A[i], alpha, sigma))
    );
  }

  // likelihood for the right censored individuals
  for (i in 1:N_B) { // vectorize?
    target += log(
      (1 - p) + p * (1 - weibull_cdf( dt1_B[i], alpha, sigma ))
    );
  }

}



generated quantities {

  int<lower=1> subject_id[N_all];
  real<lower=0> dt[N_all];
  real<lower=0> dt1[N_all];
  real<lower=0> dt2[N_all];

  real p_cond = 0.0; // buffer variable to make code more legible
  real S_t = 0.0; // buffer variable to make code more legible
  int offset = 0;
  int idx;

  // definite responder
  for (i in 1:N_A) { // definite responder, directly sample from truncated weibull
    subject_id[i] = subject_id_A[i];
    dt[i] = ttweibull_rng(alpha, sigma, dt1_A[i], dt2_A[i]);
    dt1[i] = dt1_A[i];
    dt2[i] = dt2_A[i];
  }
  offset = N_A;

  // right censored
  for (i in 1:N_B) {
    idx = i + offset;
    subject_id[idx] = subject_id_B[i];
    // Pr[ responder | no response up to t ] is not constant, need to use Bayes Theorem
    S_t = 1 - weibull_cdf(dt1_B[i], alpha, sigma);
    p_cond = S_t*p / ( (1-p) + p*S_t );
    if (bernoulli_rng(p_cond) == 1) { // responder
      dt[idx] = tweibull_rng(alpha, sigma, dt1_B[i]);
      // assuming fixed visit intervals, increment lower boundary until next visit is after dt
      dt1[idx] = 0;
      while (dt1[idx] + dvisit < dt[idx]) {
        dt1[idx] += dvisit;
      };
      dt2[idx] = dt1[idx] + dvisit; // corresponding upper bound
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
    subject_id[idx] = subject_id_C[i];
    dt[idx] = positive_infinity();
    dt1[idx] = positive_infinity();
    dt2[idx] = positive_infinity();
  }
  offset += N_C;


  // new individuals
  for (i in 1:N_D) {
    idx = i + offset;
    subject_id[idx] = subject_id_D[i];
    if (bernoulli_rng(p) == 1) { // responder
      dt[idx] = weibull_rng(alpha, sigma);
      // assuming fixed visit intervals, increment lower boundary until next visit is after dt
      dt1[idx] = 0;
      while (dt1[idx] + dvisit < dt[idx]) {
        dt1[idx] += dvisit;
      };
      dt2[idx] = dt1[idx] + dvisit; // corresponding upper bound
    } else { // non-responder
      dt[idx] = positive_infinity();
      dt1[idx] = positive_infinity();
      dt2[idx] = positive_infinity();
    }
  }

}

