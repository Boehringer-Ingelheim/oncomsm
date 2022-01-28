functions {

  // taken from https://discourse.mc-stan.org/t/rng-for-truncated-distributions/3122/7
  real tweibull_rng(real alpha, real sigma, real t) {
    real p;
    real u;
    real x;
    p = weibull_cdf(t, alpha, sigma);
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
    p2 = weibull_cdf(t2, alpha, sigma);
    u = uniform_rng(p1, p2);
    x = sigma * (-log1m(u))^(1/alpha);
    return x;
  }

}



data {

  int<lower=0> N_nonresponder;
  int<lower=0> N_interval_censored;
  int<lower=0> N_right_censored;
  int<lower=0> N_to_be_recruited;

  int<lower=1> id_nonresponder[N_nonresponder];
  int<lower=1> id_interval_censored[N_interval_censored];
  int<lower=1> id_right_censored[N_right_censored];
  int<lower=1> id_to_be_recruited[N_to_be_recruited];

  vector[N_right_censored] t1_right_censored;
  vector[N_interval_censored] t1_interval_censored;
  vector[N_interval_censored] t2_interval_censored;

  // prior hyperparameters
  // response probability (via logit)
  real prior_logor_loc;
  real<lower=machine_precision()> prior_logor_scale;
  // alpha (weibull shape)
  real prior_alpha_loc;
  real<lower=machine_precision()> prior_alpha_scale;
  // sigma (weibull scale)
  real prior_sigma_loc;
  real<lower=machine_precision()> prior_sigma_scale;

}



parameters {

  real logor_response;
  real<lower=1-machine_precision(),upper=99> alpha; // shape aka k
  real<lower=machine_precision(),upper=99> sigma; // scale aka lambda

}



transformed parameters {

  // logit(p) = logor_response <=> p = 1/(1 + exp(-logor_response))
  real<lower=0, upper=1> pr_response = 1/(1 + exp(-logor_response));

}



model {

  // prior
  logor_response ~ normal(prior_logor_loc, prior_logor_scale);
  alpha ~ normal(prior_alpha_loc, prior_alpha_scale) T[1,99]; // the "risk" of response must be increasing
  sigma ~ normal(prior_sigma_loc, prior_sigma_scale) T[machine_precision(),99];

  // likelihood for the definite non-responders (death)
  // directly increase the log likelihood (target) by log(1 - pr_response) for each
  // observed definite non-responder
  for (i in 1:N_nonresponder) { // vectorize
    target += log(1 - pr_response);
  }

  print("pr_response: ", pr_response);
  print("alpha: ", alpha);
  print("sigma: ", sigma);

  // likelihood for the interval censored individuals is F(t2) - F(t1) where F
  // is the CDF
  for (i in 1:N_interval_censored) { // vectorize?
    target += log(
      (1 - pr_response) + pr_response * (
        weibull_cdf( t2_interval_censored[i], alpha, sigma ) -
        weibull_cdf( t1_interval_censored[i], alpha, sigma )
      )
    );
  }

  // likelihood for the right censored individuals is non-response up to time
  // t1_right_censored
  for (i in 1:N_right_censored) { // vectorize?
    target += log(
      (1 - pr_response) + pr_response * (1 - weibull_cdf( t1_right_censored[i], alpha, sigma ))
    );
  }

}



generated quantities {

  real<lower=0> response_time_interval_censored[N_interval_censored];
  real<lower=0> response_time_right_censored[N_right_censored];
  real<lower=0> response_time_to_be_recruited[N_to_be_recruited];
  real pr_response_cond = 0.0; // buffer variable to make code more legible
  real S_t = 0.0; // buffer variable to make code more legible

  for (i in 1:N_interval_censored) {
    // definite responder, directly sample from truncated weibull
    response_time_interval_censored[i] = ttweibull_rng(alpha, sigma, t1_interval_censored[i], t2_interval_censored[i]);
  }

  for (i in 1:N_right_censored) {
    // Pr[ responder | no response up to t ] is not constant, need to use Bayes Theorem
    S_t = 1 - weibull_cdf(t1_right_censored[i], alpha, sigma);
    pr_response_cond = S_t*pr_response / ( (1-pr_response) + pr_response*S_t );
    if (bernoulli_rng(pr_response_cond) == 1) {
      response_time_right_censored[i] = tweibull_rng(alpha, sigma, t1_right_censored[i]);
    } else {
      // non responder, set to infinity
      response_time_right_censored[i] = positive_infinity();
    }
  }

  for (i in 1:N_to_be_recruited) {
    if (bernoulli_rng(pr_response) == 1) {
      // responder
      response_time_to_be_recruited[i] = weibull_rng(alpha, sigma);
    } else {
      // non-responder
      response_time_to_be_recruited[i] = positive_infinity();
    }
  }

}

