% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_posterior.R, R/sample_prior.R
\name{sample_posterior}
\alias{sample_posterior}
\alias{sample_prior}
\title{Sample parameters from a model}
\usage{
sample_posterior(
  model,
  data,
  now = NULL,
  nsim = 2000L,
  seed = NULL,
  warmup = 500L,
  nuts_control = list(),
  acceptable_divergent_transition_fraction = 0.1,
  ...
)

sample_prior(
  model,
  nsim = 2000L,
  seed = NULL,
  warmup = 500L,
  nuts_control = list(),
  ...
)
}
\arguments{
\item{model}{an object of class \link{srpmodel} containing prior information}

\item{data}{a data frame with variables
\verb{subject_id<chr>} (subject identifier),
\verb{group_id<chr>} (group identifier),
\verb{t<dbl>} (time of visit, relative to first visit in study),
\verb{state<chr>} (state recorded at visit).
Allowed states are "stable", "response", "progression" (or death),
and "EOF" (end of follow-up).
The EOF state marks the end of an individual's follow-up before the absorbing
state "progression".}

\item{now}{numeric, time from first visit in data if different form last
recorded visit}

\item{nsim}{integer, number of samples to draw}

\item{seed}{integer, fixed random seed; NULL for no fixed seed}

\item{warmup}{integer, number of warm-up samples for the MCMC sampler
before retaining samples; see \code{warmup} parameter in \code{\link[rstan:stan]{rstan::stan()}}.}

\item{nuts_control}{list, parameters for NUTS algorithm see \code{control}
argument in\code{\link[rstan:stan]{rstan::stan()}}}

\item{acceptable_divergent_transition_fraction, }{numeric between 0 and 1
giving the acceptable fraction of divergent transitions before throwing an
error}

\item{...}{further arguments passed to method implementations}
}
\value{
A \link[rstan:stanfit-class]{rstan::stanfit} object with posterior samples.
}
\description{
\code{sample_posterior()} draws samples from the
posterior distribution of the specified model given a data set with
visit data.

\code{sample_prior()} draws samples from the
prior distribution of the specified model object.
}
\examples{
mdl <- create_srpmodel(A = define_srp_prior())
tbl <- tibble::tibble(
  subject_id = c("A1", "A1"),
  group_id = c("A", "A"),
  t = c(0, 1.5),
  state = c("stable", "response")
)
sample_posterior(mdl, tbl, seed = 42L)

sample_prior(mdl, seed = 42L)

}
\seealso{
\code{\link[rstan:stan]{rstan::stan()}} \code{\link[=parameter_sample_to_tibble]{parameter_sample_to_tibble()}}
\code{\link[=sample_predictive]{sample_predictive()}} \code{\link[=impute]{impute()}}
}
