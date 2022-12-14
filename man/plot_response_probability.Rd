% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_response_probability.R
\name{plot_response_probability}
\alias{plot_response_probability}
\title{Plot the response probability distributions}
\usage{
plot_response_probability(
  model,
  parameter_sample = NULL,
  seed = 42L,
  nsim = 500L,
  warmup = 250,
  nuts_control = list(),
  ...
)
}
\arguments{
\item{model}{an object of class \link{srpmodel} containing prior information}

\item{parameter_sample}{a stanfit object with samples from the respective model.}

\item{seed}{integer, fixed random seed; NULL for no fixed seed}

\item{nsim}{integer, number of samples to draw}

\item{warmup}{integer, number of warm-up samples for the MCMC sampler
before retaining samples; see \code{warmup} parameter in \code{\link[rstan:stan]{rstan::stan()}}.}

\item{nuts_control}{list, parameters for NUTS algorithm see \code{control}
argument in\code{\link[rstan:stan]{rstan::stan()}}}

\item{...}{further arguments passed to method implementations}
}
\value{
a \link[ggplot2:ggplot]{ggplot2::ggplot} object
}
\description{
\code{plot_response_probability()} plots the distribution over the response
probability parameter in the specified model.
}
\examples{
mdl <- create_srpmodel(A = define_srp_prior())
plot_response_probability(mdl)

}
\seealso{
\code{\link[=plot_transition_times]{plot_transition_times()}} \code{\link[=plot_pfs]{plot_pfs()}}
}
