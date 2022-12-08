#' A stable-response-progression model
#'
#' `create_model()` takes one or more prior-specifications for an
#' SRP multi-state model and combines them into a single model object.
#' Groups are still treated as independent.
#'
#' `define_srp_prior()` specifies a prior distribution for a
#' three state model (stable, response, progression) for
#' a single group.
#'
#' @name srpmodel
#' @aliases srp-model
NULL

#' @param p_mean numeric, mean of the beta prior for the response probability
#' @param p_n numeric, beta prior equivalent sample size (a + b)
#' @param p_eta numeric, robustification parameter for beta prior; actual
#' prior is (1 - eta) beta + eta; i.e., eta is the non-informative weight.
#' @param p_min numeric, minimal response probability
#' @param p_max numeric, maximal response probability
#' @param median_t_q05 numeric of length three,
#' 5% quantiles of the log-normal distributions for the
#' median time-to-next-event for the three transitions s->r, s->p, r->p.
#' @param median_t_q95 numeric of length three,
#' 95% quantiles of the log-normal distributions for the
#' median time-to-next-event for the three transitions s->r, s->p, r->p.
#' @param shape_q05 numeric of length three,
#' 5% quantiles of the log-normal distributions for the shapes of the
#' time-to-next-event distributions for the three transitions s->r, s->p, r->p.
#' @param shape_q95 numeric of length three,
#' 95% quantiles of the log-normal distributions for the shapes of the
#' time-to-next-event distributions for the three transitions s->r, s->p, r->p.
#' @param visit_spacing numeric, fixed duration between visits
#' @param recruitment_rate numeric, constant recruitment rate
#'
#' @return `define_srp_prior()` returns an object of class `srp_prior`,
#' all inputs are accessible via
#' `$x` where `x` is the name of the input argument in the function call except
#' for the two parameters `visit_spacing` and `recruitment_rate`.
#' These two parameters are saved as attributes
#' and can be retrieved directly using `attr(mdl, "visit_spacing")` and
#' `attr(mdl, "recruitment_rate")`.
#'
#' `create_srpmodel()` returns an object of class `c("srpmodel", "list")` that
#' holds information about potentially multiple groups in a compact format and
#' can be accessed using the list operator `$name`.
#' `group_id` is a character vector with the group names,
#' `maximal_time` is the maximal follow-up time since the first visit in the
#' study, `visit_spacing` is the vector of per-group difference between visits
#' (only relevant for forward sampling), `recruitment_rate` is the vector of
#' per-group recruitment rates, `stan_model` is the pre-compiled 'stan' model
#' used for inference, `states` is the vector of state names in the multi-state
#' model, and `prior` is a list of hyperparamters for the model prior with
#' elements `p`, vector, for the response probability per group,
#' `median_t` is an `c(n_groups, 3, 2)` dimensional array where
#' `median_t[i,j,1]`
#' holds the 5% quantile of the the lognormal prior on median transition time
#' for group `i` and transition `j` and `median_t[i,j,2]` the corresponding
#' upper 95% quantile. The `shape` hyperparamter has the same format and
#' specified the corresponding quantiles for the Weibull shape parameter.
#'
#'
#' @examples
#' # a model with prior 25% response rate and variance equivalent to
#' # 10 data points (i.e. a Beta(2.5, 7.5) distribution).
#' grp <- define_srp_prior(p_mean = 0.25, p_n = 10)
#' attr(grp, "recruitment_rate")
#'
#' @rdname srpmodel
#' @export
define_srp_prior <- function(
  p_mean = 0.5,
  p_n = 3,
  p_eta = 0.0,
  p_min = 0.0,
  p_max = 1.0,
  median_t_q05 = c(1, 1, 1),
  median_t_q95 = c(36, 26, 36),
  shape_q05 = rep(0.9, 3),
  shape_q95 = rep(2.5, 3),
  visit_spacing = 1, # months
  recruitment_rate = 1
) {
  checkmate::assert_vector(median_t_q05, len = 3, any.missing = FALSE)
  checkmate::assert_vector(median_t_q95, len = 3, any.missing = FALSE)
  checkmate::assert_vector(shape_q05, len = 3, any.missing = FALSE)
  checkmate::assert_vector(shape_q95, len = 3, any.missing = FALSE)
  if (visit_spacing <= 0) stop("visit spacing must be positive") # nocov
  if (recruitment_rate <= 0) stop("recruitment_rate must be positive") # nocov
  params <- as.list(environment())
  params$visit_spacing <- NULL
  params$recrutiment_rate <- NULL
  res <- structure(
    params,
    visit_spacing = visit_spacing,
    recruitment_rate = recruitment_rate,
    class = "srp_prior"
  )
  return(res)
}



#' @param ... named `srp_prior` objects; the argument names serve as
#' group labels
#' @param maximal_time the maximal overall runtime of the trial as measured from
#' the first visit of any group. No visits past this point are sampled.
#'
#' @examples
#' # a model with two groups and different priors on the respective response
#' # probabilities
#' mdl <- create_srpmodel(
#'   A = define_srp_prior(),
#'   B = define_srp_prior(p_mean = 0.33, p_n = 10)
#' )
#' mdl$median_t
#'
#' @rdname srpmodel
#' @export
create_srpmodel <- function(
  ...,
  maximal_time = 10 * 12
) {
  group_priors <- list(...)
  group_id <- names(group_priors)
  k <- length(group_id)
  if (k == 0) {
    stop("at least one group prior must be specified") # nocov
  }
  if (any(group_id == "")) {
    stop("All arguments passed to ... must be named.") # nocov
  }
  if (length(unique(group_id)) < k) {
    stop("All names of arguments passed to ... must be unique.") # nocov
  }
  transition_labels <- c("s->r", "s->p", "r->p")
  p <- matrix(NA_real_, nrow = k, ncol = 5,
              dimnames = list(group_id, c("mean", "n", "eta", "min", "max")))
  median_t <- array(NA_real_, dim = c(k, 3, 2),
                    dimnames = list(group_id, transition_labels,
                                    c("q05", "q95")))
  shape <- array(NA_real_, dim = c(k, 3, 2),
                 dimnames = list(group_id, transition_labels, c("q05", "q95")))
  visit_spacing <- numeric(k)
  names(visit_spacing) <- group_id
  recruitment_rate <- numeric(k)
  names(recruitment_rate) <- group_id
  for (g in group_id) {
    group <- group_priors[[g]]
    visit_spacing[g] <- attr(group, "visit_spacing")
    recruitment_rate[g] <- attr(group, "recruitment_rate")
    p[g, ] <- c(group$p_mean, group$p_n, group$p_eta, group$p_min, group$p_max)
    median_t[g, , ] <- cbind(group$median_t_q05, group$median_t_q95)
    shape[g, , ] <- cbind(group$shape_mode, group$shape_q05, group$shape_q95)
  }
  res <- structure(
    list(
      group_id = group_id,
      maximal_time = maximal_time,
      visit_spacing = visit_spacing,
      recruitment_rate = recruitment_rate,
      stan_model = stanmodels$srp_model_simple,
      states = c("stable", "response", "progression"),
      prior = list(
        p = p, median_t = median_t, shape = shape
      )
    ),
    class = c("srpmodel", "list"),
    parameter_names = c("p", "median_t", "shape", "scale")
  )
  check_valid(res)
  return(res)
}
