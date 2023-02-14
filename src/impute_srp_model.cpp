#include <Rcpp.h>
using namespace Rcpp;



// declare rtruncweibul (could also use header file)
double rtruncweibull(double shape, double scale, double a, double b);



// calculate conditional response probability given stable up to at least time t
// [[Rcpp::export]]
double conditional_response_probability_srp(
  double t,
  double p,
  double shape_response, double shape_progression,
  double scale_response, double scale_progression
) {
  const double pr_survival_response = 1 - R::pweibull(
    t, shape_response, scale_response, 1, 0// lower = TRUE, log = FALSE
  );
  const double pr_survival_progression = 1 - R::pweibull(
    t, shape_progression, scale_progression, 1, 0// lower = TRUE, log = FALSE
  );
  // use formula of Bayes to calculate conditional probability of response
  const double pr_response = p * pr_survival_response / (
    p * pr_survival_response + (1 - p) * pr_survival_progression
  );
  return pr_response;
}



// [[Rcpp::export]]
DataFrame impute_srp_model(
    DataFrame df, // must be sorted by subject and t within subject
    NumericVector response_probabilities,
    NumericMatrix shapes,
    NumericMatrix scales,
    NumericVector visit_spacing,
    double max_time,
    CharacterVector states
) {
  // extract columns from df
  CharacterVector subject_id = df["subject_id"];
  IntegerVector        group = df["group_id"];
    NumericVector          t = df["t"];
  CharacterVector      state = df["state"];
  // save factor information
  CharacterVector group_ids = group.attr("levels");
  // check that df starts with "stable" (must be if properly sorted)
  if (state(0) != states(0)) {
    stop("df must start with initial state");  // # nocov
  }
  // compute maximal output length
  const int max_visits = floor(max_time/min(visit_spacing));
  const int n_subjects = unique(subject_id).length();
  const int n_rows_max = n_subjects * max_visits;
  // allocate memory for output
  CharacterVector subject_id_out (n_rows_max);
  IntegerVector        group_out (n_rows_max);
  NumericVector            t_out (n_rows_max);
  CharacterVector      state_out (n_rows_max);
  int row_idx = 0; // indicator recording rows-used in output
  // iterate over number of samples to draw
  float t_first_visit = t(0);
  NumericVector dt_response_interval = {NA_REAL, NA_REAL};
  int g; // group - 1 (index into parameter arrays)
  double p_response;
  bool response;
  double dt_response;
  double dt_progression;
  double tt; double dt;
  bool sample;
  for (int i = 0; i < df.nrows(); i++) {
    sample = false;
    if (i < df.nrows() - 1) {
      // is the next visit from the same subject?
      if (subject_id(i) == subject_id(i + 1)) {
        if (t(i + 1) < t(i)) {
          stop("t must be sorted within individuals"); // # nocov
        }
        if (group(i + 1) != group(i)) {
          stop("group assignment must be constant within individuals");
        }
        if ((state(i + 1) == "response") && (state(i) != "response")) {
          if (state(i) != states(0)) {
            stop("last visit before response must be initial state"); // # nocov
          }
          // record response interval
          dt_response_interval = {t(i) - t_first_visit, t(i + 1) - t_first_visit};
        }
      } else {
        // next visit is from different subject - check whether we need to sample
        sample = ((state(i) == states(0)) || (state(i) == states(1)));
      }
    } else {
      // no next visit
      sample = ((state(i) == states(0)) || (state(i) == states(0)));
    }
    // copy existing data
    subject_id_out(row_idx) = subject_id(i);
    group_out(row_idx) = group(i);
    t_out(row_idx) = t(i);
    state_out(row_idx) = state(i);
    row_idx += 1;
    // handle forward sampling if required
    if (sample) {
      g = group(i) - 1;
      // 1. determine response status
      if (state(i) == states(0)) { // still undecided, sample
        p_response = conditional_response_probability_srp(
          t(i) - t_first_visit, response_probabilities(g), shapes(g, 0), shapes(g, 1),
          scales(g, 0), scales(g, 1));
        response = R::rbinom(1, p_response);
        if (response == 1) {
          dt_response_interval = {t(i) - t_first_visit, INFINITY};
        }
      } else { // decided, state must be response
        response = true;
      }
      if (response) {
        // 3a. sample response time from truncated Weibull
        dt_response = rtruncweibull(shapes(g, 0), scales(g, 0),
                                    dt_response_interval(0),
                                    dt_response_interval(1));
        // 3b. sample progression time | response time
        dt_progression = rtruncweibull(shapes(g, 2), scales(g, 2),
                                       t(i) - t_first_visit - dt_response,
                                       INFINITY) + dt_response;
      } else {
        // 3c. non-responder - sample progression directly
        dt_response = 0.0; // no shift, see below
        dt_progression = rtruncweibull(shapes(g, 1), scales(g, 1),
                                       t(i) - t_first_visit,
                                       INFINITY);
      }
      // 4. fill with visits until progression
      tt = t(i);
      dt = visit_spacing(g);
      while (tt + dt < max_time) { // check whether there is enough time for another visit
        tt += dt;
        subject_id_out(row_idx) = subject_id(i);
        group_out(row_idx) = group(i);
        t_out(row_idx) = tt;
        if (tt < t_first_visit + dt_response) {
          state_out(row_idx) = states(0);
        } else {
          if (tt < t_first_visit + dt_progression) {
            if (response) {
              state_out(row_idx) = states(1);
            } else {
              state_out(row_idx) = states(0);
            }
          } else {
            state_out(row_idx) = states(2);
          }
        }
        row_idx += 1;
        if (state_out(row_idx - 1) == states(2)) {
          break; // no need to sample multiple visits form absorbing state
        }
      } // end while
    } // end if(sample)
    // reset counters for new subject
    // (if we needed to sample we next visit is new subject)
    if (i < df.nrows() - 2) {
      if (subject_id(i) != subject_id(i + 1)) {
        t_first_visit = t(i + 1);
        dt_response_interval = {NA_REAL, NA_REAL};
      }
    }
  } // end for
  // reduce to used memory only
  Rcpp::Range range = Rcpp::Range(0, row_idx - 1);
  subject_id_out = subject_id_out[range];
       group_out = group_out[range];
           t_out = t_out[range];
       state_out = state_out[range];
  // recover factor structure
  group_out.attr("class") = "factor";
  group_out.attr("levels") = group_ids;
  // combine to return DataFrame
  DataFrame res = DataFrame::create(
    Named("subject_id") = subject_id_out,
    Named("group_id") = group_out,
    Named("t") = t_out,
    Named("state") = state_out
  );
  return res;
}
