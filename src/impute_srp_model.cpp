#include <Rcpp.h>
using namespace Rcpp;



// declare rtruncweibul (could also use header file)
double rtruncweibull(double shape, double scale, double a, double b);



// calculate conditional response probability given stable up to at least time t
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
DataFrame f(
    DataFrame df, // must be sorted by subject and t within subject
    NumericVector response_probabilities,
    NumericMatrix shapes,
    NumericMatrix scales,
    NumericVector visit_spacing,
    double max_time // maximal time
) {
  // extract columns from df
  CharacterVector subject_id = df["subject_id"];
  IntegerVector        group = df["group_id"];
    NumericVector          t = df["t"];
  CharacterVector      state = df["state"];
  // save factor information
  CharacterVector group_ids = group.attr("levels");
  // check that df starts with "stable" (must be if properly sorted)
  if (state(0) != "stable") {
    stop("df must start with 'stable' state");
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
          stop("t must be sorted within individuals");
        }
        if (group(i + 1) != group(i)) {
          stop("group assignment must be constant within individuals");
        }
        if ((state(i + 1) == "response") && (state(i) != "response")) {
          if (state(i) != "stable") {
            stop("last visit before response must be 'stable'");
          }
          // record response interval
          dt_response_interval = {t(i) - t_first_visit, t(i + 1) - t_first_visit};
        }
      } else {
        // next visit is from different subject - check whether we need to sample
        sample = ((state(i) == "stable") || (state(i) == "response"));
        // make sure new subject start with 'stable'
        if (state(i + 1) != "stable") {
          stop("first visit of each individual must be 'stable' (start of treatment)");
        }
      }
    } else {
      // no next visit
      sample = ((state(i) == "stable") || (state(i) == "response"));
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
      if (state(i) == "stable") { // still undecided, sample
        p_response = conditional_response_probability_srp(
          t(i) - t_first_visit, response_probabilities(g), shapes(g, 0), shapes(g, 1),
          scales(g, 0), scales(g, 1));
        response = R::rbinom(1, p_response);
        dt_response_interval = {t(i) - t_first_visit, INFINITY};
      } else { // decided
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
          state_out(row_idx) = "stable";
        } else {
          if (tt < t_first_visit + dt_progression) {
            if (response) {
              state_out(row_idx) = "response";
            } else {
              state_out(row_idx) = "stable";
            }
          } else {
            state_out(row_idx) = "progression";
          }
        }
        row_idx += 1;
        if (state_out(row_idx - 1) == "progression") {
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


// [[Rcpp::export]]
DataFrame impute_srp_model(
    DataFrame df,
    NumericMatrix p,
    NumericVector shape_vec,
    NumericVector scale_vec,
    NumericVector visit_spacing,
    int n_sim, // number of resamples
    int n_groups
) {
  // extract columns from df
  IntegerVector subject_id = df["subject_id"];
  IntegerVector group_id = df["group_id"];
  CharacterVector from = df["from"];
  CharacterVector to = df["to"];
  NumericVector t_min = df["t_min"];
  NumericVector t_max = df["t_max"];
  NumericVector t_sot = df["t_sot"];
  // store shape/scale as multidimensional arrays
  IntegerVector dims = {n_sim, n_groups, 3};
  Function array_ind("arrayInd");
  double shape[n_sim][n_groups][3];
  double scale[n_sim][n_groups][3];
  for (int i = 0; i < shape_vec.length(); i++) {
    IntegerVector idx = array_ind(i + 1, dims);
    shape[idx[0] - 1][idx[1] - 1][idx[2] - 1] = shape_vec[i];
    scale[idx[0] - 1][idx[1] - 1][idx[2] - 1] = scale_vec[i];
  }
  IntegerVector subject_id_out (2*df.nrows()*n_sim);
  IntegerVector group_id_out (2*df.nrows()*n_sim);
  CharacterVector from_out (2*df.nrows()*n_sim);
  CharacterVector to_out (2*df.nrows()*n_sim);
  NumericVector t_min_out (2*df.nrows()*n_sim);
  NumericVector t_max_out (2*df.nrows()*n_sim);
  NumericVector t_sot_out (2*df.nrows()*n_sim);
  IntegerVector iter_out (2*df.nrows()*n_sim);
  // iterate over number of samples to draw
  int row_idx = 0;
  for (int iter = 0; iter < n_sim; iter++) {
    // iterate over rows in data
    for (int i = 0; i < df.nrows(); i++) {
      int g = group_id[i] - 1; // -1 due to 0 based indexing in C
      bool target_state_known = !CharacterVector::is_na(to[i]);
      if (target_state_known) {
        // nothing to sample - target state is observed, copy for each iteration
        subject_id_out[row_idx] = subject_id(i);
        group_id_out[row_idx] = group_id(i);
        from_out[row_idx] = from(i);
        to_out[row_idx] = to(i);
        t_min_out[row_idx] = t_min(i);
        t_max_out[row_idx] = t_max(i);
        t_sot_out[row_idx] = t_sot(i);
        iter_out[row_idx] = iter + 1;
        row_idx++;
      } else {
        if (from(i) == "stable") {
          // first sample whether response or progression happens
          double pr_response = conditional_response_probability_srp(
            t_min(i) - t_sot(i), // (min) time since start of treatment
            p(iter, g), // unconditional response probability
            // transition index -1 due to 0 based indexing in C
            shape[iter][g][0], shape[iter][g][1],
                                             scale[iter][g][0], scale[iter][g][1]
          );
          bool response = R::rbinom(1, pr_response);
          if (response) {
            double dt = rtruncweibull(
              shape[iter][g][0], scale[iter][g][0], t_min(i), R_PosInf
            );
            double dt_prog = R::rweibull(
              shape[iter][g][2], scale[iter][g][2]
            );
            // apply visit scheme
            double n_visits = floor(dt / visit_spacing(g));
            double tmin = t_min(i) + visit_spacing(g) * n_visits;
            double tmax = t_min(i) + visit_spacing(g) * (n_visits + 1);
            double n_visits_prog = floor((dt_prog + dt) / visit_spacing(g));
            double tmin_prog = t_min(i) + visit_spacing(g) * n_visits_prog;
            double tmax_prog = t_min(i) + visit_spacing(g) * (n_visits_prog + 1);
            // check corner case where transitions are in same visit interval
            if (tmin_prog == tmin) {
              // append to results vectors stable -> progression
              subject_id_out[row_idx] = subject_id(i);
              group_id_out[row_idx] = group_id(i);
              from_out[row_idx] = "stable";
              to_out[row_idx] = "progression";
              t_min_out[row_idx] = tmin;
              t_max_out[row_idx] = tmax;
              t_sot_out[row_idx] = t_sot(i);
              iter_out[row_idx] = iter + 1;
              row_idx++;
            } else {
              // append to results vectors stable -> response
              subject_id_out[row_idx] = subject_id(i);
              group_id_out[row_idx] = group_id(i);
              from_out[row_idx] = "stable";
              to_out[row_idx] = "response";
              t_min_out[row_idx] = tmin;
              t_max_out[row_idx] = tmax;
              t_sot_out[row_idx] = t_sot(i);
              iter_out[row_idx] = iter + 1;
              row_idx++;
              //append to results vectors response -> progression
              subject_id_out[row_idx] = subject_id(i);
              group_id_out[row_idx] = group_id(i);
              from_out[row_idx] = "response";
              to_out[row_idx] = "progression";
              t_min_out[row_idx] = tmin_prog;
              t_max_out[row_idx] = tmax_prog;
              t_sot_out[row_idx] = t_sot(i);
              iter_out[row_idx] = iter + 1;
              row_idx++;
            }
          } else {
            // non-responder (directly to progression)
            // sample exact time from stable to progression
            double dt = rtruncweibull(
              shape[iter][g][1], scale[iter][g][1], t_min(i), R_PosInf
            );
            // apply visit scheme
            double n_visits = floor(dt / visit_spacing(g));
            double tmin = t_min(i) + visit_spacing(g) * n_visits;
            double tmax = t_min(i) + visit_spacing(g) * (n_visits + 1);
            // append to results vectors
            subject_id_out[row_idx] = subject_id(i);
            group_id_out[row_idx] = group_id(i);
            from_out[row_idx] = "stable";
            to_out[row_idx] = "progression";
            t_min_out[row_idx] = tmin;
            t_max_out[row_idx] = tmax;
            t_sot_out[row_idx] = t_sot(i);
            iter_out[row_idx] = iter + 1;
            row_idx++;
          }
        }
        if (from(i) == "response") {
          if (from[i - 1] != "stable" ||
              subject_id[i - 1] != subject_id[i]) {
            stop("Unexpected");
          }
          double t_response = rtruncweibull(
            shape[iter][g][0], scale[iter][g][0], t_min(i-1), t_max(i-1)
          );
          double dt_progression = rtruncweibull(
            shape[iter][g][2], scale[iter][g][2],
                                             t_min(i) - t_response, R_PosInf
          );
          double n_visits_prog = floor(dt_progression / visit_spacing(g));
          double tmin_prog = t_min(i) + visit_spacing(g) * n_visits_prog;
          double tmax_prog = t_min(i) + visit_spacing(g) * (n_visits_prog + 1);
          // append results
          subject_id_out[row_idx] = subject_id(i);
          group_id_out[row_idx] = group_id(i);
          from_out[row_idx] = "response";
          to_out[row_idx] = "progression";
          t_min_out[row_idx] = tmin_prog;
          t_max_out[row_idx] = tmax_prog;
          t_sot_out[row_idx] = t_sot(i);
          iter_out[row_idx] = iter + 1;
          row_idx++;
        }
      }
    }
  }

  subject_id_out = subject_id_out[Rcpp::Range(0,(row_idx-1))];
  group_id_out = group_id_out[Rcpp::Range(0,(row_idx-1))];
  t_min_out = t_min_out[Rcpp::Range(0,(row_idx-1))];
  t_max_out = t_max_out[Rcpp::Range(0,(row_idx-1))];
  t_sot_out = t_sot_out[Rcpp::Range(0,(row_idx-1))];
  iter_out = iter_out[Rcpp::Range(0,(row_idx-1))];

  to_out = to_out[Rcpp::Range(0,(row_idx-1))];
  from_out = from_out[Rcpp::Range(0,(row_idx-1))];
  // combine to return DataFrame
  DataFrame res = DataFrame::create(
    Named("subject_id") = subject_id_out,
    Named("group_id") = group_id_out,
    Named("from") = from_out,
    Named("to") = to_out,
    Named("t_min") = t_min_out,
    Named("t_max") = t_max_out,
    Named("t_sot") = t_sot_out,
    Named("iter") = iter_out
  );
  return res;
}
