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
DataFrame impute_srp_model(
    DataFrame df,
    NumericMatrix p,
    NumericVector shape_vec,
    NumericVector scale_vec,
    NumericVector visit_spacing,
    int n_sim,
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
  // define output vectors of length 0
  IntegerVector subject_id_out (0);
  IntegerVector group_id_out (0);
  CharacterVector from_out (0);
  CharacterVector to_out (0);
  NumericVector t_min_out (0);
  NumericVector t_max_out (0);
  NumericVector t_sot_out (0);
  IntegerVector iter_out (0);
  // iterate over number of samples to draw
  for (int iter = 0; iter < n_sim; iter++) {
    // iterate over rows in data
    for (int i = 0; i < df.nrows(); i++) {
      int g = group_id[i] - 1; // -1 due to 0 based indexing in C
      bool target_state_known = !CharacterVector::is_na(to[i]);
      if (target_state_known) {
        // nothing to sample - target state is observed, copy for each iteration
        subject_id_out.push_back(subject_id(i));
        group_id_out.push_back(group_id(i));
        from_out.push_back(from(i));
        to_out.push_back(to(i));
        t_min_out.push_back(t_min(i));
        t_max_out.push_back(t_max(i));
        t_sot_out.push_back(t_sot(i));
        iter_out.push_back(iter + 1); // offset 0-based indexing
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
            double n_visits_prog = floor(dt_prog / visit_spacing(g));
            double tmin_prog = t_min(i) + visit_spacing(g) * n_visits_prog;
            double tmax_prog = t_min(i) + visit_spacing(g) * (n_visits_prog + 1);
            // append to results vectors stable -> response
            subject_id_out.push_back(subject_id(i));
            group_id_out.push_back(group_id(i));
            from_out.push_back("stable");
            to_out.push_back("response");
            t_min_out.push_back(tmin);
            t_max_out.push_back(tmax);
            t_sot_out.push_back(t_sot(i));
            iter_out.push_back(iter);
            //append to results vectors response -> progression
            subject_id_out.push_back(subject_id(i));
            group_id_out.push_back(group_id(i));
            from_out.push_back("response");
            to_out.push_back("progression");
            t_min_out.push_back(tmin_prog);
            t_max_out.push_back(tmax_prog);
            t_sot_out.push_back(t_sot(i));
            iter_out.push_back(iter);
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
            subject_id_out.push_back(subject_id(i));
            group_id_out.push_back(group_id(i));
            from_out.push_back("stable");
            to_out.push_back("progression");
            t_min_out.push_back(tmin);
            t_max_out.push_back(tmax);
            t_sot_out.push_back(t_sot(i));
            iter_out.push_back(iter);
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
          subject_id_out.push_back(subject_id(i));
          group_id_out.push_back(group_id(i));
          from_out.push_back("response");
          to_out.push_back("progression");
          t_min_out.push_back(tmin_prog);
          t_max_out.push_back(tmax_prog);
          t_sot_out.push_back(t_sot(i));
          iter_out.push_back(iter);
        }
      }
    }
  }
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
