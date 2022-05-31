#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double rtruncexp(double rate, double a, double b) {
  double res;
  double u;
  double p1;
  double p2;
  p1 = R::pexp(a, 1/rate, true, false);
  p2 = R::pexp(b, 1/rate, true, false);
  if (p2 - p1 < 1e-4) {
    u = (p2 + p1)/2.0;
  } else {
    u = Rcpp::runif(1, p1, p2)(0);
  }
  res = -log(1 - u)/rate;
  return res;
}



// [[Rcpp::export]]
double rtruncweibull(double shape, double scale, double a, double b) {
  double res;
  double u;
  double p1;
  double p2;
  p1 = R::pweibull(a, shape, scale, true, false);
  p2 = R::pweibull(b, shape, scale, true, false);
  if (p2 - p1 < 1e-4) {
    u = (p2 + p1)/2.0;
  } else {
    u = Rcpp::runif(1, p1, p2)(0);
  }
  res = scale * pow(-log(1 - u), 1/shape);
  return res;
}



// [[Rcpp::export]]
DataFrame visits_to_tte_(DataFrame data, String event, String nonevent) {

  // unpack old data for faster access
  const int n_visits = data.nrow();
  CharacterVector group_id = data["group_id"];
  CharacterVector subject_id = data["subject_id"];
  NumericVector dt = data["dt"]; // delta t since start of treatment
  CharacterVector status = data["status"];
  LogicalVector eof = data["eof"];

  // create containers for new data
  const int n_subjects = unique(subject_id).length();
  CharacterVector group_id_out(n_subjects);
  CharacterVector subject_id_out(n_subjects);
  NumericVector dt_lower(n_subjects); // interval censoring
  NumericVector dt_upper(n_subjects);

  String group_lagged = 0; // factor numbering starts at 1
  String subject_lagged = 0; // factor numbering starts at 1
  float dt_lagged = 0; // factor numbering starts at 1
  bool had_response = false;
  bool had_non_response = false;
  int subject_idx = -1;
  for (int visit = 0; visit < n_visits; visit++) {
    const String group_ = group_id(visit);
    const String subject_ = subject_id(visit);
    const float dt_ = dt(visit);
    const String status_ = status(visit);
    const bool eof_ = eof(visit);

    if (visit > 0) { // keeping lagged values
      group_lagged = group_id(visit - 1);
      subject_lagged = subject_id(visit - 1);
      dt_lagged = dt(visit - 1);
    } else { // use NA if not available
      group_lagged = NA_STRING;
      subject_lagged = NA_STRING;
      dt_lagged = NA_REAL;
    }
    if (subject_ != subject_lagged) { // new individual
      subject_idx += 1; // increment index of active individual
      if (!NumericVector::is_na(dt_lagged)) {
        // handle censoring of previous subject (if exists - non NA previous subject)
        if (!had_response & !had_non_response) {
          // at risk for response, earliest time point is last dt
          group_id_out(subject_idx - 1) = group_lagged;
          subject_id_out(subject_idx - 1) = subject_lagged;
          dt_lower(subject_idx - 1) = dt_lagged;
          dt_upper(subject_idx - 1) = NA_REAL;
        }
      } else {
        // first subject, no previous subject to handle
      }
      // reset response buffer for new subject
      dt_lagged = 0.000001; // cannot be exactly zero (response = start of treatment)
      had_response = false;
      had_non_response = false;
    }
    if (had_non_response | had_response) {
      continue; // can skip until next individual starts
    }
    if ( (status_ == nonevent) | (eof_ & !had_response) ) {
      // definite non-event, eof without response is non-response
      had_non_response = true;
      group_id_out(subject_idx) = group_;
      subject_id_out(subject_idx) = subject_;
      dt_lower(subject_idx) = R_PosInf;
      dt_upper(subject_idx) = R_PosInf;
    }
    if (status_ == event) { // an event
      had_response = true;
      group_id_out(subject_idx) = group_;
      subject_id_out(subject_idx) = subject_;
      dt_lower(subject_idx) = fmin(dt_lagged, dt_);
      dt_upper(subject_idx) = dt_;
    }
    if (visit == n_visits - 1) { // last visit, no subsequent - check if censored
      if (!had_response & !had_non_response) { // at risk
        group_id_out(subject_idx) = group_;
        subject_id_out(subject_idx) = subject_;
        dt_lower(subject_idx) = dt_;
        dt_upper(subject_idx) = NA_REAL;
      }
    }
  }
  return DataFrame::create(
    Named("group_id") = group_id_out,
    Named("subject_id") = subject_id_out,
    Named("dt1") = dt_lower,
    Named("dt2") = dt_upper
  );
}
