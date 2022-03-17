#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame recist_to_tte_(DataFrame data) {

  // unpack old data for faster access
  const int n_visits = data.nrow();
  CharacterVector subject_id = data["subject_id"];
  NumericVector dt = data["dt"]; // delta t since start of treatment
  CharacterVector status = data["status"];
  LogicalVector eot = data["eot"];

  // create containers for new data
  const int n_subjects = unique(subject_id).length();
  CharacterVector subject_id_out(n_subjects);
  NumericVector dt_lower(n_subjects); // interval censoring
  NumericVector dt_upper(n_subjects);

  String subject_lagged = 0; // factor numbering starts at 1
  float dt_lagged = 0; // factor numbering starts at 1
  bool had_response = false;
  bool had_non_response = false;
  int subject_idx = -1;
  for (int visit = 0; visit < n_visits; visit++) {
    const String subject_ = subject_id(visit);
    const float dt_ = dt(visit);
    const String status_ = status(visit);
    const bool eot_ = eot(visit);

    if (visit > 0) { // keeping lagged values
      subject_lagged = subject_id(visit - 1);
      dt_lagged = dt(visit - 1);
    } else { // use NA if not available
      subject_lagged = NA_STRING;
      dt_lagged = NA_REAL;
    }
    if (subject_ != subject_lagged) { // new individual
      subject_idx += 1; // increment index of active individual
      if (!NumericVector::is_na(dt_lagged)) {
        // handle censoring of previous subject (if exists - non NA previous subject)
        if (!had_response & !had_non_response) {
          // at risk for response, earliest time point is last dt
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
    if ( (status_ == "PD") | (eot_ & !had_response) ) {
      // definite non-responder, eot without response is non-response
      had_non_response = true;
      subject_id_out(subject_idx) = subject_;
      dt_lower(subject_idx) = R_PosInf;
      dt_upper(subject_idx) = R_PosInf;
    }
    if ((status_ == "PR") | (status_ == "CR")) { // a responder
      had_response = true;
      subject_id_out(subject_idx) = subject_;
      dt_lower(subject_idx) = fmin(dt_lagged, dt_);
      dt_upper(subject_idx) = dt_;
    }
    if (visit == n_visits - 1) { // last visit, no subsequent - check if censored
      if (!had_response & !had_non_response) { // at risk
        subject_id_out(subject_idx) = subject_;
        dt_lower(subject_idx) = dt_;
        dt_upper(subject_idx) = NA_REAL;
      }
    }
  }
  return DataFrame::create(
    Named("subject_id") = subject_id_out,
    Named("dt1") = dt_lower,
    Named("dt2") = dt_upper
  );
}
