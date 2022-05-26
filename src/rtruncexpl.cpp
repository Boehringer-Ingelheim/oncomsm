#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rtruncexp(int n, double rate, double a, double b) {
  NumericVector smpl(n);
  NumericVector u(n);
  double p1;
  double p2;
  p1 = R::pexp(a, 1/rate, true, false);
  p2 = R::pexp(b, 1/rate, true, false);
  if (p2 - p1 < 1e-4) {
    u = rep((p2 + p1)/2.0, n);
  } else {
    u = Rcpp::runif(n, p1, p2);
  }
  for (int i = 0; i < n; i++) {
    smpl[i] = -log(1 - u[i])/rate;
  }
  return smpl;
}
