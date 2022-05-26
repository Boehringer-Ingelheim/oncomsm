#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rtruncweibull(int n, double shape, double scale, double a, double b) {
  NumericVector smpl(n);
  NumericVector u(n);
  double p1;
  double p2;
  p1 = R::pweibull(a, shape, scale, true, false);
  p2 = R::pweibull(b, shape, scale, true, false);
  if (p2 - p1 < 1e-4) {
    u = rep((p2 + p1)/2.0, n);
  } else {
    u = Rcpp::runif(n, p1, p2);
  }
  for (int i = 0; i < n; i++) {
    smpl[i] = scale * pow(-log(1 - u[i]), 1/shape);
  }
  return smpl;
}
