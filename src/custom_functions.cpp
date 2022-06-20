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
