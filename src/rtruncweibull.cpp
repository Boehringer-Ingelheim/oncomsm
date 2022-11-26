#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double rtruncweibull(double shape, double scale, double a, double b) {
  double res;
  double u;
  double p1;
  double p2;
  p1 = R::pweibull(a, shape, scale, true, false);
  p2 = R::pweibull(b, shape, scale, true, false);
  u = R::runif(p1, p2);
  u = std::min(1.0 - 1e-4, std::max(1e-4, u));
  res = scale * pow(-log(1 - u), 1/shape);
  res = std::min(b, std::max(a, res));
  return res;
}
