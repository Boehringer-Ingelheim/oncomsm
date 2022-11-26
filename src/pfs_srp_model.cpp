// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;



// P(0.3 < X < 0.8), X ~ Beta(a, b)
class BetaPDF: public Func
{
private:
  double a;
  double b;
public:
  BetaPDF(double a_, double b_) : a(a_), b(b_) {}

  double operator()(const double& x) const
  {
    return R::dbeta(x, a, b, 0);
  }
};

// [[Rcpp::export]]
Rcpp::List integrate_test()
{
  const double a = 3, b = 10;
  const double lower = 0.3, upper = 0.8;
  const double true_val = R::pbeta(upper, a, b, 1, 0) -
    R::pbeta(lower, a, b, 1, 0);

  BetaPDF f(a, b);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}
