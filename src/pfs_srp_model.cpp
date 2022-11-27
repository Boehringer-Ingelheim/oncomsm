// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <Rcpp.h>
#include <RcppNumerical.h>
using namespace Rcpp;
using namespace Numer;



class SRPModelPFSIntegrand: public Func
{
private:
  double t;
  double p;
  double scale_sr, shape_sr; // stable->response
  double scale_rp, shape_rp; // response->progression
public:
  SRPModelPFSIntegrand(double t_, double p_,
                       double scale_sr_, double shape_sr_,
                       double scale_rp_, double shape_rp_) :
    t(t_), p(p_),
    scale_sr(scale_sr_), shape_sr(shape_sr_),
    scale_rp(scale_rp_), shape_rp(shape_rp_) {}

  double operator()(const double& x) const
  {
    return R::dweibull(x, shape_sr, scale_sr, 0) * // 0 = not log density
      R::pweibull(t - x, shape_rp, scale_rp, 1, 0); // 1, 0 = lower tail + not log
  }
};

// [[Rcpp::export]]
NumericVector pfs(NumericVector& t, const double& p, NumericVector& shapes, NumericVector& scales) {
  double err_est;
  int err_code;
  const int n = t.length();
  NumericVector pfs(n);
  double ind;
  double d;
  for (int i = 0; i < n; i++) {
    // define integrand for indirect path
    SRPModelPFSIntegrand f(t[i], p, scales[0], shapes[0], scales[2], shapes[2]);
    // integrate
    ind = integrate(
      f, 0.0, t[i], err_est, err_code,
      100, // max subdivisions # nocov
      1e-5, // absolute error # nocov
      1e-6, // relative error # nocov
      Integrator<double>::GaussKronrod15 // integration rule # nocov
    );
    // calculate direct progression probability | no response
    d = R::pweibull(t[i], shapes[1], scales[1], 1, 0);
    // 1 - (pr[progression | response] pr[response + pr[progression | no response] pr[no response]])
    pfs[i] = 1 - (p * ind + (1 - p) * d);
  }
  return pfs;
}
