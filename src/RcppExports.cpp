// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// impute_srp_model
DataFrame impute_srp_model(DataFrame df, NumericVector response_probabilities, NumericMatrix shapes, NumericMatrix scales, NumericVector visit_spacing, double max_time);
RcppExport SEXP _oncomsm_impute_srp_model(SEXP dfSEXP, SEXP response_probabilitiesSEXP, SEXP shapesSEXP, SEXP scalesSEXP, SEXP visit_spacingSEXP, SEXP max_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type response_probabilities(response_probabilitiesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type shapes(shapesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type scales(scalesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type visit_spacing(visit_spacingSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(impute_srp_model(df, response_probabilities, shapes, scales, visit_spacing, max_time));
    return rcpp_result_gen;
END_RCPP
}
// pfs
NumericVector pfs(NumericVector& t, const double& p, NumericVector& shapes, NumericVector& scales);
RcppExport SEXP _oncomsm_pfs(SEXP tSEXP, SEXP pSEXP, SEXP shapesSEXP, SEXP scalesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type t(tSEXP);
    Rcpp::traits::input_parameter< const double& >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type shapes(shapesSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type scales(scalesSEXP);
    rcpp_result_gen = Rcpp::wrap(pfs(t, p, shapes, scales));
    return rcpp_result_gen;
END_RCPP
}
// rtruncweibull
double rtruncweibull(double shape, double scale, double a, double b);
RcppExport SEXP _oncomsm_rtruncweibull(SEXP shapeSEXP, SEXP scaleSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rtruncweibull(shape, scale, a, b));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4srp_model_simple_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_oncomsm_impute_srp_model", (DL_FUNC) &_oncomsm_impute_srp_model, 6},
    {"_oncomsm_pfs", (DL_FUNC) &_oncomsm_pfs, 4},
    {"_oncomsm_rtruncweibull", (DL_FUNC) &_oncomsm_rtruncweibull, 4},
    {"_rcpp_module_boot_stan_fit4srp_model_simple_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4srp_model_simple_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_oncomsm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
