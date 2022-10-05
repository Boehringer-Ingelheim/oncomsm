// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

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
// impute_srp_model
double impute_srp_model(DataFrame df, NumericVector shape, int nsim, int ngroups);
RcppExport SEXP _oncomsm_impute_srp_model(SEXP dfSEXP, SEXP shapeSEXP, SEXP nsimSEXP, SEXP ngroupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    rcpp_result_gen = Rcpp::wrap(impute_srp_model(df, shape, nsim, ngroups));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4srp_model_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_oncomsm_rtruncweibull", (DL_FUNC) &_oncomsm_rtruncweibull, 4},
    {"_oncomsm_impute_srp_model", (DL_FUNC) &_oncomsm_impute_srp_model, 4},
    {"_rcpp_module_boot_stan_fit4srp_model_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4srp_model_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_oncomsm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
