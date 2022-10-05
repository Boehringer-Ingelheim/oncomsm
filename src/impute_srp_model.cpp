#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double impute_srp_model(DataFrame df, NumericVector shape, int nsim, int ngroups) {
  IntegerVector subject_id = df["subject_id"];
  IntegerVector group_id = df["group_id"];
  IntegerVector dim = {nsim, ngroups, 3};
  Function array_ind("arrayInd");
  double shape_arr[nsim][ngroups][3];

  for (int i = 0; i < shape.length(); i++) {
    IntegerVector idx = array_ind(i + 1, dim);
    shape_arr[idx[0] - 1][idx[1] - 1][idx[2] - 1] = shape[i];
  }
  return shape_arr[1][1][1];
}
