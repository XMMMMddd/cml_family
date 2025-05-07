#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int addOneCpp(int x) {
  return x + 1;
}

