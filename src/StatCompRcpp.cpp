#include <Rcpp.h>
using namespace Rcpp;
//' @title Gibbs algorithm with Rcpp.
//' @description To generate two chains 'x' and 'y' using Gibbs algorithm with Rcpp. 
//' @param N the number of samples
//' @param n0 parameter of original distribution
//' @param a0 parameter of original distribution
//' @param b0 parameter of original distribution
//' @param mu1 initial value of 'x' when generating
//' @param mu2 initial value of 'y' when generating
//' @return Return the generated chain of 'x' and 'y'
//' @importFrom stats rbeta rbinom
//' @examples
//' \dontrun{
//' N <- 5000
//' b <- 1000
//' mu1 <- 1
//' mu2 <- 0.6 
//' n0 <- 20
//' a0 <- 7
//' b0 <- 10
//' resultC <- gibbsC(N, a0, b0, n0, mu1, mu2)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, double a0, double b0, double n0, double mu1, double mu2) {
  NumericMatrix mat(N, 2);
  double x0 = 0, y0 = 0;
  mat(0, 0) = mu1;
  mat(0, 1) = mu2;
  for(int i = 1; i < N; i++) {
    y0 = mat(i-1, 1);
    mat(i, 0) = rbinom(1, n0, y0)[0];
    x0 = mat(i-1, 0);
    mat(i, 1) = rbeta(1, x0+a0, n0-x0+b0)[0];
  }
  return(mat);
}
