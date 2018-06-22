#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]


void forupdate(arma::mat & L, arma::vec & xxk, double & xkxk){
  arma::vec lk = arma::solve(arma::trimatl(L), xxk);
  double lkk   = sqrt(xkxk - sum(lk % lk));

  lk.resize(lk.size()+1);
  lk[lk.size()-1] = lkk;

  arma::vec zero(L.n_rows, arma::fill::zeros);
  arma::mat LL = arma::join_rows(L,zero);

  L = arma::join_cols(LL, lk.t());
}
