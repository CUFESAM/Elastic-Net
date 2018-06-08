#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]

void givens(arma::mat & L, arma::uword & k){
  arma::uword p = L.n_rows;
  arma::uword n = L.n_cols;
	if (p != n) 
		std::cout << "Wrong Matrix";
	if (k > p)
		std::cout << "Wrong input of k";
	L.shed_row(k);

	arma::uword mk = k;
	while(mk < p-1){
		arma::vec mx = {L(mk,mk),L(mk,mk+1)};
		double lmx   = sqrt(sum(mx % mx));
		L(mk,mk)    = lmx;
		L(mk,mk+1)  = 0;
		arma::mat gives = {{mx(0)/lmx, -mx(1)/lmx},{mx(1)/lmx,mx(0)/lmx}};

		if(mk < p-2)
			L.submat(mk+1,mk,p-2,mk+1) = L.submat(mk+1,mk,p-2,mk+1) * gives;
		mk += 1;
	}
	L.shed_col(p-1);
}
