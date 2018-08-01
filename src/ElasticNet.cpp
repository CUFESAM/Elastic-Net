#include <RcppArmadillo.h>
#include "./givens.h"
#include "./forupdate.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]


void PushBack(arma::uvec & A, arma::uword & j){
	int p = A.size();
	A.resize(p+1);
	A(p)  = j;
}

// [[Rcpp::export]]

Rcpp::List elasticnet(arma::mat & XTX, arma::vec & XTY,
						double lam2 ,double lam1 = -1){

	int p = XTX.n_rows;
	arma::vec w    = sqrt(XTX.diag(0));
	arma::vec w1   = 1/w;
	arma::vec w2   = XTX.diag();
	double alpha0  = 1/(1 + lam2);
	double alpha1  = 1 - alpha0;
	arma::mat WW   = arma::diagmat(w2);
	arma::mat XTXW = alpha0 * XTX + alpha1 * WW;

	Rcpp::LogicalVector A(p);
	A.fill(false);

	arma::vec absXTY = arma::abs(XTY) % w1;
	arma::uword j    = absXTY.index_max();
	A(j) = !A(j);

	arma::uvec VA = {j};
	arma::mat L   = {w(j)};
	double lamb   = w1(j) * std::abs(XTY(j));

	arma::vec b(p);
	b.fill(0);
	Rcpp::NumericVector relamb = {lamb};
	arma::mat reb(b.t());

	arma::uvec VAm = arma::linspace<arma::uvec>(0,p-1,p);
	VAm.shed_row(j);

	arma::vec CC(p), SCC(p);

	while(true){
		CC  = w1 % (XTY - XTXW * b);
		SCC = arma::sign(CC);

		arma::vec SCCA = SCC(VA);
		arma::vec td   = arma::solve(trimatl(L),w(VA) % SCCA);
		arma::vec d    = arma::solve(trimatu(L.t()),td);
		arma::mat XTX_ = XTX.submat(VAm, VA);
		arma::vec a    = alpha0 * w1(VAm) % (XTX_ * d);
		arma::vec gam(p, arma::fill::zeros);
		arma::vec ww   = -b(VA)/d;

		for(int i = 0; i < VA.size(); i++){
			if(ww(i) > 0 && ww(i) < lamb)
				gam(VA(i)) = ww(i);
			else
				gam(VA(i)) = lamb;
		}

		if(sum(A) < p){
			for(int i = 0; i < (p-VA.size()); i++){
				if(a(i) * lamb <= CC(VAm(i)))
					gam(VAm(i)) = (lamb - CC(VAm(i)))/(1-a(i));
				else
					gam(VAm(i)) = (lamb + CC(VAm(i)))/(1+a(i));
			}
		}

		j = gam.index_min();
		double gammin = gam(j);
		if(lam1 >= 0 && lamb - gammin <= lam1){
			b(VA) = b(VA) + (lamb - lam1) * d;
			Rcpp::List res2 = Rcpp::List::create(Rcpp::Named("b") = b);
			return(res2);
		}
		b(VA) = b(VA) + gammin * d;
		lamb  = lamb - gammin;
		relamb.push_back(lamb);
		reb   = arma::join_cols(reb,b.t());

		if(lamb == 0)
			break;

		arma::uvec jj;
		for(arma::uword i = 0; i < VA.size(); i++){
			if(VA[i] == j)
				PushBack(jj, i);
		}

		if(jj.size() == 0){
			arma::uvec j_   = {j};
			arma::vec XTXAJ = alpha0 * XTX(VA, j_);
			double XTXjj    = w2(j);

			forupdate(L,XTXAJ,XTXjj);
			PushBack(VA,j);

			arma::uvec del = find(VAm == j);
			VAm.shed_row(del(0));
		}
		else{
			givens(L,jj(0));
			PushBack(VAm, VA(jj(0)));
			VA.shed_row(jj(0));
		}
		A(j) = !A(j);
	}
	Rcpp::List res = Rcpp::List::create(Rcpp::Named("reb") = reb , Rcpp::Named("relamb") = relamb);
	return(res);
}
