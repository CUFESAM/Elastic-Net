\name{givens}
\alias{givens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get cholesky decomposition with givens transformation.}
\description{
The cholesky decomposition after deleting the k-th variable of X can be obtained by givens transformation in order to save unnecessary computation.
}

\usage{
givens(L, k)
}
\arguments{
\item{L}{an n*n lower triangular matrix,which is the cholesky decomposition of t(X)\%*\%X.}

\item{k}{the k-th row to be removed.}
}


\details{
When solve object function of linear regression, it's necessary to calculate t(X)\%*\%X, which is a lot of computation expensive with regards
to variable coming in and out in lasso. Using givens transformation to spare repeating
computation and get coefficient from the original lower triangular matrix.
}
\value{
givens gives the (n-1)*(n-1) dimensional lower triangular
coefficient matrix from a n*n lower triangular coefficient matrix removing some variable.
}
\references{
https://www.researchgate.net/publication/228964850_Generalized_Givens_Transformation_and_its_Application
}
\examples{
L <- matrix(rnorm(16),c(4,4))
L[upper.tri(L)]=0

givens(L,2)
}
