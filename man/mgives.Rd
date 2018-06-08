\name{mgives}
\alias{mgives}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Solve equation set via givens transformation}
\usage{
mgives(L, k)
}
\arguments{
\item{L}{an n*n lower triangular matrix}

\item{k}{the k-th row to be removed}
}
\value{
mgives gives the (n-1)*(n-1) dimensional lower triangular
coefficient matrix from a n*n lower triangular coefficient matrix removing some variable.
}
\description{
Using givens transformation to save unnecessary computation
}
\details{
When solve object function of linear regression, it's necessary to calculate t(X)\%*\%X, which is a lot of computation expensive with regards
to variable coming in and out in lasso. Using givens transformation to spare repeating
computation and get coefficient from the original lower triangular matrix.
}
\references{
https://www.researchgate.net/publication/228964850_Generalized_Givens_Transformation_and_its_Application
}
