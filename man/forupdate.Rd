\name{forupdate}
\alias{forupdate}
\title{
Choleski decomposition on extra element
}
\description{
forupdate takes the original choleski decomposition of x, variable L in this case, outer product of original element and extra element, variable xxk in this case, and outer product of extra element itself, variable xkxk in this case, then gives the choleski decomposition of combination of original element and extra element.
}
\usage{
forupdate(L, xxk, xkxk)
}
\arguments{
  \item{L}{lower triangular matrix in choleski decomposition of original element x*x^T}
  \item{xxk}{outer product of original element and extra element}
  \item{xkxk}{outer product of extra element itself}
}
\value{
  gives returns a (p + px)*(p + px) dimensional lower triangular matrxi.
}
\examples{
#x = matrix(rnorm(5), ncol = 1)
#xk = matrix(rnorm(1), ncol = 1)
#xxk = x\%*\%t(xk)
#xkxk = xk\%*\%t(xk)
#L = chol(x\%*\%t(x))
#chol(rbind(x, xk)\%*\%t(rbind(x, xk)))
#forupdate(L, xxk, xkxk)
}
