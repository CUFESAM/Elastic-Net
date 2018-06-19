\name{forupdate}

\alias{forupdate}

\title{

Cholesky decomposition of a matrix added to a new column

}

\description{

forupdate takes the Cholesky factor(the lower triangular matrix) of the multiplication of  the original matrix x and its transpose, the multiplication result of x transpose and the coloumn vector to be added, and the inner product of the column vector as parameters, and returns the Cholesky factor of the updated matrix.

}

\usage{

forupdate(L, xxk, xkxk)

}

\arguments{

  \item{L}{a lower triangular matrix which is the Cholesky factor of original x^T*x}

  \item{xxk}{multiplication result of original matrix and the column to be added}

  \item{xkxk}{inner product of the new column vector}

}

\value{

The lower triangular Cholesky factor of the multiplication result of the updated matrix(X) and its transpose, i.e. \eqn{L'L=X^T*X}.

}

\examples{

x <- matrix(rnorm(60),10)  #the original matrix x

##cholesky decomposition of t(x)\%*\%x

xtx <- t(x)\%*\%x

L <- t(chol(xtx))

xk <- matrix(rnorm(10),ncol = 1) #the column to be added

##parameters of function 'forupdate'

xxk <- t(x)\%*\%xk

xkxk <- sum(xk*xk)

##obtain Cholesky factor of the updated matrix t(X)\%*\%X using forupdate

L_new_f <- forupdate(L, xxk, xkxk)

##obtain Cholesky factor by  performing Cholesky decomposition on the new matrix

X <- cbind(x,xk)

XTX <- t(X)\%*\%X

L_new_c <- t(chol(XTX))

}
