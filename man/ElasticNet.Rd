\name{elasticnet}
\alias{elasticnet}
\title{
A fast way fitting elastic net using RcppArmadillo
}
\description{
Elastic net is a regularization and variable selection method which linearly combines the L1 penalty of the lasso and L2 penalty of ridge methods. Based on this method, elastic- net is designed to return the trace of finding the best linear regression model. Compared with the existed R version of ElasticNet, our version speeds up the algorithm by using Cholesky decomposition, Givens rotation and RcppArmadillo.
}
\usage{
elasticnet(XTX, XTY, lam2, lam1 = -1)
}
\arguments{
  \item{XTX}{
  The product of the transpose of independent variable X and itself.
  }
  \item{XTY}{
  The product of the transpose of independent variable X and response variable Y
  }
  \item{lam1}{
  Penalty of L1-norm. No L1 penalty when lam1 = -1
  }
  \item{lam2}{
  Penalty of L2-norm, a  hyper-paramater
  }
}
\details{
When only lambda2 is given, elasticnet will return the trace of variable selection with lambda1 decreasing from lambda1_0 to zero. lambda1_0 is a value for lambda1 when there is only one predictor (the one most correlated with the response variable) in the model.

If lambda1 and lambda2 are both given, it will also return a trace. But in this case, the trace will stop when lambda1 and lambda2 reach the given ones.

To speed up the algorithm, we use some calculational tricks:

In the consideration of the low efficiency of R dealing with high-dimensional matrix, we use lower triangular matrices during the iteration of the algorithm to avoid massive matrix calculations. When adding one predictor into the model, we update XTX by recalcuting the lower triangular matrix in the Cholesky decomposition of it. While re- moving one predictor from the model, we update the lower triangular matrix with the help of Givens rotations.

Furthermore, due to the low efficiency of R dealing with loops, we rewrite the entire algorithm with RcppArmadillo, a C++ linear algebra library.
}
\value{
A list will be returned. When only lambda2 is given, the returned list contains the trace of lambda1 (relamb) and the corresponding coefficients of the predictors (reb). If both lambda1 and lambda2 are given, the corresponding coefficients of the predictors will be returned.
}

\examples{
    #Use R built-in datasets mtcars for a model fitting
    x <- as.matrix(mtcars[,-1])
    y <- as.matrix(mtcars[, 1])

    XTX <- t(x) \%*\% x
    XTY <- t(x) \%*\% y

    #Prints the output of elastic net model with lambda2 = 0
    res <- elasticnet(XTX,XTY,lam2 = 0)
}
