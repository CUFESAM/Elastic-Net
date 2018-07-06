\name{elasticnet}
\alias{elasticnet}
\title{
A fast way fitting elastic net using RcppArmadillo
}
\description{
ElasticNet is a regression using regularized method which linearly combines the L1 and L2 penalties of the lasso and ridge methods. For given variables and ratio of L2-norm, sometimes with givens ratio of L1-norm, Elastic will return the trace of choosing the best linear regression model. Conparing with the existed R version ElasticNet, we using givens rotation and forward update to speed up the method which is slow when dealing with big data.
}
\usage{
elasticnet(XTX, XTY, lam2, lam1 = -1)
}
\arguments{
  \item{XTX}{
  The product of the transpose of independent variable x and itself. 
  }
  \item{XTY}{
  The product of the transpose of independent variable x and dependent variable y
  }
  \item{lam1}{
  Ratio of L1-norm, which controls the degree of penalty of L1-norm. lam1 = -1 means not it is not given.
  }
  \item{lam2}{
  The hyper-parameter which controls the degree of penalty of L2-norm.
  }
}
\details{
If there is only lambda2 we put in, the ELastic Net will return a clearly entire trace of variable selection. Variable selection begin with lambda1 is as large as it only alow one variable appear in the regression model and end with lambda1 equals 0.

But if lambda1 and lambda2 are all given, it will also return the trace but will stop at where lambda1 and lambda2 are given.

To speed up the algorithm, we use some calculational tricks. Considering R is not good at solving matrix especially when dimensions of the matrix is big. So we try to build up a lower triangular matrix during the iteration of the algorithm to avoid mass matrix operation. To solving the variable's in and out based on a lower triangular matrix,  we use forward update to put variable in and givens ratation to help variable out when it is needed. In this way, we succeed in keeping matrix in lower triangular formation to speed up.

Due to the low efficiency R dealing with loops, we choose to build the entire algorithm using RcppArmadillo, a C++ linear algebra library.
}
\value{
A list will be returned, which contains a set of lambda1 where the variable included in the model will change. Also a set of coefficient corresponding to lambda1 will be contained in the list.
If lambda1/2 is assigned, the coefficients under which will also be returned.
}

\examples{
  \dontrun{
    install.packages(ElemStatLearn)
    library(ElemStatLearn)
    data(porstate)            #provide data

    x <- as.matrix(prostate[,1:8])
    y <- prostate[,9]
    XTX <- t(x) \%*\% x
    XTY <- t(x) \%*\% y
    res <- elasticnet(XTX,XTY,lam2 = 0)
  }
}