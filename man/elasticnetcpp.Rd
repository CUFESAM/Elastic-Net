\name{elasticnetcpp}
\alias{elasticnetcpp}
\title{
A fast way fitting elastic net
}
\description{
ElasticNet is used to fit lasso, ridge and elastic net. For given independent and dependent variable, the best model will be chosen by MSE.
}
\usage{
elasticnetcpp(XTX, XTY, lam2, lam1 = -1)
}
\arguments{
  \item{XTX}{
  The product of the transpose of independent variable x and itself.
  }
  \item{XTY}{
  The product of the transpose of independent variable x and dependent variable y
  }
  \item{lam1}{
  The hyper-parameter which controls the degree of penalty of L1 norm.
  }
  \item{lam2}{
  The hyper-parameter which controls the degree of penalty of L2 norm.
  }
}
\details{
If lambda1 or lambda2 is set different from the default value, the explicit coefficient vector will be returned under given lambda1/2
}
\value{
A list will be returned, which contains a set of lambda1 where the variable included in the model will change. Also a set of coefficient corresponding to lambda1 will be contained in the list.
If lambda1/2 is assigned, the coefficient under which will also be returned.
}
