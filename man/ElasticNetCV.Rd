\name{ElasticNetCV}
\alias{ElasticNetCV}
\title{
Cross validation
}
\description{
Perform the cross validation of the elastic net.
}
\usage{
ElasticNetCV(x, y)
ElasticNetCV(x, y)$Elasticnet(lambda1 = -1, lambda2)
ElasticNetCV(x, y)$cv.choosemodel(lambda2 = NULL, k)
ElasticNetCV(x, y)$output()
}
\arguments{
  \item{x}{ A data.frame or matrix of features }
  \item{y}{ A vector of response}
}

\details{
This function only reads data into its environment, and returns a list of three functions. To perform elastic net or cross validation of elastic net, access the corresponding element of the returned list. See examples below.
}
\value{
A list of three functions will be returned containing:
  \item{Elasticnet}{Given lambda1 (optional) and lambda2, performs elastic net and returns the coefficients of each variable. If lambda1 is not passed, prints out the complete trajectory of lambda1 and corresponding coefficient of each variable.}
  \item{cv.choosemodel}{Given the parameter k folds and lambda2 (optional), performs cross validation to select the optimal lambda1 and coefficient of each variable. If lambda2 is not passed, selects the optimal lambda2 from 0 to 0.01 by the increment of 0.001, as well as the corresponding optimal lambda1 and coefficient of each variable. }
  \item{output}{Prints the cross validation output, including minimum MSE, coefficients of each variable, lambda1 and lambda2.}
}

\examples{
data(housing)
x <- housing[,-140]
y <- housing[,140]
cvfit <- ElasticNetCV(x,y)

# perform cross validation and print output
cvfit$cv.choosemodel(k = 20)
cvfit$output()

# perform elasticnet
cvfit$Elasticnet(lambda1 = 30,lambda2 = 0.005)

# return the trajectory of lambda1 and corresponding coefficient of each variable
cvfit$Elasticnet(lambda2 = 0.005)
}

