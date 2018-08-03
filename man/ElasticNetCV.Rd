\name{ElasticNetCV}
\alias{ElasticNetCV}
\title{
Cross validation
}
\description{
Computes k-fold cross-validation for elastic net.
}
\usage{
ElasticNetCV(x, y)

}
\arguments{
  \item{x}{ A data.frame or matrix of predictors}
  \item{y}{ A vector of response variables}
}

\details{
This function reads data into its environment and returns a list of three outcomes. To perform elastic net or cross-validation of elastic net, use the corresponding element of the returned list. See examples below.
The penalty of L1-norm and L2-norm is denoted by lambda1 and lambda2 respectively.
}
\value{
  \item{cv.choosemodel}{Given the parameter k folds and lambda2 (optional),  cv.choosemodel performs cross-validation  to  select  the  opti-  mal lambda1 and computes the corresponding coefficient of each variable. If lambda2 is NULL, cv.choosemodel selects the optimal lambda2 from a sequence going from 0 to 1 in steps of 0.1 and the corresponding optimal lambda1, then it returns the coefficient of each variable.}
A list of three outcomes will be returned:
  \item{Elasticnet}{Given lambda1 (optional) and lambda2, Elasticnet_ calculates an elastic net-regularized regression and returns the coefficients of each variable. If lambda1 is NULL, Elasticnet_ prints out the trace of lambda1 and the corresponding coefficient of each variable.}
  \item{output}{Prints the cross-validation outputs,  including the minimum   MSE, the coefficient of each variable, lambda1 and lambda2.}
  \item{predict}{Reads a data.frame of the testing data set and returns  predictions using the trained model.}
}

\examples{
  #Use R built-in datasets mtcars for a model fitting
  x <- mtcars[,-1]
  y <- mtcars[, 1]

  #fit model
  model <- ElasticNetCV(x,y)

  #fit a elastic net with lambda2 = 1
  model$Elasticnet_(lambda2 = 1)

  #choose model using cv
  model$cv.choosemodel(k = 31)    #Leave-one-out cross validation
  model$output()				  #See the output

  #predict
  pre <- mtcars[1:3,-1]
  model$predict(pre)
}

