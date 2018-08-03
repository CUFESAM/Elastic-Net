\name{fasterElasticNet-package}
\alias{fasterElasticNet-package}
\alias{fasterElasticNet}
\docType{package}
\title{
  Fitting ElasticNet in a fast way.
}
\description{
  FasterElasticNet uses some math algorithm such as cholesky decomposition and forward solve etc. to reduce the amount of computation. We also use Rcpp with Armadillo to improve our algorithm by speeding up almost 5 times compared by the R version.
}
\details{
  To use fasterElasticNet, dataset x(mxn) and y(mx1) should be put into the function to fit the model. Then, a completely trace of lambda1 and lambda2 can be computed if no lambda1 and lambda2 were input by using ElasticNet. Using cv.choosemodel with the number of folds will returns a best model with smallest MSE after cross-validation. Using output to print the output and predict function will return the prediction based on a new dataset.
}
\author{
Jingyi Ma

Maintainer: Linyu Zuo <zuozhe5959@gmail.com>
}
\references{
  BRADLEY, EFRON, TREVOR, HASTIE, IAIN, JOHNSTONE, AND, ROBERT, TIBSHIRANI. LEAST ANGLE REGRESSION[J]. The Annals of Statistics, 2004, 32(2): 407-499
}
\keyword{ package }
\seealso{
  https://github.com/CUFESAM/Elastic-Net
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
