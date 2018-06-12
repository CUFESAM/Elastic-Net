\name{givens}
\alias{givens}
\title{
Construct a givens transformation matrix
}
\description{
In givens transformation, if we want a vector to rotate a angle \eqn{\theta} clockwise, we are supposed to construct a matrix with format of:
  \deqn{|cos(\theta) sin(\theta)|}
  \deqn{|-sin(\theta) cos(\theta)|}

}
\usage{
givens(L, k)
}
\arguments{
  \item{L}{two numbers}
  \item{k}{square root of the sumption of the square of the mx}
}
