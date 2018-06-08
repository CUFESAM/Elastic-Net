\name{msig}
\alias{msig}
\title{
Construct a special symmetric matrix
}
\description{
Creat a p*p dimensional symmetric matrix where the diagonal elements are one and
  \eqn{
    X_ij = \rho^| i - j |
  }
}
\usage{
msig(p, rho)
}
\arguments{
  \item{p}{
  Demension of the matrix you want.
  }
  \item{rho}{
  Base number.
  }
}
\value{
A p*p dimensional matrix will be returned.
}
\examples{
msig(1, 1)
msig(1, 5)
msig(5, 1)
msig(5, 5)
}

