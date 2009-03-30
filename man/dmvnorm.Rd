\name{dmvnorm}
\alias{dmvnorm}
\alias{logdmvnorm}
\title{The Multivariate Normal Density}
\description{
  Density and log-density for the multivariate normal distribution 
  with mean equal to \code{mu}
  and variance matrix equal to \code{sigma}.
}
\usage{
dmvnorm(y, mu=NULL, sigma=NULL)
logdmvnorm(y, mu=NULL, sigma=NULL)
}
\arguments{
  \item{y}{Either a d-vector or an \eqn{n\times d}{n x d} matrix, where \eqn{d}
  is the dimension of the normal distribution and \eqn{n} is the number of
  points at which the density is to be evaluated.}
  \item{mu}{\eqn{d}-vector:  Mean of the normal distribution (or NULL uses
  the origin as default)}
  \item{sigma}{\eqn{d\times d}{d x d} matrix:  Variance matrix of the normal
  distribution (or NULL uses the identity matrix as default)}
}
\value{
  \code{dmvnorm} gives the densities,
  \code{logdmvnorm} gives the logarithm of the densities
}
\details{
This code is written to be efficient, using the qr-decomposition of the 
covariance matrix (and using it only once, rather than recalculating it
for both the determinant and the inverse of \code{sigma}).
}
\seealso{
  \code{\link{qr}}, \code{\link{qr.solve}}, \code{\link{dnorm}}, 
  \code{\link{rmvnorm}}
}
\keyword{distribution}
