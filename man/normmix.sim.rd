\name{normmix.sim}
\title{Simulate from Mixtures of Normals}
\alias{normmix.sim}
\usage{
normmix.sim(n, lambda, mu, sigma, m = 1)
}

\description{
Simulate from a mixture of univariate normal distributions.
}
\arguments{
  \item{n}{Number of cases to simulate.}
  \item{lambda}{Vector of mixture probabilities, with length equal to
    desired number of components.  This is assumed to sum
    to 1; if not, it is normalized.}
  \item{mu}{Vector of means.}
  \item{sigma}{Vector of standard deviations.}
  \item{m}{Number of repeated measurements per case.}
}
\value{
  \code{normmix.sim} returns an nxm matrix in which each row is
  an i.i.d. sample from one of the components of a mixture of univariate
  normals.  Every entry of the matrix has a marginal distribution equal
  to a mixture of normals, though there is dependence among observations
  in the same row due to the fact that the component is held fixed in
  each row.
}
 
\seealso{
\code{\link{makemultdata}}
}
\examples{
##Generate data from a 2-component mixture of normals.

n<-500
lambda<-rep(1, 2)/2
mu<-c(0, 5)
sigma<-rep(1, 2)
mixnorm.data<-normmix.sim(n, lambda, mu, sigma)

##A histogram of the simulated data.

hist(mixnorm.data)
}

\keyword{file}
