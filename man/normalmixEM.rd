\name{normalmixEM}
\title{EM Algorithm for Mixtures of Univariate Normals}
\alias{normalmixEM}
\usage{
normalmixEM(x, lambda = NULL, mu = NULL, sigma = NULL, k = 2,
            arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
            maxit = 1000, maxrestarts=20, verb = FALSE, fast=FALSE)
}
\description{
  Return EM algorithm output for mixtures of normal distributions.
}
\arguments{
  \item{x}{A vector of length n consisting of the data.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1.  This determines number of components.  If NULL, then \code{lambda} is
    random from uniform Dirichlet and number of
    components is determined by \code{mu}.}
  \item{mu}{A k-vector of initial values for the mean parameters.  
    If NULL, then the vector is generated from a normal distribution according to a binning method
    done on the data. If both \code{lambda} and \code{mu} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{A k-vector of initial values for the standard deviation parameters.  
    If NULL, then 1/\code{sigma}$^2$ has random exponential entries according to a binning method done on the data.  
    If \code{lambda}, \code{mu}, and \code{sigma} are
    NULL, then number of components is determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless \code{lambda}, \code{mu}, and \code{sigma}
    are all NULL.}
  \item{arbmean}{If TRUE, then the component densities are allowed to have different \code{mu}s. If FALSE, then
  a scale mixture will be fit.}
  \item{arbvar}{If TRUE, then the component densities are allowed to have different \code{sigma}s. If FALSE, then
  a location mixture will be fit.}
  \item{epsilon}{The convergence criterion.  Convergence is declared when the change in 
  the observed data log-likelihood increases by less than epsilon.}
  \item{maxit}{The maximum number of iterations.}
  \item{maxrestarts}{The maximum number of restarts allowed in case of a problem
  with the particular starting values chosen (each restart uses randomly chosen
  starting values).}
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
  \item{fast}{If TRUE and k==2 and arbmean==TRUE, then use 
  \code{\link{normalmixEM2comp}}, which is a much faster version of the EM 
  algorithm for this case.
  This version is less protected against certain kinds of underflow
  that can cause numerical problems and it does not permit any restarts.  If
  k>2, \code{fast} is ignored.}
}
\value{
  \code{normalmixEM} returns a list of class \code{mixEM} with items:
  \item{x}{The raw data.}
  \item{lambda}{The final mixing proportions.}
  \item{mu}{The final mean parameters.}
  \item{sigma}{The final standard deviations. If \code{arbmean} = FALSE, then only the smallest standard
   deviation is returned. See \code{scale} below.}
  \item{scale}{If \code{arbmean} = FALSE, then the scale factor for the component standard deviations is returned.
   Otherwise, this is omitted from the output.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.  This vector
  includes both the initial and the final values; thus, the number of iterations 
  is one less than its length.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
  \code{\link{mvnormalmixEM}}, \code{\link{normalmixEM2comp}}
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
}
\examples{
##Analyzing the Old Faithful geyser data with a 2-component mixture of normals.

data(faithful)
attach(faithful)
system.time(out<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03))
out
system.time(out2<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03, fast=TRUE))
out2 # same thing but much faster
}

\keyword{file}
