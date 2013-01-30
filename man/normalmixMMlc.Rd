\name{normalmixMMlc}
\title{EC-MM Algorithm for Mixtures of Univariate Normals 
	with linear constraints}
\alias{normalmixMMlc}
\usage{
normalmixMMlc  (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2,
		mean.constr = NULL,
	  	mean.lincstr = NULL, mean.constant = NULL,
	  	var.lincstr = NULL, gparam = NULL,
          	epsilon = 1e-08, maxit = 1000, maxrestarts=20, 
          	verb = FALSE) 
}
\description{
  Return EC-MM algorithm output for mixtures of normal distributions
  with linear constraints on the means and variances parameters,
  as in Chauveau and Hunter (2011).
  The linear constraint for the means is of the form
  \eqn{\mu = M \beta + C}, where \eqn{M} and \eqn{C} are matrix 
  and vector specified as parameters.
  The linear constraints for the variances are actually specified on
  the inverse variances, by \eqn{\pi = A \gamma}, where \eqn{\pi}
  is the vector of inverse variances, and \eqn{A} is a matrix
  specified as a parameter (see below).
}
\details{
  This is a specific "EC-MM" algorithm for normal mixtures 
  with linear constraints on the means and variances parameters.
  EC-MM here means that this algorithm is similar to 
  an ECM algorithm as in Meng and Rubin (1993), 
  except that it uses conditional MM-steps 
  instead of simple M-steps. Conditional means that it 
  alternates between maximizing with respect to the \code{mu}
  and \code{lambda} while holding \code{sigma} fixed, and maximizing with
  respect to \code{sigma} and \code{lambda} while holding \code{mu}
  fixed.  This ECM generalization of EM is forced in the case of linear constraints because there is no closed-form EM algorithm.
}
\arguments{
  \item{x}{A vector of length n consisting of the data.}

  \item{lambda}{Initial value of mixing proportions.  Automatically 
  repeated as necessary 
  to produce a vector of length \code{k}, then normalized to sum to 1.
  If \code{NULL}, then \code{lambda} is random from a uniform Dirichlet
  distribution (i.e., its entries are uniform random and then it is 
  normalized to sum to 1).}

  \item{mu}{Starting value of vector of component means.
  	  If non-NULL and a vector,
  \code{k} is set to \code{length(mu)}.  If NULL, then the initial value
  is randomly generated from a normal distribution with center(s) determined
  by binning the data.}

  \item{sigma}{Starting value of vector of component standard deviations 
  for algorithm.  Obsolete for linear constraint on the inverse variances,
  use \code{gparam} instead to specify a starting value.
  Note: This needs more precision}

  \item{k}{Number of components.  Initial value ignored unless \code{mu} 
  and \code{sigma} are both NULL.}

  \item{mean.constr}{First, simplest way to define 
  	  equality constraints on the mean parameters, given as
  a vector of length \code{k}, like in \code{\link{mvnormalmixEM}}.  
  Each vector entry helps specify the constraints,
  if any, on the corresponding mean parameter:  If \code{NA}, the corresponding
  parameter is unconstrained.  If numeric, the corresponding
  parameter is fixed at that value.  If a character string consisting of
  a single character preceded by a coefficient, such as \code{"0.5a"}
  or \code{"-b"}, all parameters using the same single character in their
  constraints will fix these parameters equal to the coefficient times
  some the same free parameter.  For instance, if 
  \code{mean.constr = c(NA, 0, "a", "-a")}, then the first mean parameter
  is unconstrained, the second is fixed at zero, and the third and forth
  are constrained to be equal and opposite in sign.
  Note: if there are no linear constraints for the variances, it is
  more efficient to use directly \code{\link{mvnormalmixEM}}.}
  
  \item{mean.lincstr}{\eqn{k\times p}
  	  matrix \eqn{M} in the linear constraint for the means
  	equation \eqn{\mu = M \beta + C}, with \eqn{p \le k}.}

  \item{mean.constant}{\eqn{k}-vector of constants \eqn{C} 
  			in the linear constraint for the means
  			equation \eqn{\mu = M \beta + C}.}

  \item{var.lincstr}{\eqn{k\times q}
  	  matrix \eqn{A} in the linear constraint for the
  	  inverse variances equation \eqn{\pi = A \gamma}, 
  	with \eqn{q \le k}.}
  
  	\item{gparam}{\eqn{q}-vector of starting values
  		for the \eqn{\gamma} parameter in the
  	  linear constraint for the inverse variances.
  	  This parameter is required.}
  
  \item{epsilon}{The convergence criterion.  
  	  Convergence is declared when the change in 
 	 the observed data log-likelihood increases by less than epsilon.}

  \item{maxit}{The maximum allowed number of iterations.}

  \item{maxrestarts}{The maximum number of restarts allowed in case of a problem
  with the particular starting values chosen due to one of the variance
  estimates getting too small
  (each restart uses randomly chosen
  starting values).  It is well-known that when each component of a normal
  mixture may have its own mean and variance, the likelihood has no maximizer;
  in such cases, we hope to find a "nice" local maximum with this algorithm
  instead, but occasionally the algorithm finds a "not nice" solution and
  one of the variances goes to zero, driving the likelihood to infinity.}

  \item{verb}{If TRUE, then various updates are printed during each 
  iteration of the algorithm.} 

}
\value{
  \code{normalmixMMlc} returns a list of class \code{mixEM} with items:
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
  \item{beta}{The final \eqn{\beta} parameter estimate.}
  \item{gamma}{The final \eqn{\gamma} parameter estimate.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
  \code{\link{mvnormalmixEM}}, \code{\link{normalmixEM2comp}}
}
\references{
  \itemize{
  \item McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, 
  John Wiley \& Sons, Inc.
  \item Meng, X.-L. and Rubin, D. B. (1993) Maximum Likelihood Estimation
  Via the ECM Algorithm:  A General Framework, \emph{Biometrika} 80(2):
  267-278.
  \item Chauveau, D. and Hunter, D.R. (2011)
  	ECM and MM algorithms for mixtures with constrained parameters,
  	\emph{HAL preprint number ???}.
  \item Thomas, H., Lohaus, A., and Domsch, H. (2011), 
  	Extensions of Reliability Theory, in 
  	Nonparametric Statistics and Mixture Models: A Festschrift in
  	Honor of Thomas Hettmansperger \emph{Singapore:  World Scientific},
  	pp. 309-316.
  }
}
\examples{
##Analyzing the Habituation data From Thomas et al (2011)
## with a 3-component mixture of normals with linear constraints.

data(Habituationdata)
attach(Habituationdata)

# to be continued...
}

\keyword{file}
