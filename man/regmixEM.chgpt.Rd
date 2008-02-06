\name{regmixEM.chgpt}
\title{EM Algorithm for Mixtures of Regressions with a Changepoint}
\alias{regmixEM.chgpt}
\usage{
regmixEM.chgpt(y, x, lambda = NULL, gamma = NULL, beta = NULL,
               sigma = NULL, k = 2, T = 3, t = NULL, 
               epsilon = 1e-08, maxit = 10000, verb = FALSE)
}

\description{
  Returns EM algorithm output for a mixture of a regression with a changepoint and a simple
  linear regression.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An n-vector of predictor values.  A column of ones is automatically appended to x.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1. If NULL, then \code{lambda} is
    random from uniform Dirichlet.}
  \item{gamma}{Initial value of \code{gamma} parameters for the changepoint component.  Should be a 3-dimenstional vector.
    If NULL, then \code{gamma} has standard normal entries according to a binning method done on the data.}
  \item{beta}{Initial value of \code{beta} parameters for the simple linear regression component.  
    Should be a 2-dimenstional vector. If NULL, then \code{beta} has standard normal entries according to a binning method done on the data.}
  \item{sigma}{A vector of standard deviations.  If NULL, then 1/\code{sigma}$^2$ has
    random standard exponential entries according to a binning method done on the data.}
  \item{k}{Number of components.  Currently, this value must be set equal to 2.}
  \item{T}{The number of values to leave off for the range of all possible changepoints to be tested.}
  \item{t}{Initial value of the changepoint to consider.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{regmixEM.chgpt} returns a list of class \code{mixEM} with items:
  \item{lambda}{The final mixing proportions.}
  \item{gamma}{The final regression coefficients for the changepoint component.}
  \item{beta}{The final regression coefficients for the simple linear regression component.}
  \item{sigma}{The final standard deviations.}
  \item{cutpoint}{The estimated value for the changepoint in the changepoint component.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{A nx2 matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{regmixEM}}
}
\examples{
## EM output for simulated data.
 
w<-rbinom(100, 1, .5)
cpt<-50
x<-sort(runif(100, 0, 10))
x1<-cbind(1, x)
xt<-cbind(x1, (x-x[cpt])*(x>x[cpt]))
beta<-c(5, -1)
gamma<-c(15, -1, 2)
y<-w*rnorm(100, mean = xt\%*\%gamma, sd = .1) + 
   (1-w)*rnorm(100, mean = x1\%*\%beta, sd = .1)
out<-regmixEM.chgpt(y = y, x = x, t = cpt, beta = beta, 
                    gamma = gamma, verb = TRUE, epsilon = 1e-03)
out
}


\keyword{file}
