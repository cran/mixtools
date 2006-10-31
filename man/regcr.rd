\name{regcr}
\title{Add a Bayesian Credible Region for Regression Lines to a Scatterplot}
\alias{regcr}
\usage{
regcr(beta, x, alpha = .05, nonparametric = FALSE, plot = FALSE,
      xyaxes = TRUE, ...)
}

\description{
Produce a credible region for regression lines based on
a sample of beta parameters (e.g., a sample from the posterior
distribution).  The beta parameters are the intercept and slope
from a simple linear regression.
}
\arguments{
  \item{beta}{An nx2 matrix of regression parameters.  The first column
gives the intercepts and the second column gives the slopes.} 
  \item{x}{An n-vector of the predictor variable which is necessary when
  nonparametric = TRUE.}
  \item{alpha}{The proportion of the beta sample to remove. In other
words, 1-alpha is the level of the credible region.}
  \item{nonparametric}{
If nonparametric = TRUE, then the region is based on the convex
hull of the remaining beta after trimming, which is accomplished
using a data depth technique.
If nonparametric = FALSE, then the region is based on the
asymptotic normal approximation.}
  \item{plot}{If plot = TRUE, lines are added to the existing plot.
  The type of plot created depends on the value of xyaxes.}
  \item{xyaxes}{If xyaxes = TRUE and plot = TRUE, then a credible region
    for the regression lines is plotted on the x-y axes, presumably
    overlaid on a scatterplot of the data.  If xyaxes = FALSE and
    plot = TRUE, the (convex) credible region for the regression line is
    plotted on the beta, or intercept-slope, axes, presumably overlaid
    on a scatterplot of beta.}
  \item{...}{Graphical parameters passed to \code{lines} or \code{plot}
    command.}   
}
\value{
  \code{regcr} returns a list containing the following items:
  \item{boundary}{A matrix of points in beta, or intercept-slope, space
    arrayed along the boundary of the credible region.}
  \item{upper}{A matrix of points in x-y space arrayed along the upper
    crebible limit for the regression line.}
  \item{lower}{A matrix of points in x-y space arrayed along the lower
    credible limit for the regression line.}
}
 
\seealso{
\code{\link{regmixEM}}, \code{\link{regmixMH}}
}
\examples{
## Nonparametric credible regions for mixtures of regressions fit to NOdata. 

data(NOdata)
attach(NOdata)
beta<-matrix(c(1.3, -0.1, 0.6, 0.1), 2, 2)
sigma<-c(.02, .05)
MH.out<-regmixMH(Equivalence, NO, beta = beta, s = sigma, 
                 sampsize = 10000, omega = .0013)
attach(data.frame(MH.out$theta))
beta.c1<-cbind(beta0.1[9950:9999], beta1.1[9950:9999])
beta.c2<-cbind(beta0.2[9950:9999], beta1.2[9950:9999])
plot(NO, Equivalence)
regcr(beta.c1, x = NO, nonparametric = TRUE, plot = TRUE, 
      col = 2)
regcr(beta.c2, x = NO, nonparametric = TRUE, plot = TRUE, 
      col = 3)

}


\keyword{file}
