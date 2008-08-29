\name{compCDF}
\title{Plot the Component CDF}
\alias{compCDF}
\usage{
compCDF(x, weights)
}
\description{
Plot the components' CDF via the posterior probabilities.
}
\arguments{
  \item{x}{A matrix containing the raw data. Rows are subjects and columns
   are repeated measurements.}
  \item{weights}{The weights to compute the empirical CDF; however, most of
   time they are the posterior probabilities.}
}
\value{
  \code{compCDF} returns an object which is a list with components:
  \item{result}{The component means and standard deviations for a k-component mixture.}
  \item{plot}{The plotted component CDF.}
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
  
  Elmore, R. T., Hettmansperger, T. P. and Xuan, F. (2004) The Sign Statistic, One-Way Layouts
  and Mixture Models, \emph{Statistical Science} \bold{19(4)}, 579--587.
}
\seealso{
\code{\link{makemultdata}}, \code{\link{multmixmodel.sel}}, \code{\link{multmixEM}}.
}
\examples{
## The sulfur content of the coal seams in Texas

A<-c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B<-c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59) 
C<-c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49) 
D<-c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32) 
E<-c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)

## dis.coal<-makemultdata(A, B, C, D, E, 
##                        cuts = median(c(A, B, C, D, E)))
## temp<-multmixEM(dis.coal$y, lambda = dis.coal$lambda, 
##                 theta = dis.coal$theta)

## Now plot the components' CDF via the posterior probabilities

## compCDF(dis.coal$x, temp$posterior)
}

\keyword{file}

