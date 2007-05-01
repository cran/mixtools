\name{sqrt.mat}
\title{Calculates the Square Root of a Diagonalizable Matrix}
\alias{sqrt.mat}
\usage{
sqrt.mat(x)
}
\description{
Returns the square root of a diagonalizable matrix.
}
\arguments{
  \item{x}{An nxn diagonalizable matrix.}
}
\value{
  \code{sqrt.mat} returns the square root of \code{x}.  
}
\details{
This function is called by \code{regcr}.
}
\seealso{
\code{\link{regcr}}
}
\examples{
a<-matrix(c(1, -0.2, -0.2, 1), 2, 2)
sqrt.mat(a)

}

\keyword{internal}
