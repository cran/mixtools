\name{npEMindrep}
\title{Nonparametric EM-like Algorithm for Mixtures of Independent Repeated Measurements}
\alias{npEMindrep}
\usage{
npEMindrep(x, mu0, blockid = 1:ncol(x), bw=bw.nrd0(as.matrix(x)), h=bw,
           eps=1e-8, maxiter = 100, stochastic = FALSE, verb = TRUE)
}

\description{
  Returns nonparametric EM algorithm output (Benaglia et al, 2007) for mixtures
  of multivariate (repeated measures) data where the coordinates of a row (case)
  in the data matrix are assumed to be independent, conditional on the mixture
  component (subpopulation) from which they are drawn.
}
\arguments{
  \item{x}{An \eqn{n\times r}{n x r} matrix of data.  Each of the \eqn{n} rows is a case,
  and each case has \eqn{r} repeated measurements.  These measurements are assumed
  to be conditionally independent, conditional on the mixture component (subpopulation) 
  from which the case is drawn.}
  \item{mu0}{Either an \eqn{m\times r}{m x r} matrix specifying the initial
      centers for the \link{kmeans} function, or an integer \eqn{m} specifying the 
      number of initial centers, which are then choosen randomly in
  \link{kmeans}}
  \item{blockid}{A vector of length \eqn{r} identifying coordinates 
  (columns of \code{x}) that are
  assumed to be identically distributed (i.e., in the same block).  For instance,
  the default has all distinct elements, indicating that no two coordinates 
  are assumed identically distributed and thus a separate set of \eqn{m} 
  density estimates is produced for each column of \eqn{x}.  On the other hand,
  if \code{blockid=rep(1,ncol(x))}, then the coordinates in each row 
  are assumed conditionally i.i.d.}
  \item{bw}{Bandwidth for density estimation, equal to the standard deviation 
  of the kernel density.}
  \item{h}{Alternative way to specify the bandwidth, to provide backward 
  compatibility.}
  \item{eps}{Tolerance limit for declaring algorithm convergence.  Convergence
  is declared whenever the maximum change in any coordinate of the 
  \code{lambda} vector (of mixing proportion estimates) does not exceed 
  \code{eps}.}

  \item{maxiter}{The maximum number of iterations allowed, for both 
  stochastic and non-stochastic versions; 
  for non-stochastic algorithms (\code{stochastic = FALSE}), convergence
  may be declared before \code{maxiter} iterations (see \code{eps} above).}
  
  \item{stochastic}{Flag, if FALSE (the default), runs the non-stochastic version
      of the npEM algorithm, as in Benaglia et al (2007). Set to TRUE to
      run a stochastic version which simulates the posteriors at each
      iteration, and runs for \code{maxiter} iterations.}

  \item{verb}{If TRUE, print updates for every iteration of the algorithm as
  it runs}
}
\value{
\code{npEMindrep} returns a list of class \code{npEM} with the following items:
  \item{data}{The raw data (an \eqn{n\times r}{n x r} matrix).}
  \item{posteriors}{An \eqn{n\times m}{n x m} matrix of posterior probabilities for
   observation. If \code{stochastic = TRUE}, this matrix is computed 
   from an average over the \code{maxiter} iterations.}
  \item{bandwidth}{Same as the \code{bw} input argument, returned because this
  information is needed by any method that produces density estimates from the
  output.}
  \item{blockid}{Same as the \code{blockid} input argument, but recoded to have
  positive integer values.  Also needed by any method that produces density 
  estimates from the output.}
  \item{lambda}{The sequence of mixing proportions over iterations.}
  \item{lambdahat}{The final mixing proportions if \code{stochastic = FALSE}, 
  or the average mixing proportions if \code{stochastic = TRUE}.}
}
\seealso{
\code{\link{plot.npEM}}, \code{\link{normmixrm.sim}}, \code{\link{spEMsymloc}},
\code{\link{plotseq.npEM}}
}
\references{
\item Bordes, L., Chauveau, D., and Vandekerkhove, P. (2007),
An EM algorithm for a semiparametric mixture model, 
Computational Statistics and Data Analysis, 51: 5429-5443.

\item Benaglia, T., Chauveau, D., and Hunter, D. R. (2007), An EM-like algorithm
for semi- and non-parametric estimation in multivariate mixtures, Penn State
Department of Statistics Technical Report 07-01.
}
\examples{
## Examine and plot water-level task data set.

## First, try a 3-component solution where no two coordinates are
## assumed i.d.
data(Waterdata)
a <- npEMindrep(Waterdata, mu0=3, bw=4) # Assume indep but not iid
plot(a) # This produces 8 plots, one for each coordinate

## Next, same thing but pairing clock angles that are directly opposite one
## another (1:00 with 7:00, 2:00 with 8:00, etc.)
b <- npEMindrep(Waterdata, mu0=3, blockid=c(4,3,2,1,3,4,1,2), bw=4) # iid in pairs
plot(b) # Now only 4 plots, one for each block
}


\keyword{file}
