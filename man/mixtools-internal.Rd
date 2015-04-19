\name{mixtools-internal}
\alias{kern.B}
\alias{kern.C}
\alias{kern.G}
\alias{kern.O}
\alias{kern.T}
\alias{kfoldCV}
\alias{npMSL_old}
\alias{splitsample}
\alias{wbw.kCV}

\title{Internal 'mixtools' Functions}
\description{
  Internal kernel and semiparametric-related functions for the package \code{mixtools}.
}
\usage{
kern.B(x, xi, h, g = 0)
kern.C(x, xi, h)
kern.G(x, xi, h)
kern.O(x, xi, h)
kern.T(x, xi, h)
kfoldCV(h, x, nbsets = 2, w = rep(1, length(x)), 
        lower = mean(x) - 5*sd(x), upper = mean(x) + 5*sd(x))
npMSL_old(x, mu0, blockid = 1:ncol(x),
          bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
          h=bw, eps=1e-8, maxiter=500, bwiter = maxiter,
          ngrid = 200, post = NULL, verb = TRUE)
splitsample(n, nbsets = 2)
wbw.kCV(x, nbfold = 5, w = rep(1, length(x)), 
        hmin = 0.1*hmax, hmax = NULL)
}

\arguments{
  \item{x}{A vector of values to which local modeling techniques are applied.}
  \item{xi}{An n-vector of data values.}
  \item{h}{The bandwidth controlling the size of the window used for the
  local estimation around \code{x}.}
  \item{g}{A shape parameter required for the symmetric beta kernel.  The default
  is \code{g} = 0 which yields the uniform kernel.  Some common values are \code{g} = 1 for the
  Epanechnikov kernel, \code{g} = 2 for the biweight kernel, and \code{g} = 3 for the triweight kernel.}
  \item{mu0}{See updated arguments in the \code{npMSL} function.}
  \item{blockid}{See updated arguments in the \code{npMSL} function.}
  \item{bw}{See updated arguments in the \code{npMSL} function.}
  \item{samebw}{See updated arguments in the \code{npMSL} function.}
  \item{h}{See updated arguments in the \code{npMSL} function.}
  \item{eps}{See updated arguments in the \code{npMSL} function.}
  \item{maxiter}{See updated arguments in the \code{npMSL} function.}
  \item{bwiter}{See updated arguments in the \code{npMSL} function.}
  \item{ngrid}{See updated arguments in the \code{npMSL} function.}
  \item{post}{See updated arguments in the \code{npMSL} function.}
  \item{verb}{See updated arguments in the \code{npMSL} function.}
  \item{n}{See updated arguments in the \code{npMSL} function.}
  \item{nbsets}{See updated arguments in the \code{npMSL} function.}
  \item{w}{See updated arguments in the \code{npMSL} function.}
  \item{lower}{See updated arguments in the \code{npMSL} function.}
  \item{upper}{See updated arguments in the \code{npMSL} function.}
  \item{nbfold}{See updated arguments in the \code{npMSL} function.}
  \item{hmin}{See updated arguments in the \code{npMSL} function.}
  \item{hmax}{See updated arguments in the \code{npMSL} function.}
}

\details{
  These are usually not to be called by the user.
}

\seealso{
  \code{\link{npMSL}}
}


\keyword{internal}
