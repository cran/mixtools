\name{mixtools-internal}
\alias{kern.B}
\alias{kern.C}
\alias{kern.G}
\alias{kern.O}
\alias{kern.T}

\title{Kernels for 'mixtools' Functions}
\description{
  Internal kernel functions for the package \code{mixtools}.
}
\usage{
kern.B(x, xi, h, g = 0)
kern.C(x, xi, h)
kern.G(x, xi, h)
kern.O(x, xi, h)
kern.T(x, xi, h)
}

\arguments{
  \item{x}{A vector of values to which local modeling techniques are applied.}
  \item{xi}{An n-vector of data values.}
  \item{h}{The bandwidth controlling the size of the window used for the
  local estimation around \code{x}.}
  \item{g}{A shape parameter required for the symmetric beta kernel.  The default
  is \code{g} = 0 which yields the uniform kernel.  Some common values are \code{g} = 1 for the
  Epanechnikov kernel, \code{g} = 2 for the biweight kernel, and \code{g} = 3 for the triweight kernel.}
}

\details{
  These are usually not to be called by the user.
}

\keyword{internal}
