\name{mixtools initializations}
\alias{logisregmix.init}
\alias{multmix.init}
\alias{mvnormalmix.init}
\alias{normalmix.init}
\alias{poisregmix.init}
\alias{regmix.init}
\alias{regmix.mixed.init}
\alias{repnormmix.init}

\title{Initializations for Various EM Algorithms in 'mixtools'}
\description{
  Internal intialization functions for EM algorithms in the package \code{mixtools}.
}
\usage{
logisregmix.init(y, x, N, lambda = NULL, beta = NULL, k = 2)
multmix.init(y, lambda = NULL, theta = NULL, k = 2)
mvnormalmix.init(x, lambda = NULL, mu = NULL, sigma = NULL, 
                 k = 2, arbmean = TRUE, arbvar = TRUE)
normalmix.init(x, lambda = NULL, mu = NULL, s = NULL, k = 2, 
               arbmean = TRUE, arbvar = TRUE)
poisregmix.init(y, x, lambda = NULL, beta = NULL, k = 2)
regmix.init(y, x, lambda = NULL, beta = NULL, s = NULL, k = 2, 
            addintercept = TRUE, arbmean = TRUE, arbvar=TRUE)
regmix.mixed.init(y, x, w = NULL, sigma = NULL, 
                  arb.sigma = TRUE, alpha = NULL, lambda = NULL, 
                  mu = NULL, R = NULL, arb.R = TRUE, k = 2, 
                  mixed = FALSE, addintercept.fixed = FALSE,
                  addintercept.random = TRUE)
repnormmix.init(x, lambda = NULL, mu = NULL, s = NULL, k = 2, 
                arbmean = TRUE, arbvar = TRUE)
}

\details{
  These are usually not to be called by the user. Definitions of the arguments appear in the respective EM algorithms.
}

\keyword{internal}
