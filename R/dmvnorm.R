# Alternative version of dmvnorm to eliminate dependence of mixtools
# on additional package 'mvtnorm'
# Written (hopefully) to be more efficient than mvtnorm version, which uses both
# a call to "eigen" and a call to "mahalanobis", by using only a single
# call to the more efficient "qr" (read "Note" under ?qr)

# Note:  These functions assume that each ROW of y is a separate position vector.
# i.e., y is assumed to be nxd, where d=dimension
dmvnorm <- function(y, mean=NULL, sigma=NULL) {
  exp(logdmvnorm(y, mean=mean, sigma=sigma))
}

logdmvnorm <- function(y, mean=NULL, sigma=NULL) {
  if (is.vector(y)) 
    y <- matrix(y, nrow=1)
  d <- ncol(y)
  if (is.null(mean))
    mean <- rep(0, d)
  if (is.null(sigma))
    sigma <- diag(d)
  y <- sweep(y, 2, mean)
  k <- d * 1.8378770664093454836 # that constant is log(2*pi)
  a <- qr(sigma)
  logdet <- sum(log(abs(diag(a$qr))))
  if(nrow(y)==1) 
    mahaldist <- as.vector(y %*% qr.solve(a,t(y)))
  else
    mahaldist <- rowSums((y %*% qr.solve(a)) * y)
  -0.5*(mahaldist + logdet + k)
}


