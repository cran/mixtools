normalmix.init = function (x, lambda = NULL, mu = NULL, s = NULL, k = 2, 
                           arbmean = TRUE, arbvar = TRUE) {
  n = length(x)
  x = sort(x)
  x.bin = list()
  for (j in 1:k) {
    x.bin[[j]] <- x[max(1, floor((j - 1) * n/k)):ceiling(j * n/k)]
  }
  if (is.null(s)) {
    s.hyp = as.vector(sapply(x.bin, sd))
    if (arbvar) {
      s = 1/rexp(k, rate = s.hyp)
    } else {
      s = 1/rexp(1, rate = mean(s.hyp))
    }
  }
  if (is.null(s) == FALSE && arbvar == TRUE) {
    k = length(s)
  }
  if (is.null(mu)) {
    mu.hyp <- as.vector(sapply(x.bin, mean))
	  if (arbmean) {
      mu = rnorm(k, mean = mu.hyp, sd = s)
	  } else {
      mu = rnorm(1, mean = mean(mu.hyp), sd = mean(s))
	  }
  }
  if (is.null(mu) == FALSE && arbmean == TRUE) {
    k = length(mu)
  }
  if (is.null(lambda)) {
    lambda = runif(k)
    lambda = lambda/sum(lambda)
  } else {
    k = length(lambda)
  }
  list(lambda = lambda, mu = mu, s = s, k = k)
}