## EM-like algorithm for a nonparametric mixture model with
## independent repeated measures - some ID, some not
npEMindrep <- function(x, mu0, blockid = 1:ncol(x),
                       bw=bw.nrd0(as.matrix(x)), h=bw, eps=1e-8,
                       maxiter=100, stochastic = FALSE,
                       verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject
  if (is.matrix(mu0)) m <- dim(mu0)[1] # mu0=centers
  	else m <- mu0 # when mu0=number of clusters
  z.hat <- matrix(0, nrow=n, ncol=m)
  u <- match(blockid, unique(blockid)) 
  tt0 <-  proc.time() # for total time
  ## Initial Values
  kmeans <- kmeans(x, mu0)
  for(j in 1:m)
    z.hat[kmeans$cluster==j, j] <- 1

  iter <- 0
  if (stochastic) {
    sumpost <- matrix(0, n, m)
  }
  finished <- FALSE
  lambda <- matrix(0,maxiter,m)
  while(!finished) {
    iter <- iter + 1
    t0 <- proc.time()

    ## M-Step
    lambda[iter,] <- colMeans(z.hat) 

    ## E-step
    if(stochastic){
      z <- t(apply(z.hat, 1, function(prob) rmultinom(1, 1, prob)))
      z.tmp <- sweep(z, 2, colSums(z), "/")
    }
    else { 
      z.tmp <- sweep(z.hat, 2, colSums(z.hat), "/")
    }
    fkernel <- matrix(1, n, m)
    for (k in 1:max(u)) {
      r2 = sum(u==k)
      ans <- .C("kernelsmoothrepeated", n=as.integer(n),
                m=as.integer(m), r=as.integer(r2), x=as.double(x[,u==k]),
                bw=as.double(bw), z=as.double(z.tmp), f=double(n*m), 
              PACKAGE="mixtools")
      fkernel <- fkernel * matrix(ans$f, ncol=m)
    }
    lambda.f <- sweep(fkernel, 2, lambda[iter,], "*")
    z.hat <- lambda.f/rowSums(lambda.f)

    finished <- iter >= maxiter
    if (stochastic) {
      sumpost <- sumpost + z.hat
    } else if (iter>1) { # This convergence criterion is too simplistic:
      change <- lambda[iter,] - lambda[iter-1,]
      finished <- finished | (max(abs(change)) < eps)
    }
    if (verb) {
    t1 <- proc.time()
    cat("iteration", iter, "  lambda ", round(lambda[iter,], 4))
    cat(" time", (t1-t0)[3], "\n")
    }
  }
  if (verb) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter,], 4))
    cat(", total time", (tt1-tt0)[3], "s\n")    
  }
  if(stochastic) {
    return(structure(list(data=x, posteriors=sumpost/iter, lambda=lambda,
                          bandwidth=bw, blockid=u, lambdahat=colMeans(lambda)), 
                        class="npEM"))
  } else {
    return(structure(list(data=x, posteriors=z.hat, lambda=lambda[1:iter,],
                          bandwidth=bw, blockid=u, lambdahat=lambda[iter,]), 
                        class="npEM"))
  }
}



