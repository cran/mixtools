altnpEM <- function(x, mu0, blockid = 1:ncol(x),
                       bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                       h=bw, eps=1e-8,
                       maxiter=30, verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject
  u <- match(blockid, unique(blockid))
  if (is.matrix(mu0)) 
    m <- dim(mu0)[1]  # mu0=centers
  else 
    m <- mu0  # mu0=number of clusters
  if(!samebw && !is.matrix(bw)) { 
    bw <- matrix(bw, nrow=max(u), ncol=m)
  }
  z.hat <- matrix(0, nrow = n, ncol = m)
  tt0 <-  proc.time() # for total time
  ## Initial Values
  if(m == 1) z.hat <- matrix(1, nrow = n, ncol = m)
  else{
    kmeans <- kmeans(x, mu0)
    for(j in 1:m)
      z.hat[kmeans$cluster==j, j] <- 1
  }
  iter <- 0
  finished <- FALSE
  lambda <- matrix(0, nrow = maxiter, ncol = m)
  loglik <- NULL

  tmp <- 1:n
  ngrid <- 200
  xtra <- (max(x)-min(x))/10
  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  f <- array(1/m/diff(grid[1:2]), c(ngrid, m, r))
  Delta <- mean(diff(grid))
  oldloglik <- -Inf

  while (!finished) {
    iter <- iter + 1
    bw.old <- bw
    t0 <- proc.time()

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)

    ## density estimation step
    z=.C("altnpEM", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), as.double(bw),
       as.double(x), as.double(grid), old.f=as.double(f),
       new.f=as.double(f), as.double(lambda[iter,]), as.double(z.hat),
       conv=double(n*m*r))

#    newf <- f
#    for (k in 1:r) {
#      cat(" ", k, " ")
#      for (j in 1:m) {
#        for (a in 1:ngrid) {
#          for (i in 1:n) {
#            tmp [i] <- z.hat[i,j] * dnorm((x[i,k] - grid[a])/bw)/bw /
#            (sum(dnorm((x[i,k]-grid)/bw)/bw * f[,j,k]) * Delta)
#            #             cnvlv (x[i,k], grid, f[,j,k], bw, Delta)
#            #cnvlv <- function(x, u, fjk, h, delta) {
#              #  sum(dnorm((x-u)/h)/h * fjk) * delta
#            #}
#          }
#          newf[a,j,k] = f[a,j,k] / n / lambda[iter, j] * sum(tmp)
#          if (lambda[iter,j]==0 || sum(tmp)==0) browser()
#        }
#      }
#    }
#    f <- newf

    f <- array(z$new.f, c(ngrid, m, r))
#    print(apply(f,2:3,sum) * Delta)
#print(max(abs(f-f2)))
    
    plot(b, block=1, title=
         paste("Block 1, Alternative algorithm with bandwidth=", bw))
    for(i in 1:3) {
      lines(grid, f[,i,1]*lambda[iter,i], lwd=3, lty=2, col=i+1)
    }
    legend ("topleft", legend=round(lambda[iter,],3), col=2:4, lty=2, lwd=2)

    
#    browser()

    ## E-step (for next iteration)
    z=.C("altnpEM_Estep", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), as.double(bw),
       as.double(x), as.double(grid), f=as.double(f),
       as.double(lambda[iter,]), post=as.double(z.hat),
       loglik = double(1))

#    for (i in 1:n) {
#      tmp2 <- lambda[iter, ]
#      for (j in 1:m) {
#        for (k in 1:r) {
#          tmp2[j] <- tmp2[j] * cnvlv(x[i,k], grid, f[,j,k], bw, Delta)
#        }
#      }
#      z.hat[i, ] <- tmp2 / sum(tmp2)
#    }

    z.hat <- matrix(z$post, n, m)
#    print(max(abs(z.hat-z.hat2)))
    loglik <- z$loglik
    loglikchange <- loglik - oldloglik
    oldloglik <- loglik
    finished <- iter >= maxiter
    if (iter>1 && max(abs(lambda[iter, ] - lambda[iter-1, ])) < eps)
      finished <- TRUE
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, "  lambda ", round(lambda[iter, ], 4))
      cat(" obj change", round(loglikchange,4))
      cat(" time", (t1 - t0)[3], "\n")
    }
  }
  f <- array(z$f, c(ngrid, m, r))

  if (verb) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter, ], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }
  return(structure(list(data = x, posteriors = z.hat, 
                        lambda = lambda[1:iter,], bandwidth = bw, 
                        blockid = u, lambdahat = lambda[iter,], f=f,
                        grid=grid,
                        loglik = loglik),
                    class="npEM"))
}


cnvlv <- function(x, u, fjk, h, delta) {
  sum(dnorm((x-u)/h)/h * fjk) * delta
}


