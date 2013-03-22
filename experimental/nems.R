## successive versions of nems algorithms under development
## nems.v1 : initial version
## nems.v2 : nems with underflow handling in the E-step
## current = nems.v3 : v2 plus block structure

nems.v1 <- function(x, mu0, blockid = 1:ncol(x),
                       bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                       h=bw, eps=1e-8,
                       maxiter=500, ngrid=200, verb = TRUE){
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
#  ngrid <- 200
  xtra <- (max(x)-min(x))/10
  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  f <- array(1/m/diff(grid[1:2]), c(ngrid, m, r)) 
  # nb: = 1/(m*diff(grid[1:2]))
  # f not normalized for uniform over grid ?	
  Delta <- mean(diff(grid)) # = diff(grid[1:2]) ?
  oldloglik <- -Inf

  while (!finished) {
    iter <- iter + 1
    bw.old <- bw
    t0 <- proc.time()

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)

    ## density estimation step
    
    z=.C("nems_Mstep_v1", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), as.double(bw),
       as.double(x), as.double(grid), #old.f=as.double(f),
       new.f=as.double(f), 
       #as.double(lambda[iter,]), 
       as.double(z.hat)#,
       #conv=double(n*m*r)
       )
    f <- array(z$new.f, c(ngrid, m, r))
#    print(apply(f,2:3,sum) * Delta)
# print(max(abs(f-f2)))
    
#    plot(b, block=1, title=
#         paste("Block 1, NEMS-like algorithm, bandwidth=", bw, ", iter=",iter))
#    for(i in 1:3) {
#     lines(grid, f[,i,1]*lambda[iter,i], lwd=3, lty=2, col=i+1)
#   }
#    legend ("topleft", legend=round(lambda[iter,],3), col=2:4, lty=2, lwd=2)

    
#    browser()

    ## E-step (for next iteration)
    z=.C("nems_Estep_v1", as.integer(ngrid), as.integer(n),
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

####################################################
## version 2 with underflow checking
nems.v2 <- function(x, mu0, blockid = 1:ncol(x),
                       bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                       h=bw, eps=1e-8,
                       maxiter=500, ngrid=200, verb = TRUE){
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
  xtra <- (max(x)-min(x))/10
  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  f <- array(1/m/diff(grid[1:2]), c(ngrid, m, r)) 
  # nb: = 1/(m*diff(grid[1:2]))
  # f not normalized for uniform over grid ?	
  Delta <- mean(diff(grid)) # = diff(grid[1:2]) ?
  oldloglik <- -Inf

  while (!finished) {
    iter <- iter + 1
    bw.old <- bw
    t0 <- proc.time()
	nb_udfl=0;	# nb of underflows log(0) cancelled in nems_Estep_v2.c
	nb_nan=0;	# nb of nonzero K()*log(0) cancelled in nems_Estep_v2.c

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)

    ## density estimation step - v2 
    
    z=.C("nems_Mstep_v2", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), as.double(bw),
       as.double(x), as.double(grid), #old.f=as.double(f),
       new.f=as.double(f), 
       #as.double(lambda[iter,]), 
       as.double(z.hat))
    f <- array(z$new.f, c(ngrid, m, r)) # check sum(f == 0)
# print(apply(f,2:3,sum) * Delta)
# print(max(abs(f-f2)))
    
#    browser()

    ## E-step (for next iteration)
    z=.C("nems_Estep_v2", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), as.double(bw),
       as.double(x), as.double(grid), f=as.double(f),
       as.double(lambda[iter,]), post=as.double(z.hat),
       loglik = double(1),
       nb_udfl=as.integer(nb_udfl), nb_nan=as.integer(nb_nan))

	nb_udfl = z$nb_udfl; nb_nan = z$nb_nan; 
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
    if (sum(is.nan(z.hat)) > 0) cat("Error!! NaN in z.hat")
#    print(max(abs(z.hat-z.hat2)))
    loglik <- z$loglik
    loglikchange <- loglik - oldloglik
    oldloglik <- loglik
    finished <- iter >= maxiter
    if (iter>1 && max(abs(lambda[iter, ] - lambda[iter-1, ])) < eps)
      finished <- TRUE
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, ": lambda ", round(lambda[iter, ], 4))
      cat(" obj change", round(loglikchange,4))
      cat(" time", (t1 - t0)[3], "\n")
      if (nb_udfl > 0) cat("average underflows=", nb_udfl/(n*m*r)," ")
      if (nb_nan > 0) cat("average NaNs=", nb_nan/(n*m*r),"\n")
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



#######################################################
## current version = version 3 
## implementing block structure
nems <- nems.v3 <- function(x, mu0, blockid = 1:ncol(x),
                       bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                       h=bw, eps=1e-8,
                       maxiter=500, ngrid=200, verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject

  u <- match(blockid, unique(blockid))
  B <- max(u) 		# nb of blocks
  BlS <- rep(0,B)	# block sizes = C_ell in JCGS paper
  for (ell in 1:B) {
      BlS[ell] <- sum(u == ell)}      
      
  if (is.matrix(mu0)) # mu0=centers
    m <- dim(mu0)[1]  else m <- mu0  # mu0=number of clusters
    
  if(!samebw && !is.matrix(bw)) { 
    bw <- matrix(bw, nrow=max(u), ncol=m)  
  }
  z.hat <- matrix(0, nrow = n, ncol = m)
  tt0 <-  proc.time() # for total time
  ## Initial Values
  if(m == 1) z.hat <- matrix(1, nrow = n, ncol = m) else {
    kmeans <- kmeans(x, mu0)
    for(j in 1:m)
      z.hat[kmeans$cluster==j, j] <- 1
  }
  iter <- 0
  finished <- FALSE
  lambda <- matrix(0, nrow = maxiter, ncol = m)
  loglik <- NULL
  total_udfl <- 0; total_nan <- 0 # for keeping track and output result

  tmp <- 1:n
  xtra <- (max(x)-min(x))/10
  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  # f stored on a ngrid by m by B array 
  # f_{g,j,ell} = f_{j ell}(u_g)
  # f <- array(1/m/diff(grid[1:2]), c(ngrid, m, B)) 
  # this f was not normalized for being uniform over grid 
  	
  Delta <- diff(grid[1:2]) 
  
  f <- array(1/((ngrid-1)*Delta), c(ngrid, m, B)) 
  oldloglik <- -Inf

  while (!finished) {
    iter <- iter + 1
    bw.old <- bw
    t0 <- proc.time()
    nb_udfl=0;  # nb of underflows log(0) cancelled in nems_Estep_v3.c
    nb_nan=0;  # nb of nonzero K()*log(0) cancelled in nems_Estep_v3.c

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)
    ## density estimation M-step   
    z=.C("nems_Mstep_v3", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), 
       as.integer(B), as.integer(BlS), as.integer(u),
       as.double(bw), as.double(x), as.double(grid), 
       new.f=as.double(f), 
       as.double(lambda[iter,]), 
       as.double(z.hat))
       
    f <- array(z$new.f, c(ngrid, m, B)) # check sum(f == 0)
# print(apply(f,2:3,sum) * Delta)
# print(max(abs(f-f2)))
    
#    browser()

    ## E-step (for next iteration)
    z=.C("nems_Estep_v3", as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), as.integer(u),
       as.double(bw),
       as.double(x), as.double(grid), f=as.double(f),
       as.double(lambda[iter,]), post=as.double(z.hat),
       loglik = double(1),
       nb_udfl = as.integer(nb_udfl), nb_nan = as.integer(nb_nan))

    nb_udfl = z$nb_udfl; nb_nan = z$nb_nan; 
    total_udfl <- total_udfl + nb_udfl
    total_nan <- total_nan + nb_nan
    
    z.hat <- matrix(z$post, n, m)
    if (sum(is.nan(z.hat)) > 0) cat("Error!! NaN in z.hat") # obsolete ?
    loglik <- z$loglik
    loglikchange <- loglik - oldloglik
    oldloglik <- loglik
    finished <- iter >= maxiter
    if (iter>1 && max(abs(lambda[iter, ] - lambda[iter-1, ])) < eps)
      finished <- TRUE
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, ": lambda ", round(lambda[iter, ], 4))
      cat(" obj change", round(loglikchange,4))
      cat(" time", (t1 - t0)[3])
      if (nb_udfl > 0) cat("average underflows=", nb_udfl/(n*m*r)," ")
      if (nb_nan >0) cat("average NaNs=", nb_nan/(n*m*r))
# Note: these average mean nb of nan over ngrid convolution
      cat("\n")
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
                        grid = grid, loglik = loglik,
			meanUdfl = total_udfl/(n*m*r*iter),# average underflow
			meanNaN = total_nan/(n*m*r*iter)), # average nan
                    class="npEM")) # define a "NEMS" class ?
}

