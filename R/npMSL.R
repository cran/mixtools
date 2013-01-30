## nonparametric algorithm for Smoothed Likelihood Maximization 
## implementing block structure
#npMSL.old <- function(x, mu0, blockid = 1:ncol(x),
#                       bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
#                       h=bw, eps=1e-8,
#                       maxiter=500, ngrid=200, verb = TRUE){
#  bw <- h # h is alternative bandwidth argument, for backward compatibility
#  x <- as.matrix(x)
#  n <- nrow(x)      # number of subjects
#  r <- ncol(x)      # number of measurements in each subject

#  u <- match(blockid, unique(blockid))
#  B <- max(u) 		# nb of blocks
#  BlS <- rep(0,B)	# block sizes = C_ell in JCGS paper
#  for (ell in 1:B) {
#      BlS[ell] <- sum(u == ell)}

#  if (is.matrix(mu0)) # mu0=centers
#    m <- dim(mu0)[1]  else m <- mu0  # mu0=number of clusters
    
#  if(!samebw && !is.matrix(bw)) {
#    bw <- matrix(bw, nrow=max(u), ncol=m)
#  }
#  z.hat <- matrix(0, nrow = n, ncol = m)
#  tt0 <-  proc.time() # for total time
  ## Initial Values
#  if(m == 1) z.hat <- matrix(1, nrow = n, ncol = m) else {
#    kmeans <- kmeans(x, mu0)
#    for(j in 1:m)
#      z.hat[kmeans$cluster==j, j] <- 1
#  }
#  iter <- 0
#  finished <- FALSE
#  lambda <- matrix(0, nrow = maxiter, ncol = m)
#  loglik <- NULL
#  total_udfl <- 0; total_nan <- 0 # for internal checks

#  tmp <- 1:n
#  xtra <- (max(x)-min(x))/10
#  grid <- seq(min(x)-xtra, max(x)+xtra, length=ngrid)
  # f stored on a ngrid by m by B array 
  # f_{g,j,ell} = f_{j ell}(u_g)
  # f <- array(1/m/diff(grid[1:2]), c(ngrid, m, B)) 
  # this f was not normalized for being uniform over grid 
  	
#  Delta <- diff(grid[1:2])
  
#  f <- array(1/((ngrid-1)*Delta), c(ngrid, m, B))
#  oldloglik <- -Inf

#  while (!finished) {
#    iter <- iter + 1
#    bw.old <- bw
#    t0 <- proc.time()
#    nb_udfl=0;  # nb of underflows log(0) cancelled in nems_Estep_v3.c
#    nb_nan=0;  # nb of nonzero K()*log(0) cancelled in nems_Estep_v3.c

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
#    lambda[iter, ] <- colMeans(z.hat)
    ## density estimation M-step   
#    z=.C("npMSL_Mstep", as.integer(ngrid), as.integer(n),
#       as.integer(m), as.integer(r),
#       as.integer(B), as.integer(BlS), as.integer(u),
#       as.double(bw), as.double(x), as.double(grid),
#       new.f=as.double(f),
#       as.double(lambda[iter,]),
#       as.double(z.hat),
#       PACKAGE="mixtools")
       
#    f <- array(z$new.f, c(ngrid, m, B)) # check sum(f == 0)
## print(apply(f,2:3,sum) * Delta)
## print(max(abs(f-f2)))
##    browser()

    ## E-step (for next iteration)
#    z=.C("npMSL_Estep", as.integer(ngrid), as.integer(n),
#       as.integer(m), as.integer(r),
#       as.integer(B), as.integer(u),
#       as.double(bw),
#       as.double(x), as.double(grid), f=as.double(f),
#       as.double(lambda[iter,]), post=as.double(z.hat),
#       loglik = double(1),
#       nb_udfl = as.integer(nb_udfl), nb_nan = as.integer(nb_nan),
#       PACKAGE="mixtools")

#    nb_udfl = z$nb_udfl; nb_nan = z$nb_nan;
#    total_udfl <- total_udfl + nb_udfl
#    total_nan <- total_nan + nb_nan
    
#    z.hat <- matrix(z$post, n, m)
#    if (sum(is.nan(z.hat)) > 0) cat("Error!! NaN in z.hat") # obsolete ?
#    loglik <- z$loglik
#    loglikchange <- loglik - oldloglik
#    oldloglik <- loglik
#    finished <- iter >= maxiter
#    if (iter>1 && max(abs(lambda[iter, ] - lambda[iter-1, ])) < eps)
#      finished <- TRUE
#    if (verb) {
#      t1 <- proc.time()
#      cat("iteration", iter, ": lambda ", round(lambda[iter, ], 4))
#      cat(" obj change", round(loglikchange,4))
#      cat(" time", (t1 - t0)[3])
#      if (nb_udfl > 0) cat("average underflows=", nb_udfl/(n*m*r)," ")
#      if (nb_nan >0) cat("average NaNs=", nb_nan/(n*m*r))
# Note: these average mean nb of nan over ngrid convolution
#      cat("\n")
#    }
#  }
## f <- array(z$f, c(ngrid, m, r)) # obsolete in block version

#  if (verb) {
#    tt1 <- proc.time()
#    cat("lambda ", round(lambda[iter, ], 4))
#    cat(", total time", (tt1 - tt0)[3], "s\n")
#  }
#  return(structure(list(data = x, posteriors = z.hat,
#                        lambda = lambda[1:iter,], bandwidth = bw,
#                        blockid = u, lambdahat = lambda[iter,], f=f,
#                        grid = grid, loglik = loglik,
#			meanUdfl = total_udfl/(n*m*r*iter),# average underflow
#			meanNaN = total_nan/(n*m*r*iter)), # average nan
#                    class="npEM")) # define a "NEMS" class ?
#}






################################################################ 
################################################################ 
# development of version with adaptive bandwidth
################################################################ 
################################################################ 
npMSL <- function(x, mu0, blockid = 1:ncol(x),
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
    
  if(!samebw && !is.matrix(bw)) { # create initial bandwidth matrix 
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
  total_udfl <- 0; total_nan <- 0 # eventual NaN and underflow in C code

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

  orderx <- xx <- list() # preparation for adaptive bandwidth
  for(k in 1:B) {
    xx[[k]] <- as.vector(x[, u==k])
    if (!samebw) {
      orderx[[k]] = order(xx[[k]]) # only needed for IQR calculation for bw
    }
  }
  CftEstep <- ifelse(samebw, "npMSL_Estep", "npMSL_Estep_bw")
  # CftEstep <- "npMSL_Estep_bw" # temporary, for testing only the M-step  
  CftMstep <- ifelse(samebw, "npMSL_Mstep", "npMSL_Mstep_bw")
  
  while (!finished) {
    iter <- iter + 1
    bw.old <- bw	# is this needed?
    t0 <- proc.time()
    nb_udfl=0;  # nb of underflows log(0) cancelled in nems_Estep.c
    nb_nan=0;  # nb of nonzero K()*log(0) cancelled in nems_Estep.c

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)
    
    ## density estimation in M-step: WKDE-step  
      cs <- colSums(z.hat)
      z.tmp <- sweep(z.hat, 2, cs, "/")
      z.tmp[, cs==0] <- 1/NROW(z.tmp) # Just in case

    ## adaptive bandwidth update
    for (ell in 1:B) {
      r2 <- BlS[ell]	# block size = nb of coordinates
      if (!samebw) {
        wts <- apply(z.tmp, 2, function(z) rep(z/r2, r2))
        variances <- colSums(wts*outer(xx[[ell]], colSums(wts*xx[[ell]]),
        			'-')^2)
        iqr <- apply(as.matrix(wts[orderx[[ell]],]), 2, wIQR, 
			xx[[ell]][orderx[[ell]]],
			already.sorted=TRUE, already.normalized=TRUE)
        h <- bw[ell, ] <- 0.9 * pmin(sqrt(variances), iqr/1.34) * 
                   pmax(1,r2*n*lambda[iter, ])^(-1/5) 
                   # Note:  Doesn't allow "sample size" < 1.
#      browser()
      }
    }
        
    z=.C(CftMstep, as.integer(ngrid), as.integer(n),
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
    z=.C(CftEstep, as.integer(ngrid), as.integer(n),
       as.integer(m), as.integer(r), 
       as.integer(B), as.integer(u),
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
      if ((nb_udfl > 0) || (nb_nan >0)) cat("\n  ")
      if (nb_udfl > 0) cat("average underflows=", round(nb_udfl/(n*m*r),3)," ")
      if (nb_nan >0) cat("average NaNs=", round(nb_nan/(n*m*r),3))
# Note: these average mean nb of nan over ngrid convolution
      cat("\n")
    }
  }
  # f <- array(z$f, c(ngrid, m, r)) # obsolete in block version

 if (!samebw) {
    rownames(bw) <- paste("block", 1:max(u))
    colnames(bw) <- paste("component", 1:m)
  	}

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
			meanNaN = total_nan/(n*m*r*iter)), # average NaN's
                    class="npEM")) # define a "npMSL" class ?
}

