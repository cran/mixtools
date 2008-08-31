multmixEM <- function (y, lambda = NULL, theta = NULL, k = 2, maxit = 10000, 
    epsilon = 1e-08, verb = FALSE) {
  n <- nrow(y)
  p <- ncol(y)
  m <- colSums(y)
  r <- rowSums(y)
  if (any(r!=r[1])) stop("All row sums of y must be equal")  
  mncoeffs <- lgamma(1+r[1]) - rowSums(lgamma(1+y))

  tmp <- multmix.init(y=y, lambda=lambda, theta=theta, k=k)
  lambda <- tmp$lambda
  theta <- tmp$theta
  k <- tmp$k
  restarts<-0
  mustrestart <- FALSE

  while (restarts < 15) {  
    oldll <- -Inf  
    num <- exp(y %*% t(log(pmax(theta,1e-100))) + outer(mncoeffs,log(lambda),"+"))
    dens <- rowSums(num)
    postz <- num/dens
    newll <- sum(log(dens))
    diff <- newll - oldll
    ll <- NULL
    iter <- 0
    while ((iter < maxit) && diff > epsilon) {
      iter <- iter + 1
      oldll <- newll
      ll <- c(ll, oldll)
      zty <- t(postz) %*% y
      theta <- zty/apply(zty, 1, sum)
      theta <- cbind(theta[,1:(p-1)],1-apply(as.matrix(theta[,1:(p-1)]),1,sum))
      lambda <- apply(postz, 2, mean)
      
      num <- exp(y %*% t(log(pmax(theta,1e-100))) + outer(mncoeffs,log(lambda),"+"))
      dens <- rowSums(num)
      postz <- num/dens
      newll <- sum(log(dens))
      diff <- newll - oldll      
      if (diff<0 || is.na(newll) || is.infinite(newll) || is.nan(newll)) {
        mustrestart <- TRUE
        break
      }      
      if (verb) {
        cat("iteration=", iter, "diff=", diff, "log-likelihood", 
            ll[iter], "lambda", lambda, "\n") 
      }
    }
    if (mustrestart) {
      cat("Restarting due to numerical problem.\n")
      mustrestart <- FALSE
      restarts <- restarts + 1
      tmp <- multmix.init(y=y, k=k)
      lambda <- tmp$lambda
      theta <- tmp$theta
      k <- tmp$k
    } else {
      if (iter == maxit) {
        cat("Warning: Not convergent after", maxit, "iterations\n")
      }
      theta[,p] <- 1-apply(as.matrix(theta[,1:(p-1)]),1,sum)
      colnames(theta) <- c(paste("theta", ".", 1:p, sep = ""))
      rownames(theta) <- c(paste("comp", ".", 1:k, sep = ""))
      colnames(postz) <- c(paste("comp", ".", 1:k, sep = ""))
      cat("number of iterations=", iter, "\n")
      return(list(y=y, lambda = lambda, 
                  theta = theta, loglik = sum(log(dens)), posterior = postz, 
                  all.loglik=ll, restarts=restarts, ft="multmixEM"))
    }
  }
  stop("Too many restarts!")
}
