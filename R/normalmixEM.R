normalmixEM <-
function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbmean = TRUE, arbvar = TRUE, 
          epsilon = 1e-08, maxit = 1000, maxrestarts=20, verb = FALSE, fast=FALSE) {
  x <- as.vector(x)
  tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
                        k = k, arbmean = arbmean, arbvar = arbvar)
  lambda <- tmp$lambda 
  mu <- tmp$mu 
  sigma <- tmp$s
  k <- tmp$k 
  arbvar <- tmp$arbvar 
  arbmean <- tmp$arbmean
  if (fast==TRUE && k==2 && arbmean==TRUE) {
    a <- normalmixEM2comp (x, lambda=lambda[1], mu=mu, sigsqrd=sigma^2, 
                           eps=epsilon, maxit=maxit, verb=verb)
  } else {
    n <- length(x)
    const <- n*log(2*pi)/2
    notdone <- TRUE
    while(notdone) {
      # Initialize everything
      notdone <- FALSE
      tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
                            k = k, arbmean = arbmean, arbvar = arbvar)
      lambda <- tmp$lambda
      mu <- tmp$mu
      k <- tmp$k
      sigma <- tmp$s
      var <- sigma^2
      diff <- epsilon+1
      iter <- 0
      postprobs <- matrix(nrow = n, ncol = k)
      restarts <- 0
      mu <- rep(mu, k)[1:k]
      sigma <- rep(sigma,k)[1:k]
      res <- outer(x, mu, "-")^2
      normres <- sweep(res, 2, 2*sigma^2, "/")
      dnorms <- sweep(exp(-sweep(res, 2, 2*sigma^2, "/")), 2, lambda/sigma, "*")
      ll <- obsloglik <- sum(log(apply(dnorms, 1, sum))) - const
      while (diff > epsilon && iter < maxit) {
        # Main EM loop.  First, calculate posteriors (E-step)
        tmp <- apply(dnorms, 1, sum)
        if(any(tmp==0)) { 
          if (verb) {
            cat(paste("warning:  Some posteriors may be 0/0 in iteration ",iter,".\n",
                      "Using slower but more stable calculation.\n")) 
          }
          ratio <- lambda/sigma
          for(j in 1:(k-1)) {
            postprobs[,j] = 1/(apply(sweep(exp(as.matrix(normres[,j]-normres[,-j])), 
                                           2, ratio[-j], "*"), 1, sum)/ratio[j] + 1)
          }
          postprobs[,k]=1-apply(as.matrix(postprobs[,(1:(k-1))]),1,sum)
        } else {
          postprobs <- sweep(dnorms, 1, tmp, "/")
        }
        # Next, update parameters (M-step)
        lambda.new <- apply(postprobs, 2, mean)
        sing <- any(is.nan(postprobs)) | any(is.na(lambda.new)) | any(lambda.new<1e-08)
        if(!sing) {
          if (arbmean) {
            mu.new <- apply(postprobs * x, 2, mean)/lambda.new
          } else {
            mu.new <- sum(postprobs*x)/n
          }
          mu.new <- rep(mu.new,k)[1:k]      
          if (arbvar) {
            sigma2.new <- apply(postprobs*res,2,mean)/lambda.new
          } else {
            sigma2.new <- sum(postprobs*res)/n
          }
          res <- outer(x, mu.new, "-")^2
          lambda <- lambda.new
          mu <- mu.new
          sigma <- sqrt(sigma2.new)
          sing <- sum(sigma < 1e-08)
          sigma <- rep(sigma, k)[1:k]
          normres <- sweep(res, 2, 2*sigma^2, "/")
          dnorms <- sweep(exp(-sweep(res, 2, 2*sigma^2, "/")), 2, lambda/sigma, "*")                                                       
          newobsloglik <- sum(log(apply(dnorms, 1, sum))) - const
        }
        sing <- sing | (abs(newobsloglik)==Inf) | is.na(newobsloglik)# | sum(postprobs) != n
        if (sing) {
          notdone <- TRUE
          cat("Need new starting values due to singularity...\n")
          restarts <- restarts + 1
          lambda <- mu <- sigma <- NULL
          if(restarts>maxrestarts) stop("Too many tries!")
        } else {
          diff <- newobsloglik - obsloglik
          obsloglik <- newobsloglik
          ll <- c(ll, obsloglik)
          iter <- iter + 1
          if (verb) {
            cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", 
                obsloglik, "\n")
          }
        }
      }
    }
    if (iter == maxit) {
      cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if(arbmean == FALSE){
      scale.order = order(sigma)
      sigma.min = min(sigma)
      postprobs = postprobs[,scale.order]
      colnames(postprobs) <- c(paste("comp", ".", 1:k, sep = ""))
      a=list(x=x, lambda = lambda[scale.order], mu = mu, sigma = sigma.min, 
             scale = sigma[scale.order]/sigma.min, loglik = obsloglik, 
             posterior = postprobs, all.loglik=ll, restarts=restarts, 
             ft="normalmixEM")
    } else {
      colnames(postprobs) <- c(paste("comp", ".", 1:k, sep = ""))
      a=list(x=x, lambda = lambda, mu = mu, sigma = sigma, loglik = obsloglik, 
             posterior = postprobs, all.loglik=ll, restarts=restarts, 
             ft="normalmixEM")
    }
  }
  class(a) = "mixEM"
  a
}

