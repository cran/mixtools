multmixEM = function (y, lambda = NULL, theta = NULL, k = 2, maxit = 10000, 
    epsilon = 1e-08, verb = FALSE) 
{
    n <- nrow(y)
    p <- ncol(y)
    m <- apply(y, 2, sum)
    tmp <- multmix.init(y=y, lambda=lambda, theta=theta, k=k)
    lambda <- tmp$lambda
    theta <- tmp$theta
    k <- tmp$k
        if(is.null(theta) == FALSE && sum(apply(theta,1,sum)!=1)>0){
        cat("Rows of theta must sum to 1!","\n")
} 


    postz = matrix(0, n, k)
    loglik <- numeric(maxit)
    oldll = -Inf
    iter <- 0
    diff <- epsilon * 2
    ll <- NULL
	restarts<-0
    while ((iter < maxit) && diff > epsilon) {
         iter <- iter + 1
        for (j in 1:k) {
            postz[, j] = lambda[j] * exp(apply(y, 1, function(x, 
                th) ldmult(x, th), th = theta[j, ]))
        }

            dens = apply(postz, 1, sum)
            newll = loglik[iter] = sum(log(dens))
            postz <- postz/dens
#        postz <- cbind(postz[,1:(k-1)],1-apply(as.matrix(postz[,1:(k-1)]),1,sum))

    if(sum(is.nan(postz))>0){
        cat("Need new starting values due to underflow...","\n")
		restarts<-restarts+1
		if(restarts>15) stop("Too many tries!")
            tmp <- multmix.init(y=y, k=k)
            lambda <- tmp$lambda
        theta <- tmp$theta
            k <- tmp$k
            loglik <- numeric(maxit)
            oldll = -Inf
            iter <- 0
            diff <- epsilon * 2
        ll <- NULL
        } else {
        diff = newll - oldll
            oldll = newll
        ll <- c(ll, oldll)
        zty = t(postz) %*% y
        newtheta = zty/apply(zty, 1, sum)
      newtheta = cbind(newtheta[,1:(p-1)],1-apply(as.matrix(newtheta[,1:(p-1)]),1,sum))
        newlambda = apply(postz, 2, mean)
        theta <- newtheta
        lambda <- newlambda
    sing <- sum(lambda<1e-08)*(k>1)
            if(diff<0 || is.na(newll) || sing>0 || abs(newll)==Inf || sum(postz)!= n){
                cat("Need new starting values due to singularity...","\n")
			restarts<-restarts+1
			if(restarts>15) stop("Too many tries!")
                tmp <- multmix.init(y=y, k=k)
                    lambda <- tmp$lambda
                theta <- tmp$theta
                    k <- tmp$k
                    loglik <- numeric(maxit)
                    oldll = -Inf
                    iter <- 0
                    diff <- epsilon * 2
                ll <- NULL
            } else {        if (verb) {
                        cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                         loglik[iter], "\n") }
        
        }
        }
    }
#        for (j in 1:k) {
#            postz[, j] = lambda[j] * exp(apply(y, 1, function(x, 
#                th) ldmult(x, th), th = theta[j, ]))
#        }
#
#            dens = apply(postz, 1, sum)
#            newll = loglik[iter] = sum(log(dens))
#            postz <- postz/dens
#        postz <- cbind(postz[,1:(k-1)],1-apply(as.matrix(postz[,1:(k-1)]),1,sum))
#        ll <- c(ll, newll)


        if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    theta[,p] <- 1-apply(as.matrix(theta[,1:(p-1)]),1,sum)
colnames(theta) <- c(paste("theta", ".", 1:p, sep = ""))
rownames(theta) <- c(paste("comp", ".", 1:k, sep = ""))
colnames(postz) <- c(paste("comp", ".", 1:k, sep = ""))
    cat("number of iterations=", iter, "\n")
    list(y=y, lambda = lambda, 
        theta = theta, loglik = sum(log(dens)), posterior = postz, all.loglik=ll, restarts=restarts, ft="multmixEM")
}