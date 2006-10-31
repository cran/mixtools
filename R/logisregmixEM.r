logisregmixEM = function (y, x, N = NULL, lambda = NULL, beta = NULL, k = 2, addintercept = TRUE, 
    epsilon = 1e-08, maxit = 10000, verb=FALSE) 
{
    logit <- function(x) 1/(1 + exp(-x))
    if (addintercept) {
        x = cbind(1, x)
    } else x = as.matrix(x)
    n <- length(y)
    p <- ncol(x)
    if(is.null(N)){
    N = rep(1,n)
    }
    tmp <- logisregmix.init(y=y, x=x, N=N, lambda=lambda, beta=beta, k=k)
    lambda <- tmp$lambda
    beta <- tmp$beta
    k <- tmp$k

    xbeta <- x %*% beta
    z <- matrix(0, n, k)
    diff <- 1
    iter <- 0
    comp <- t(t(dbinom(y, size = N, prob = logit(xbeta))) * lambda)
    compsum <- apply(comp, 1, sum)
    obsloglik <- sum(log(compsum))
    ll<-obsloglik
    while (diff > epsilon && iter < maxit) {
        j.star = apply(xbeta, 1, which.max)
        for (i in 1:n) {
            for (j in 1:k) {
                z[i, j] = lambda[j]/lambda[j.star[i]] * (logit(xbeta[i, 
                  j])/logit(xbeta[i, j.star[i]]))^y[i] * ((1 - 
                  logit(xbeta[i, j]))/(1 - logit(xbeta[i, j.star[i]])))^(N[i] - 
                  y[i])
            }
        }
        z = z/apply(z, 1, sum)
    if(sum(is.nan(z))>0){
        cat("Need new starting values due to underflow...","\n")
        tmp <- logisregmix.init(y=y, x=x, N=N, k=k)
        lambda <- tmp$lambda
        beta <- tmp$beta
        k <- tmp$k
          diff <- 1
            iter <- 0
        xbeta <- x %*% beta
        comp <- t(t(dbinom(y, size = N, prob = logit(xbeta))) * lambda)
        compsum <- apply(comp, 1, sum)
        obsloglik <- sum(log(compsum))
        ll<-obsloglik
            
    } else{
            lambda <- apply(z, 2, mean)
            lm.out = lapply(1:k, function(j) try(glm.fit(x, cbind(y, N - y), weights = z[, 
                j], family = binomial()), silent=TRUE))
        beta = sapply(lm.out,coef)        
      xbeta <- x %*% beta
        comp <- t(t(dbinom(y, size = N, prob = logit(xbeta))) * 
            lambda)
        compsum <- apply(comp, 1, sum)
        newobsloglik <- sum(log(compsum))

    if(abs(newobsloglik)==Inf || is.nan(newobsloglik) || sum(z)!=n){
        cat("Need new starting values due to singularity...","\n")
        tmp <- logisregmix.init(y=y, x=x, N=N, k=k)
        lambda <- tmp$lambda
        beta <- tmp$beta
        k <- tmp$k
          diff <- 1
            iter <- 0
        xbeta <- x %*% beta
        comp <- t(t(dbinom(y, size = N, prob = logit(xbeta))) * lambda)
        compsum <- apply(comp, 1, sum)
        obsloglik <- sum(log(compsum))
        ll <- obsloglik
    } else{
      diff <- newobsloglik - obsloglik
        obsloglik <- newobsloglik
    ll <- c(ll, obsloglik)
        iter <- iter + 1
        if (verb) {
            cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                obsloglik, "\n") }

    }
}
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
rownames(beta) <- c(paste("beta", ".", 0:(p-1), sep = ""))
colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    cat("number of iterations=", iter, "\n")
    a=list(x=x, y=y, lambda = lambda, beta = beta, loglik = obsloglik, posterior = z, all.loglik=ll, 
    ft="logisregmixEM")
    class(a) = "mixEM"
    a
}
