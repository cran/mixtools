regmixEM = function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2, 
    addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, maxit = 10000, 
    verb = FALSE) 
{
    if(arbmean == FALSE && arbvar == FALSE){
    stop(paste("Must change constraints on beta and/or sigma!","\n"))
    }
    s = sigma
    if (addintercept) {
        x = cbind(1, x)
    }
    n <- length(y)
    p <- ncol(x)
    tmp <- regmix.init(y = y, x = x, lambda = lambda, beta = beta, 
        s = s, k = k, addintercept = addintercept, arbmean = arbmean, arbvar = arbvar)
    lambda <- tmp$lambda
    beta <- tmp$beta
    s <- tmp$s
    k <- tmp$k
    diff <- 1
    iter <- 0
    xbeta <- x %*% beta
    res <- (y - xbeta)^2
    if(arbmean == FALSE){
    res <- sapply(1:k,function(i) res)
    }
#    comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[,i * arbmean + (1 - arbmean)], 
#        s[i * arbvar + (1 - arbvar)]))
#    comp <- sapply(comp, cbind)
#    compsum <- apply(comp, 1, sum)
#    obsloglik <- sum(log(compsum))
    comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
        s^2)))))
    obsloglik <- sum(log(apply(comp, 1, sum)))
    ll <- obsloglik
    z = matrix(nrow = n, ncol = k)
    while (diff > epsilon && iter < maxit) {
        for (i in 1:n) {
            for (j in 1:k) {
                z.denom = c()
                for (h in 1:k) {
                  z.denom = c(z.denom, (lambda[h]/lambda[j]) * 
                    (s[j * arbvar + (1 - arbvar)]/s[h * arbvar + 
                      (1 - arbvar)]) * exp(-0.5 * ((1/s[h * arbvar + 
                    (1 - arbvar)]^2) * res[i, h] - (1/s[j * arbvar + 
                    (1 - arbvar)]^2) * res[i, j])))
                }
                z[i, j] = 1/sum(z.denom)
            }
        }
        lambda.new <- apply(z, 2, mean)
        if (sum(lambda.new < 1e-08)>0) {
            sing <- 1
        }
        else {
            if (addintercept) {
                lm.out <- lapply(1:k, function(i) lm(y ~ x[, 
                  -1], weights = z[, i]))
            }
            else lm.out <- lapply(1:k, function(i) lm(y ~ x - 
                1, weights = z[, i]))
            beta.new <- sapply(lm.out, coef)
        if (arbmean == FALSE) {
        beta.new <- apply(t(apply(z,2,sum)*t(beta.new)),1,sum)/n
        }
        xbeta.new <- x %*% beta.new
        res <- (y - xbeta.new)^2
            if(arbmean == FALSE){
            res <- sapply(1:k,function(i) res)
            }
            if (arbvar) {
                s.new <- sqrt(sapply(1:k, function(i) sum(z[, 
                  i] * (res[, i]))/sum(z[, i])))
            }
            else s.new <- sqrt(sum(z * res)/n)
            lambda <- lambda.new
            beta <- beta.new
          xbeta <- x%*%beta
            s <- s.new
            sing <- sum(s < 1e-08)
    comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[,i * arbmean + (1 - arbmean)], 
        s[i * arbvar + (1 - arbvar)]))
    comp <- sapply(comp, cbind)
    compsum <- apply(comp, 1, sum)
    newobsloglik <- sum(log(compsum))
#            comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
#                s^2)))))
#            newobsloglik <- sum(log(apply(comp, 1, sum)))
        }
        if (newobsloglik < obsloglik || sing > 0 || abs(newobsloglik) == 
            Inf || is.nan(newobsloglik) || sum(z) != n) {
            cat("Need new starting values due to singularity...", 
                "\n")
            tmp <- regmix.init(y = y, x = x, k = k, addintercept = addintercept, 
                arbmean = arbmean, arbvar = arbvar)
            lambda <- tmp$lambda
            beta <- tmp$beta
            s <- tmp$s
            k <- tmp$k
            diff <- 1
            iter <- 0
            xbeta <- x %*% beta
            res <- (y - xbeta)^2
    if(arbmean == FALSE){
    res <- sapply(1:k,function(i) res)
    }
#    comp <- lapply(1:k, function(i) lambda[i] * dnorm(y, xbeta[,i * arbmean + (1 - arbmean)], 
#        s[i * arbvar + (1 - arbvar)]))
#    comp <- sapply(comp, cbind)
#    compsum <- apply(comp, 1, sum)
#    obsloglik <- sum(log(compsum))
            comp <- t((lambda/sqrt(2 * pi * s^2)) * t(exp(-t(t(res)/(2 * 
                s^2)))))
            obsloglik <- sum(log(apply(comp, 1, sum)))
    ll <- obsloglik

        }
        else {
            diff <- newobsloglik - obsloglik
            obsloglik <- newobsloglik
          ll <- c(ll, obsloglik)
            iter <- iter + 1
            if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  obsloglik, "\n")
            }
        }
    }
    scale.order = order(s)
    sigma.min = min(s)
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if(arbmean == FALSE){
    z=z[,scale.order]
names(beta) <- c(paste("beta", ".", 0:(p-1), sep = ""))
colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, y=y, lambda = lambda[scale.order], beta = beta, sigma = sigma.min, scale = s[scale.order]/sigma.min, loglik = obsloglik, 
        posterior = z[,scale.order], all.loglik=ll, ft="regmixEM")
    class(a) = "mixEM"
    a
    } else {
rownames(beta) <- c(paste("beta", ".", 0:(p-1), sep = ""))
colnames(beta) <- c(paste("comp", ".", 1:k, sep = ""))
colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, y=y, lambda = lambda, beta = beta, sigma = s, loglik = obsloglik, 
        posterior = z, all.loglik=ll, ft="regmixEM")
    class(a) = "mixEM"
    a
    }
}
