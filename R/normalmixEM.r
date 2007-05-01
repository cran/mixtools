normalmixEM = 
function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbmean = TRUE, arbvar = TRUE, 
    epsilon = 1e-08, maxit = 10000, verb = FALSE) 
{
    if(arbmean == FALSE && arbvar == FALSE){
    stop(paste("Must change constraints on mu and/or sigma!","\n"))
    }
    x <- as.vector(x)
    n <- length(x)
    tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
        k = k, arbmean = arbmean, arbvar = arbvar)
    lambda <- tmp$lambda
    mu <- tmp$mu
    k <- tmp$k
    sigma <- tmp$s
    var <- sigma^2
    diff <- 1
    iter <- 0
    comp <- lapply(1:k, function(i) lambda[i] * dnorm(x, mu[i * arbmean + (1 - arbmean)], 
        sigma[i * arbvar + (1 - arbvar)]))
    comp <- sapply(comp, cbind)
    compsum <- apply(comp, 1, sum)
    obsloglik <- sum(log(compsum))
    ll <- obsloglik
    z = matrix(nrow = n, ncol = k)
	restarts <- 0
    while (diff > epsilon && iter < maxit) {
        res <- outer(x, mu, "-")^2
        for (i in 1:n) {
            for (j in 1:k) {
                z.denom = c()
                for (h in 1:k) {
                  z.denom = c(z.denom, (lambda[h]/lambda[j]) * 
                    (sigma[j * arbvar + (1 - arbvar)]/sigma[h * 
                      arbvar + (1 - arbvar)]) * exp(-0.5 * ((1/sigma[h * 
                    arbvar + (1 - arbvar)]^2) * res[i, h * arbmean + (1 - arbmean)] - (1/sigma[j * 
                    arbvar + (1 - arbvar)]^2) * res[i, j * arbmean + (1 - arbmean)])))
                }
                z[i, j] = 1/sum(z.denom)
            }
        }
	
#	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)
	z=z/apply(z,1,sum)

        lambda.new <- apply(z, 2, mean)
        sing <- sum(is.nan(z))
        if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
            sing <- 1
        }
        else {
            if (arbmean) {
                mu.new <- apply(z * x, 2, sum)/apply(z, 2, sum)
            }
            else {
                mu.new <- sum(apply(z * x, 2, sum))/n
            }
            if (arbvar) {
                sigma2.new <- sapply(1:k, function(i) sum(z[, 
                  i] * (x - mu.new[i * arbmean + (1 - arbmean)])^2)/apply(z, 2, sum)[i])
            }
            else {
                sigma2.new = sum(sapply(1:k, function(i) sum(z[, 
                  i] * (x - mu.new[i * arbmean + (1 - arbmean)])^2)/n))
            }
            lambda <- lambda.new
            mu <- mu.new
            sigma <- sqrt(sigma2.new)
            sing <- sum(sigma < 1e-08)
            comp <- lapply(1:k, function(i) lambda[i] * dnorm(x, 
                mu[i * arbmean + (1 - arbmean)], sigma[i * arbvar + (1 - arbvar)]))
            comp <- sapply(comp, cbind)
            compsum <- apply(comp, 1, sum)
            newobsloglik <- sum(log(compsum))
        }
        if (sing > 0 || abs(newobsloglik) == Inf || is.na(newobsloglik) || 
            sum(z) != n) {
            cat("Need new starting values due to singularity...", 
                "\n")
		restarts <- restarts + 1
		if(restarts>15) stop("Too many tries!")
            tmp <- normalmix.init(x = x, k = k, arbmean = arbmean, arbvar = arbvar)
            lambda <- tmp$lambda
            mu <- tmp$mu
            k <- tmp$k
            sigma <- tmp$s
            var <- sigma^2
            diff <- 1
            iter <- 0
            comp <- lapply(1:k, function(i) lambda[i] * dnorm(x, 
                mu[i * arbmean + (1 - arbmean)], sigma[i * arbvar + (1 - arbvar)]))
            comp <- sapply(comp, cbind)
            compsum <- apply(comp, 1, sum)
            obsloglik <- sum(log(compsum))
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
    scale.order = order(sigma)
    sigma.min = min(sigma)
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if(arbmean == FALSE){
    z = z[,scale.order]
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, lambda = lambda[scale.order], mu = mu, sigma = sigma.min, scale = sigma[scale.order]/sigma.min, loglik = obsloglik, 
        posterior = z, all.loglik=ll, restarts=restarts, ft="normalmixEM")
    class(a) = "mixEM"
    a
    } else {
    colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, lambda = lambda, mu = mu, sigma = sigma, loglik = obsloglik, 
        posterior = z, all.loglik=ll, restarts=restarts, ft="normalmixEM")
    class(a) = "mixEM"
    a
    }
}
