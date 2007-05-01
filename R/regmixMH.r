regmixMH=function (y, x, lambda = NULL, beta = NULL, s = NULL, k = 2, 
    addintercept = TRUE, mu = NULL, sig = NULL, sampsize = 1000, 
    omega = 0.01, thin = 1) 
{
    if (addintercept) {
        x = cbind(1, x)
    }
    n <- length(y)
    p <- ncol(x)
    if (is.null(s)) {
        s = sqrt(1/rexp(k))
    }
    else k = length(s)
    if (is.null(beta)) {
        beta = matrix(rnorm(p * k), p, k)
    }
    else k = ncol(beta)
    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    }
    else k = length(lambda)
    if (is.null(mu)) {
        mu = 0 * beta
    }
    if (is.null(sig)) {
        sig = 5 * sqrt(var(y)) + 0 * beta
    }
    L.theta <- matrix(nrow = n, ncol = k)
    pi.beta <- matrix(nrow = p, ncol = k)
    pi.sigma <- c()
    pi.lambda <- c()
    new.L.theta <- matrix(nrow = length(y), ncol = k)
    new.pi.beta <- matrix(nrow = p, ncol = k)
    new.pi.sigma <- c()
    new.pi.lambda <- c()
    accepts = 0
    theta <- matrix(c(beta, s, lambda), nrow = 1)
    thetalist <- NULL
    for (i in 2:sampsize) {
        pi.beta <- dnorm(beta, mu, sig)
        pi.sigma <- dexp(s)
        pi.lambda <- dbeta(lambda, 1, 1)
        L.theta <- dnorm(y - x %*% beta, 0, matrix(s, n, k, byrow = TRUE)) %*% 
            matrix(lambda, k, 1)
        Lik.theta <- prod(apply(L.theta, 1, sum))
        prior <- prod(pi.beta) * prod(pi.sigma) * prod(pi.lambda)
        f.theta <- Lik.theta * prior
        new.beta <- beta + omega * matrix(rcauchy(k * p), p, 
            k)
        new.sigma <- log(s) + omega * rcauchy(k)
        new.sigma <- exp(new.sigma)
#        temp <- exp(log(lambda[1:(k - 1)]) - log(lambda[k]) + 
#            omega * rcauchy(k - 1))
#        new.lambda <- c(temp, 1)/(1 + sum(temp))
	new.lambda <- lambda.pert(lambda,omega*rcauchy(k-1))
        new.pi.beta <- dnorm(new.beta, mu, sig)
        new.pi.sigma <- dexp(new.sigma)
        new.pi.lambda <- dbeta(new.lambda, 1, 1)
        new.L.theta <- dnorm(y - x %*% new.beta, 0, matrix(new.sigma, 
            n, k, byrow = TRUE)) %*% matrix(new.lambda, k, 1)
        new.Lik.theta <- prod(apply(new.L.theta, 1, sum))
        new.prior <- prod(new.pi.beta) * prod(new.pi.sigma) * 
            prod(new.pi.lambda)
        new.f.theta <- new.Lik.theta * new.prior
        new.theta <- c(new.beta, new.sigma, new.lambda)
        a <- new.f.theta/f.theta
        r <- runif(1)
        if (a > r) {
            theta <- new.theta
            beta <- new.beta
            s <- new.sigma
            lambda <- new.lambda
            accepts = accepts + 1
        }
        if (i%%thin == 0) 
            thetalist = rbind(thetalist, theta)
    }
    cat(paste("Acceptance rate: ", 100 * accepts/sampsize, "%\n", 
        sep = ""))
    colnames(thetalist) <- c(paste("beta", rep(0:(p - 1), k), 
        ".", rep(1:k, rep(p, k)), sep = ""), paste("s.", 1:k, 
        sep = ""), paste("lambda.", 1:k, sep = ""))
    invisible(thetalist)
	a=list(x=x, y=y, theta=thetalist, components=k)
	class(a)="mixMCMC"
	a
}