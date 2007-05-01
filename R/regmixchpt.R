regmixEM.chgpt <- function (y, x, lambda = NULL, gamma = NULL, beta = NULL, sigma = NULL, k = 2, T = 3, t=NULL,
    epsilon = 1e-08, maxit = 10000, verb=FALSE) 
{

    x = cbind(1, x)
    n <- length(y)
    p <- ncol(x)

	if(is.null(t)) {
		stop(paste("You must specify a changepoint value!", "\n")) }


tmp<-regmix.chgpt.init(y=y, x=x, lambda=lambda, beta=beta, gamma=gamma, sigma=sigma, t=t, k=k)
lambda<-tmp$lambda
beta<-as.vector(tmp$beta)
gamma<-as.vector(tmp$gamma)
sigma<-tmp$sigma
    diff <- 1
    iter <- 0
	xbeta <- x%*%beta
	comp2 <- lambda[2]*dnorm(y,mean=xbeta,sd=sigma[2])

	obsloglik <- NULL
	seq <- T:(n-(T+1))
	for(i in seq){
		x.t <- (x[,2]-x[i,2])*(x[,2]>x[i,2])
		xgamma <- cbind(x,x.t)%*%t(t(gamma))
		comp1 <- lambda[1]*dnorm(y,mean=xgamma,sd=sigma[1])
		comp <- cbind(comp1,comp2)
    		obsloglik[i-(T-1)] <- sum(log(apply(comp, 1, sum)))
 	}

	if(sum(obsloglik==-Inf)==length(obsloglik)){
		cutpt <- seq[floor(mean(c(T,(n-T-1))))]}
	else cutpt <- seq[which.max(obsloglik)]
		x.t <- (x[,2]-x[cutpt,2])*(x[,2]>x[cutpt,2])
		xgamma <- cbind(x,x.t)%*%t(t(gamma))
		comp1 <- lambda[1]*dnorm(y,mean=xgamma,sd=sigma[1])
	comp <- cbind(comp1,comp2)
	compsum <- apply(comp,1,sum)
	obsloglik <- max(obsloglik,na.rm=T)
	ll<-obsloglik
	restarts<-0

    while (diff > epsilon && iter < maxit) {
	  z <- comp/compsum
	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)

	lambda.new <- apply(z,2,mean)

		gamma.new <- try(solve(t(cbind(x,x.t)) %*% (as.vector(z[,1]) * 
            cbind(x,x.t))) %*% (t(cbind(x,x.t)) %*% (as.vector(z[,1]) * y)),silent=TRUE)

	
		beta.new <- try(solve(t(x) %*% (as.vector(z[,2]) * 
            x)) %*% (t(x) %*% (as.vector(z[,2]) * y)),silent=TRUE)



	if(sum(sigma < 1e-08)>0 || sum(lambda < 1e-08)>0 || class(gamma.new)=="try-error" || class(beta.new)=="try-error"
		|| sum(is.na(gamma.new)>0) || sum(is.na(beta.new)>0)){
            	cat("Need new starting values due to singularity...", "\n")
		restarts <- restarts + 1
		if(restarts>15) stop("Too many tries!")	
		tmp<-regmix.chgpt.init(y=y, x=x, t=t, k=k)
		lambda<-tmp$lambda
		beta<-tmp$beta
		gamma<-tmp$gamma
		sigma<-tmp$sigma
		diff <- 1
    		iter <- 0
		xbeta <- x%*%beta
		comp2 <- lambda[2]*dnorm(y,mean=xbeta,sd=sigma[2])

		obsloglik <- NULL
		seq <- T:(n-(T+1))
		for(i in seq){
			x.t <- (x[,2]-x[i,2])*(x[,2]>x[i,2])
			xgamma <- cbind(x,x.t)%*%t(t(gamma))
			comp1 <- lambda[1]*dnorm(y,mean=xgamma,sd=sigma[1])
			comp <- cbind(comp1,comp2)
    			obsloglik[i-(T-1)] <- sum(log(apply(comp, 1, sum)))
 		}

		if(sum(obsloglik==-Inf)==length(obsloglik)){
				cutpt <- seq[floor(mean(c(T,(n-T-1))))]}
			else cutpt <- seq[which.max(obsloglik)]
		cutpt <- seq[which.max(obsloglik)]
		x.t <- (x[,2]-x[cutpt,2])*(x[,2]>x[cutpt,2])
		xgamma <- cbind(x,x.t)%*%t(t(gamma))
		comp1 <- lambda[1]*dnorm(y,mean=xgamma,sd=sigma[1])
		comp <- cbind(comp1,comp2)
		compsum <- apply(comp,1,sum)
		obsloglik <- max(obsloglik,na.rm=T)
		ll<-obsloglik


	} else{
	sigma.new <- NULL
		sigma.new[1] <- sqrt(sum(as.vector(z[,1])*(y-cbind(x,x.t)%*%gamma.new)^2)/sum(z[,1]))
		sigma.new[2] <- sqrt(sum(as.vector(z[,2])*(y-x%*%beta.new)^2)/sum(z[,2]))

	gamma <- gamma.new
	beta <- beta.new
	sigma <- sigma.new
	lambda <- lambda.new

	xbeta <- x%*%beta
	comp2 <- lambda[2]*dnorm(y,mean=xbeta,sd=sigma[2])

	newobsloglik <- NULL
	for(i in seq){
		x.t <- (x[,2]-x[i,2])*(x[,2]>x[i,2])
		xgamma <- cbind(x,x.t)%*%t(t(gamma))
		comp1 <- lambda[1]*dnorm(y,mean=xgamma,sd=sigma[1])
		comp <- cbind(comp1,comp2)
    		newobsloglik[i-(T-1)] <- sum(log(apply(comp, 1, sum)))
 	}
		if(sum(newobsloglik==-Inf)==length(newobsloglik)){
				cutpt <- seq[floor(mean(c(T,(n-T-1))))]}
			else cutpt <- seq[which.max(newobsloglik)]
		x.t <- (x[,2]-x[cutpt,2])*(x[,2]>x[cutpt,2])
		xgamma <- cbind(x,x.t)%*%t(t(gamma))
		comp1 <- lambda[1]*dnorm(y,mean=xgamma,sd=sigma[1])
	comp <- cbind(comp1,comp2)
	compsum <- apply(comp,1,sum)
	newobsloglik <- max(newobsloglik,na.rm=T)
	diff <- newobsloglik-obsloglik
        obsloglik <- newobsloglik
	ll <- c(ll, obsloglik)
	iter<-iter+1
           if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  obsloglik, "\n")
            }
	}

    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")

    a=list(lambda = lambda, gamma = gamma, beta = beta, sigma = sigma, cutpoint=cutpt, loglik = obsloglik,
        posterior = z, all.loglik=ll, restarts=restarts, ft="regmixEM.chgpt")
    class(a) = "mixEM"
    a
}












