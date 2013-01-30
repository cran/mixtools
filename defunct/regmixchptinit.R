regmix.chgpt.init = function (y, x, lambda = NULL, beta = NULL, gamma = NULL, sigma = NULL, t=NULL, k=2) 
{
    x <- x[, -1]
    n <- length(y)
    p <- ncol(x)

	if(is.null(lambda)){
        lambda = runif(k)
        lambda = lambda/sum(lambda)
	} else k=length(lambda)

    A = round(lambda * n)
	while(min(A)<4){
	       lambda = runif(k)
             lambda = lambda/sum(lambda)
		 A = round(lambda * n)
	}



x.1<-x
x.1.t<-(x.1-x[t])*(x.1>x[t])
y.1<-y
A.samp<-sample(1:n,A[2])
x.2<-x[A.samp]
y.2<-y[A.samp]


out.1<-lm(y.1~x.1+x.1.t)
out.2<-lm(y.2~x.2)

gamma.sim <- out.1$coef
beta.sim <- out.2$coef

s.1=sqrt(anova(out.1)$Mean[2])
s.2=sqrt(anova(out.2)$Mean[2])

s.sim=c(s.1,s.2)


    if (is.null(beta)){
	beta <- rnorm(2,mean=beta.sim,sd=s.2)
	}

    if (is.null(gamma)){
	gamma <- rnorm(3,mean=gamma.sim,sd=s.1)
        }


	if(is.null(sigma)){
		sigma<-1/rexp(2, rate = s.sim)
	}


list(lambda = lambda, beta = matrix(beta,ncol=1), gamma=matrix(gamma,ncol=1), sigma = sigma)
}


