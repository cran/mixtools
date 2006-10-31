# Simulate from a normal mixture.  If m>1, then each row in the result
# is sampled from the same component (as in repeated measures).

normmix.sim = function(n,lambda,mu,sigma,m=1) {
  k=length(lambda)
  comp=sample(1:k,n,rep=TRUE,prob=lambda)
  matrix(rnorm(n*m,mean=mu[comp],sd=sigma[comp]),n,m)
}

