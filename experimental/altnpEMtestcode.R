library(mixtools)
data(Waterdata)
b=npEM(Waterdata, mu0=3, bw=4, maxit=100)
#a=altnpEM(Waterdata, mu0=3, bw=4, maxit=50)
a=nems(Waterdata, mu0=3, bw=4, maxit=100)
f=function(j) {
 plot(b, block=j)
 for(i in 1:3) lines(a$grid, a$f[,i,j]*a$lambdahat[i], lwd=3, lty=2, col=i+1)
}

