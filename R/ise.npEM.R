ise.npEM <- function(npEMout, component=1, block=1, truepdf, lower=-Inf,
                     upper=Inf, plots = TRUE, ...){
	# returns the Integrated Squared Error
	# between f_{comp,block} as estimated by npEM,
	# and truepdf
  #true2 <- function(u) truepdf(u, ...)
	coords <- npEMout$blockid == block
	bs <- sum(coords) # block size
	xx <- as.vector(npEMout$data[,coords]) # flatten data 
	wts <- rep(npEMout$post[,component],bs) # duplicate weights
	integrand = function(u,...) {
    (wkde(xx,u,wts,npEMout$bandwidth) - truepdf(u,...))^2
  }
	numint <- integrate(integrand,lower,upper, ...)
	if (plots) {
    # plot of estimated and truepdf
    ise <- paste(round(numint$value,4))
    temp=paste(component, block, sep="")
    title = substitute(expression(paste("Integrated Squared Error for ",
                                        f[temp]," = ",ise,sep="")))
    if (!is.finite(lower)) {
      lower <- min(xx)
    }
    if (!is.finite(upper)) {
      upper <- max(xx)
    }    
    u <- seq(lower,upper, 0.01)
    fhat <- wkde(xx,u,wts,npEMout$bandwidth)
    ymax <- max(max(truepdf(u)),max(fhat))
    plot(u,fhat,type="l",ylim=c(0,ymax),
         main=eval(title),ylab="")
    lines(u,truepdf(u, ...),lwd=2,col=2)
	}
	numint
}

