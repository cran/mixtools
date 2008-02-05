plot.npEM <- function(x, blocks = NULL, hist=TRUE, 
                      scale = TRUE, title=NULL, breaks="Sturges", 
                      ylim=NULL, dens.col, ...) {
  ask <- par(ask=TRUE)
  r <- NCOL(x$data)
  m <- NCOL(x$posteriors)
  blockid <- x$blockid
  if (is.null(blocks)) {
    if(!is.null(blockid)) {
      blocks <- 1:max(blockid)
    } else {
      blocks <- blockid <- 1:r
    }
  }
  ylim.orig <- ylim
  out <- list(x=list(), y=list())
  for(i in 1:length(blocks)) {
    coords <- blockid == blocks[i]
    ylim <- ylim.orig
    if (is.null(title)) {
      tt <- paste(which(coords), collapse=",")
      tt <- paste("Coordinate", ifelse(sum(coords)>1, "s ", " "), tt, sep="")
    } else { 
      tt <- rep(title,length(blocks))[i]                  
    }
    dx <- dy <- NULL
    for (j in 1:m) {
      d <- density(x, component=j, block=i, scale=scale)
      dx <- cbind(dx, d$x)
      dy <- cbind(dy, d$y)
    }
    if (is.null(ylim)) {
      ylim=range(dy)
      if (hist) {
        xx <- as.vector(as.matrix(x$data)[,coords])
        ylim[2] <- max(ylim[2], hist(xx, breaks=breaks, plot=FALSE)$density)
      }
    }
    pf <- plot # Use plot or hist as plotting fn the 1st time only, then lines
    if (hist) {
      pf <- lines
      hist(xx, breaks=breaks, prob=TRUE, main=tt, ylim=ylim, ...)
    }
    if (missing(dens.col)) 
      dens.col <- 2:(m+1)
    dens.col <- rep(dens.col, length.out=m)
    for (j in 1:m) {
      pf(dx[,j],dy[,j], type="l", lwd=2, col=dens.col[j], main=tt, ylim=ylim, ...)
      pf <- lines
    }
    legend("topright", legend=round(x$lambdahat,3), fill=dens.col)
    out$x[[i]]<-dx
    out$y[[i]]<-dy
  }
  par(ask=ask)
  invisible(out)
}

    

