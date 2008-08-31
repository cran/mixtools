wIQR <- function(wt=rep(1,length(x)), x, already.sorted=FALSE, already.normalized=FALSE) {
  if(!already.sorted) {
    wt <- wt[o<-order(x)]
    x <- x[o]
  }
  if(!already.normalized) {
    wt <- wt/sum(wt)
  }
  diff(x[findInterval(c(.25,.75),cumsum(wt))])
}
