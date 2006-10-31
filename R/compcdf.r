

compCDF<-function(x,weights){
   x<-t(x)
   k<-ncol(weights)

   w.cdf<-function(pt,x,w){

    w<-outer(rep(1,nrow(x)),w)

    weighted.mean(x<=pt,w,na.rm=TRUE)

    }
   tdex<-seq(min(na.omit(c(x)))-0.5,max(na.omit(c(x)))+0.2,len=500)

   weights<-matrix(weights,ncol=k)

   sqk<-floor(sqrt(k)+.99)

   n<-length(na.omit(c(x)))

   mean.data<-1:k

   var.data<-1:k

   mean.data[1]<-weighted.mean(x,outer(rep(1,nrow(x)),weights[,1]),na.rm=TRUE)

   var.data[1]<-weighted.mean((x-mean.data[1])^2,

               outer(rep(1,nrow(x)),weights[,1]),na.rm=TRUE)

   temp.cdf<-apply(matrix(tdex,ncol=1),1,w.cdf,x=x,w=weights[,1])

   temp.sd<-sqrt(temp.cdf*(1-temp.cdf)/n)

   plot(tdex,temp.cdf,xlab=" ",ylab="CDF",type="l",col=1,lty=1,main="Components' CDF")

   if(k>1){

      for(i in 2:k){

         temp.cdf<-apply(matrix(tdex,ncol=1),1,x=x,w.cdf,w=weights[,i])

         temp.sd<-sqrt(temp.cdf*(1-temp.cdf)/n)

         lines(tdex,temp.cdf,col=1,lty=i)

          mean.data[i]<-weighted.mean(x,outer(rep(1,nrow(x)),

                                              weights[,i]),na.rm=TRUE)

         var.data[i]<-weighted.mean((x-mean.data[i])^2,

         outer(rep(1,nrow(x)),weights[,i]),na.rm=TRUE)

      }

   }

   result=rbind(mean.data,sqrt(var.data))
   rownames(result) = c("Mean", "StDev")

   list(result=result)

}



