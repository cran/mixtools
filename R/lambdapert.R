lambda.pert <- function(lambda,pert){

k=length(lambda)


if(k>2){
temp=c(lambda[1],sapply(2:(k-1),function(i) lambda[i]/(1-sum(lambda[1:(i-1)]))))
} else{
temp=lambda[1]
}


temp2=logit(temp)+pert

temp3=inv.logit(temp2)


new.lambda=NULL

if(k>2){
new.lambda[1]=temp3[1]
for(i in 1:(k-1)){
	new.lambda[i]=temp3[i]*(1-sum(new.lambda[1:(i-1)]) )
}
new.lambda[k]=1-sum(new.lambda)
} else{
new.lambda=c(temp3,1-temp3)
}


new.lambda

}