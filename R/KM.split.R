#KMplot function
#' @title Kaplan-Meier estimator of binary splitting
#' @description Given a cut-off-point and selected covariate, return the survival curve for binary classification and the P-value of two sample log-rank test.
#' @param t.vec :Vector of survival times (time to either death or censoring)
#' @param d.vec :Vector of censoring indicators (1=death, 0=censoring)
#' @param X.mat :n by p matrix of covariates, where n is the sample size and p is the number of covariates
#' @param x.name :the name of covariate
#' @param cutoff :cut-off-point
#' @return P-value of two sample logrank test and a plot of two KM estimates
#' @examples
#' data(Lung,package="compound.Cox")
#' train_Lung=Lung[which(Lung[,"train"]==TRUE),] #select training data
#' t.vec=train_Lung[,1]
#' d.vec=train_Lung[,2]
#' x.mat=train_Lung[,-c(1,2,3)]
#' KM.split(t.vec,d.vec,x.mat,x.name="ANXA5",cutoff=1)
#' @import survival
#' @import plot
#' @export

KM.split=function(t.vec,d.vec,X.mat,x.name,cutoff){
  group=(X.mat[,x.name]<=cutoff)*1
  left.time=t.vec[which(group==1)]
  left.cens=d.vec[which(group==1)]
  left.fit=survfit(Surv(left.time,left.cens)~1)
  plot(left.fit,conf.int = F,mark.time = T,lwd=3,col="red",xlim=c(0,max(t.vec)),main="Kaplan Meier Estimator")
  par(new=TRUE)
  right.time=t.vec[which(group==0)]
  right.cens=d.vec[which(group==0)]
  right.fit=survfit(Surv(right.time,right.cens)~1)
  plot(right.fit,conf.int = F,mark.time=T,lwd=3,col="blue",xlim=c(0,max(t.vec)))
  legend("bottomright",legend=c(paste(x.name,"<=",cutoff,sep = ""),paste(x.name,">",cutoff,sep = ""))
         ,col=c("red","blue"),lty=1:1,cex=1)
  p=pchisq(survdiff(Surv(t.vec,d.vec)~group)$chisq,1,lower.tail = F)
  names(p)=paste("P-value of"," H0:h(t|",x.name,"<=",cutoff,")=h(t|",x.name,">",cutoff,")",sep="")
  text(quantile(t.vec)[4],0.9,paste("P-value=",as.character(round(p,5))),cex=1.3)
  print(p)
}
