KM.split=function(t.vec,d.vec,X.mat,x.name,cutoff){
  group=(X.mat[,x.name]<=cutoff)*1
  left.time=t.vec[which(group==1)]
  left.cens=d.vec[which(group==1)]
  left.fit=survfit(Surv(left.time,left.cens)~1)
  plot(left.fit,conf.int=F,mark.time=T,lwd=3,col="red",xlim=c(0,max(t.vec)))
  par(new=TRUE)
  right.time=t.vec[which(group==0)]
  right.cens=d.vec[which(group==0)]
  right.fit=survfit(Surv(right.time,right.cens)~1)
  plot(right.fit,conf.int=F,mark.time=T,lwd=3,col="blue",xlim=c(0,max(t.vec)))
  legend("bottomright",legend=c(paste(x.name,"<=",cutoff,sep=""),
                                paste(x.name,">",cutoff,sep=""))
         ,col=c("red","blue"),lty=1:1,cex=1)
  p=pchisq(survdiff(Surv(t.vec,d.vec)~group)$chisq,1,lower.tail=F)
  names(p)=paste("P-value of",
                 " H0:h(t|",x.name,"<=",cutoff,")=h(t|",x.name,">",cutoff,")",sep="")
  text(quantile(t.vec)[4],0.9,paste("P-value=",as.character(round(p,5))),cex=1.3)
  print(p)
}
