uni.tree=function(t.vec,d.vec,X.mat,P.value=0.01,d0=0.01,
                  S.plot=FALSE,score=TRUE){
  c.vec=min(X.mat):(max(X.mat)-1) ## cut-off values
  m=length(c.vec) ## no. of cut-off values
  recursive=function(t.vec,d.vec,X.mat,condition=NULL,Risk.score=0,rank.weight=1){
    n=length(t.vec)
    if(n>1){
      if(score==TRUE){
        p=ncol(X.mat)
        w.mat=1*(X.mat[,rep(1:p,rep(m,p))]>matrix(rep(c.vec,p),1,p*m)[rep(1,n),])
        uni.score.res=uni.score(t.vec,d.vec,w.mat,d0=d0)
        P.values=uni.score.res$P
        cut.value=rep(c.vec,p)[which.min(P.values)]
        cut.name=colnames(X.mat)[floor((which.min(P.values)-1)/m)+1]
        cut.P = min(P.values)
      }else{
        uni.logrank.res=uni.logrank(t.vec,d.vec,X.mat)
        cut.value=uni.logrank.res["cut_off_point",1]
        cut.name=colnames(uni.logrank.res)[1]
        cut.P=uni.logrank.res["Pvalue",1]
      }
      #stopping criterion
      if(cut.P <= P.value ){
        tmp=X.mat[,cut.name]<=cut.value
        if(score == TRUE){
          Zvalue = uni.score.res$Z[which.min(P.values)]
          if(Zvalue <= 0){
            t.left = t.vec[tmp];t.right=t.vec[!tmp]
            d.left = d.vec[tmp];d.right=d.vec[!tmp]
            X.left = X.mat[tmp,];X.right=X.mat[!tmp,]
          }else{
            t.left=t.vec[!tmp];t.right=t.vec[tmp]
            d.left=d.vec[!tmp];d.right=d.vec[tmp]
            X.left=X.mat[!tmp,];X.right=X.mat[tmp,]
          }
        }else{
          x=X.mat[,cut.name]<=cut.value
          lrsummary = survdiff(Surv(t.vec,d.vec)~x)
          Zvalue = (lrsummary$obs - lrsummary$exp)[1]
          if(Zvalue <= 0){
            t.left=t.vec[tmp];t.right=t.vec[!tmp]
            d.left=d.vec[tmp];d.right=d.vec[!tmp]
            X.left=X.mat[tmp,];X.right=X.mat[!tmp,]
          }else{
            t.left=t.vec[!tmp];t.right=t.vec[tmp]
            d.left=d.vec[!tmp];d.right=d.vec[tmp]
            X.left=X.mat[!tmp,];X.right=X.mat[tmp,]
          }
        }
        #prepare for subtree
        if(is.null(condition)){
          if(Zvalue <= 0){
            condition_left=paste0(cut.name,"<=",cut.value,sep="")
            condition_right=paste0(cut.name,">",cut.value,sep="")
          }else{
            condition_left=paste0(cut.name,">",cut.value,sep="")
            condition_right=paste0(cut.name,"<=",cut.value,sep="")
          }
        }else{
          if(Zvalue <= 0){
            condition_left=paste0(condition," & ",paste0(cut.name,"<=",cut.value))
            condition_right=paste0(condition," & ",paste0(cut.name,">",cut.value))
          }else{
            condition_left=paste0(condition," & ",paste0(cut.name,">",cut.value))
            condition_right=paste0(condition," & ",paste0(cut.name,"<=",cut.value))
          }
        }
        if(S.plot == TRUE){
          KM.split(t.vec,d.vec,X.mat,cut.name,cut.value)
        }
        return(
          list(NODE=data.frame(Information=c("Inner.node",
                                             Risk.score,
                                             cut.P,
                                             cut.name,
                                             cut.value,
                                             round(Zvalue,2)
          ),
          row.names = c("node_status:","Risk score:","P-value:","name:","cut_value:","Z-value")
          )
          ,Left=recursive(t.left,d.left,X.left,condition=condition_left,
                          Risk.score=Risk.score+rank.weight,rank.weight=0.1*rank.weight)
          ,Right=recursive(t.right,d.right,X.right,condition=condition_right,
                           Risk.score=Risk.score-rank.weight,rank.weight=0.1*rank.weight)
          )
        )
      }else{
        return(list(NODE=data.frame(Information=c("terminal node",
                                                  Risk.score,
                                                  cut.P,
                                                  n,
                                                  if(is.null(condition)){"NULL"}else{condition}
        ),
        row.names=c("node_status:","Risk score:","P-value:","samplesize:","condition:")
        )
        )
        )

      }
    }else{
      return(list(NODE=data.frame(Information=c("terminal node",
                                                Risk.score,
                                                n,
                                                if(is.null(condition)){"NULL"}else{condition}
      ),
      row.names=c("node_status:","Risk score:","samplesize:","condition:")
      )
      )
      )
    }
  }
  return(recursive(t.vec,d.vec,X.mat))
}
