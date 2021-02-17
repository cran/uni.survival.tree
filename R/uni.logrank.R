#uni.logrank#
#' @title Univariate selection under log-rank test
#' @description List out all the possible variables. The cut-off-point has been optimal.
#' @param t.vec :Vector of survival times (time to either death or censoring)
#' @param d.vec :Vector of censoring indicators (1=death, 0=censoring)
#' @param X.mat :n by p matrix of covariates, where n is the sample size and p is the number of covariates
#' @return A table contains an information of binary splitting, including P-value of the test, cut-off-point, sample sizes of two nodes, for each covariate
#' @return Pvalue: P-value of two sample logrank test
#' @return cut_off_point: to transfer x into binary case
#' @return left.sample.size: sample size of left child node
#' @return right.sample.size: sample size of right child node
#' @details the output of the covariates has been ordered by the P-value from the smallest to the largest
#' @examples
#' data(Lung,package="compound.Cox")
#' train_Lung=Lung[which(Lung[,"train"]==TRUE),] #select training data
#' t.vec=train_Lung[,1]
#' d.vec=train_Lung[,2]
#' x.mat=train_Lung[,-c(1,2,3)]
#' uni.logrank(t.vec,d.vec,x.mat)
#' @import survival
#' @export
uni.logrank=function(t.vec,d.vec,X.mat){
  n=length(t.vec) #sample size
  if(n == 1){
    warning('there is only one sample')
  }else{
    ################this step is to avoid some column that contains only one value##################
    #tmp=sapply(x.mat,function(y){length(unique(y))})>1
    tmp=c()
    for(k in 1:ncol(X.mat)){
      tmp[k]=length(unique(X.mat[,k]))>1
    }
    X.mat = X.mat[,tmp]
    ################################################################################################
    ####initialize the vector of output
    Pvalue = left.sample.size = right.sample.size = cut_off_point =
      array(rep(0,ncol(X.mat)),dimnames = list(colnames(X.mat)))
    for(i in 1:ncol(X.mat)){
      Pvalue_comparison.array = c()
      x.level = c()
      x.unique = unique(X.mat[,i])
      for(j in 1:(length(x.unique)-1)){
        x.level[j] = sort(x.unique)[j]
        group=1*(X.mat[,i] <= x.level[j])
        logrank_res = survdiff(Surv(t.vec,d.vec)~group)
        Pvalue_comparison.array[j] = pchisq(logrank_res$chisq,1,lower.tail = F)
      }
      Pvalue[i] = min(Pvalue_comparison.array)
      cut_off_point[i] = x.level[which.min(Pvalue_comparison.array)]
      left.sample.size[i] = sum(X.mat[,i] <= cut_off_point[i])
      right.sample.size[i] = sum(X.mat[,i] > cut_off_point[i])
    }
    res_table=rbind(Pvalue,cut_off_point,left.sample.size,right.sample.size)
    return(res_table[,order(res_table["Pvalue",])])
  }
}
