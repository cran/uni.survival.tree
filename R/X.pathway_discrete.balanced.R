#X.pathway_discrete.balanced
#' @title Generate a matrix of gene expressions (discrete version of X.pathway() against to Emura (2012)) in the presence of gene pathways
#' @description Generate a matrix of gene expressions in the presence of gene pathways, we first produce the matrix by X.pathway(Emura et al. 2012), then we change each value to 1 ~ 4 depend on the quantile.
#' @param n :the number of individuals (sample size)
#' @param p :the number of genes
#' @param q1 :the number of genes in the first pathway
#' @param q2 :the number of genes in the second pathway
#' @param rho1 :the correlation coefficient for the first pathway
#' @param rho2 :the correlation coefficient for the second pathway
#' @return X	n by p matrix of gene expressions
#' @examples
#' ## generate 6 gene expressions from 10 individuals
#' X.pathway_discrete.balanced(n=10,p=6,q1=2,q2=2,rho1=0.5,rho2=0.5)
#' @import compound.Cox
#' @references Emura T, Chen YH, Chen HY (2012). Survival Prediction Based on Compound Covariate under Cox Proportional Hazard Models. PLoS ONE 7(10): e47627. doi:10.1371/journal.pone.0047627
#' @export
X.pathway_discrete.balanced = function (n, p, q1, q2, rho1 = 0.5, rho2 = 0.5){
  X=X.pathway(n, p, q1, q2, rho1 , rho2 )
  r1 = 1/(1 + sqrt((1 - rho1)/rho1))
  r2 = 1/(1 + sqrt((1 - rho2)/rho2))
  SD0 = sqrt(3/4)
  SD1 = sqrt((3/4) * (2 * r1^2 - 2 * r1 + 1))
  SD2 = sqrt((3/4) * (2 * r2^2 - 2 * r2 + 1))
  for(i in 1:n){
    for(j in 1:p){
      if( j <= q1){
        xrange=seq(-1.5/SD1,1.5/SD1,0.05)
        quant25 = quantile(xrange)[2];quant50 = quantile(xrange)[3];quant75 = quantile(xrange)[4]
        if(X[i,j] < quant25){X[i,j]=1}
        else if(X[i,j] >= quant25 && X[i,j] < quant50){X[i,j]=2}
        else if(X[i,j] >= quant50 && X[i,j] < quant75){X[i,j]=3}
        else{X[i,j]=4}
      }
      else if( j >q1 && j<=q1+q2){
        xrange=seq(-1.5/SD2,1.5/SD2,0.05)
        quant25 = quantile(xrange)[2];quant50 = quantile(xrange)[3];quant75 = quantile(xrange)[4]
        if(X[i,j] < quant25){X[i,j]=1}
        else if(X[i,j] >= quant25 && X[i,j] < quant50){X[i,j]=2}
        else if(X[i,j] >= quant50 && X[i,j] < quant75){X[i,j]=3}
        else{X[i,j]=4}
      }
      else{
        xrange=seq(-1.5/SD0,1.5/SD0,0.05)
        quant25 = quantile(xrange)[2];quant50 = quantile(xrange)[3];quant75 = quantile(xrange)[4]
        if(X[i,j] < quant25){X[i,j]=1}
        else if(X[i,j] >= quant25 && X[i,j] < quant50){X[i,j]=2}
        else if(X[i,j] >= quant50 && X[i,j] < quant75){X[i,j]=3}
        else{X[i,j]=4}
      }
    }
  }
  return(X)
}
