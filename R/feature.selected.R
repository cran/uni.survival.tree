#  feature selection #
#' @title Features which are used to build a tree
#' @description The name of covariates which use to partition the covariate space.
#' @param tree :Tree structure made by getTree() function
#' @return An array of characters which is the name of those covariate used in the tree
#' @details The sequence only contain unique name(Only count one time even the single covariate use to split twice or more)
#' @examples
#' data(Lung,package="compound.Cox")
#' train_Lung=Lung[which(Lung[,"train"]==TRUE),] #select training data
#' t.vec=train_Lung[,1]
#' d.vec=train_Lung[,2]
#' x.mat=train_Lung[,-c(1,2,3)]
#' res=uni.tree(t.vec,d.vec,x.mat,P.value=0.01,d0=0.01,S.plot=FALSE,score=TRUE)
#' feature.selected(res)
#' @export
feature.selected=function(tree){
  #split.covariate function can record the covariate use to split in current node
  split.covariate=function(tree,covariate.seq=NULL){
    if(as.character(tree[[1]]$Information)[1] == "terminal node"){
      return(covariate.seq)
    }else{
      covariate.name=as.character(tree[[1]]$Information)[4]  #[4] is the name which record the name of covariate
      covariate.seq=c(covariate.seq,covariate.name)  #combine new selected covariate
      left_subtree = tree[[2]]  #[[2]] call out the left child
      right_subtree = tree[[3]] #[[3]] call out the right child
      return(c(split.covariate(left_subtree,covariate.seq=covariate.seq),split.covariate(right_subtree,covariate.seq=covariate.seq)))
    }
  }
  return(unique(split.covariate(tree,covariate.seq=NULL)))
}


