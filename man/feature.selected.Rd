% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature.selected.R
\name{feature.selected}
\alias{feature.selected}
\title{The names of features that are selected in a tree}
\usage{
feature.selected(tree)
}
\arguments{
\item{tree}{:an object made from the "uni.tree" function}
}
\value{
An array of characters that are the names from those covariates selected in the tree
}
\description{
 The function returns the names of features (covariates) that are selected as the internal nodes of a tree.
 Only the names of the covariates are shown by excluding the cutt-off values.  
}
\details{
 The outputs show important features for predicting survival outcomes. 
}
\examples{
data(Lung,package="compound.Cox")
train_Lung=Lung[which(Lung[,"train"]==TRUE),] #select training data
t.vec=train_Lung[,1]
d.vec=train_Lung[,2]
x.mat=train_Lung[,-c(1,2,3)]
res=uni.tree(t.vec,d.vec,x.mat,P.value=0.01,d0=0.01,S.plot=FALSE,score=TRUE)
feature.selected(res)
}
