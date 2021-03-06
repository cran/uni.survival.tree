% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uni.logrank.R
\name{uni.logrank}
\alias{uni.logrank}
\title{Univariate binary splits by the logrank test}
\usage{
uni.logrank(t.vec, d.vec, X.mat)
}
\arguments{
\item{t.vec}{:Vector of survival times (time to either death or censoring)}

\item{d.vec}{:Vector of censoring indicators (1=death, 0=censoring)}

\item{X.mat}{:n by p matrix of covariates, where n is the sample size and p is the number of covariates}
}
\value{
A dataframe containing:

Pvalue: the P-value of the two-sample logrank test, where the cut-off value is optimized

cut_off_point: the optimal cutt-off values of the binary splits given a feature

left.sample.size: the sample size of a left child node

right.sample.size: the sample size of a right child node
}
\description{
The output is the summary of significance tests for binary splits, where the cut-off values are optimized for each covariate.
}
\details{
The output can be used to construct a logrank tree.
}
\examples{
data(Lung,package="compound.Cox")
train_Lung=Lung[which(Lung[,"train"]==TRUE),] #select training data
t.vec=train_Lung[,1]
d.vec=train_Lung[,2]
x.mat=train_Lung[,-c(1,2,3)]
uni.logrank(t.vec,d.vec,x.mat)
}
