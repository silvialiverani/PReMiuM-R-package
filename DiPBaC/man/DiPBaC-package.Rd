\name{DiPBaC-package}
\alias{DiPBaC-package}
\alias{DiPBaC}
\docType{package}
\title{
Dirichlet Process Bayesian Clustering
}
\description{
Dirichlet Process Bayesian Clustering and functions for the post-processing of its output. 
}
\details{
\tabular{ll}{
Package: \tab DiPBaC\cr
Type: \tab Package\cr
Version: \tab 2.0.10\cr
Date: \tab 2012-12-13\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
David Hastie and Silvia Liverani

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{
Hastie et al (2012) blah blah
}
\examples{

# generation of data for clustering
## generation of fixed effects
fe1<-rnorm(200,0,1)
fe2<-runif(200,0,1)
## generation of the outcome 
beta<-c(2,3)
W <- cbind(fe1,fe2) %*% beta
theta<- c(-7,0,3)
clusterIndex<-c(rep(1,80),rep(2,60),rep(3,60))
mu<-theta[clusterIndex]+W
p<-1/(1+exp(-mu))
outcome<-vector()
for (i in 1:200){
    if(runif(1)<p[i]){
        outcome[i]<-1
    }else{
        outcome[i]<-0
    }
}
## generation of the covariates
covariateProbs<-list(list(c(0.8,0.1,0.1),
    c(0.8,0.1,0.1),
    c(0.8,0.1,0.1),
    c(0.8,0.1,0.1),
    c(0.8,0.1,0.1)),
    list(c(0.1,0.8,0.1),
    c(0.1,0.8,0.1),
    c(0.1,0.8,0.1),
    c(0.1,0.8,0.1),
    c(0.8,0.1,0.1)),
    list(c(0.8,0.1,0.1),
    c(0.1,0.1,0.8),
    c(0.1,0.1,0.8),
    c(0.1,0.1,0.8),
    c(0.1,0.1,0.8)))
X<-data.frame(Var1=rep(NA,200),Var2=rep(NA,200),
    Var3=rep(NA,200),Var4=rep(NA,200),Var5=rep(NA,200))
for (i in 1:200){
    for (j in 1:5){
        u<-runif(1)
        for(kk in 1:3){
            if(u<cumsum(covariateProbs[[clusterIndex[i]]][[j]])[kk]){
                X[i,j]<-kk-1
                break
            }
        }
    }	
}

inputData<-data.frame(cbind(outcome,X,fe1,fe2))

runInfoObj<-profRegr(yModel="Bernoulli", xModel="Discrete", nSweeps=100, 
    nBurn=100, data=inputData, output="output", 
    covNames=c("Var1","Var2","Var3","Var4","Var5"),
    fixedEffectsNames=c("fe1","fe2"))

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))


}
