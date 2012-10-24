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
Version: \tab 2.0.2\cr
Date: \tab 2012-10-23\cr
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

outcome<-c(2,2,1,1,1,2,1,1,1,1,1,1)
Variable1<-c(4,1,2,0,2,2,rep(4,6))
Variable2<-c(0,0,2,0,2,2,rep(2,6))
Variable3<-c(0,1,0,0,2,2,rep(0,6))
Variable4<-c(0,1,1,0,2,2,rep(0,6))
Variable5<-c(0,1,2,1,1,1,rep(0,6))

inputData<-data.frame(cbind(outcome,Variable1,Variable2,Variable3,Variable4,Variable5))

profRegr(yModel="Bernoulli", xModel="Discrete", nSweeps=100, nBurn=100, data=inputData, output="output", covNames=c("Variable1","Variable2","Variable3","Variable4","Variable5"))

runInfoObj<-readRunInfo(fileStem="output",directoryPath="output")
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"output/summary.png",whichCovariates=c(1,2))
}
