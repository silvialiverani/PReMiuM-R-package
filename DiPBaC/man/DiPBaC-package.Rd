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

# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=1000, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")

# example with Bernoulli outcome and Mixed covariates
inputs <- generateSampleDataFile(clusSummaryBernoulliMixed())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=1000, 
    nBurn=1000, data=inputs$inputData, output="output", 
    discreteCovs = inputs$discreteCovs,
    continuousCovs = inputs$continuousCovs, nClusInit=10)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2,4,5))


}
