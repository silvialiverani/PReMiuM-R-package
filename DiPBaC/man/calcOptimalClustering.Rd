\name{calcOptimalClustering}
\alias{calcOptimalClustering}
\title{Calculation of the optimal clustering}
\description{Calculates the optimal clustering.}
\usage{
calcOptimalClustering(disSimObj, maxNClusters=NULL, useLS=F)
}
\arguments{
\item{disSimObj}{An object of type disSimObj.}
\item{maxNClusters}{Set the maximum number of clusters allowed. This is set to the maximum number explored.}
\item{useLS}{This is set to FALSE by default. If it is set to TRUE then the least-squares method is used for the calculation of the optimal clustering, as described in Molitor et al (2010).}
}
\value{
Need to write this 
\item{none}{it creates files outside}
}
\author{David Hastie and Silvia Liverani}
\examples{
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

runInfoObj<-readRunInfo(fileStem="output",directoryPath=".")
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))
}

