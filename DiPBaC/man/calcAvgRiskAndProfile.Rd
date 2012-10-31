\name{calcAvgRiskAndProfile}
\alias{calcAvgRiskAndProfile}
\title{Calculation of the average risks and profiles}
\description{Calculation of the average risks and profiles.}
\usage{
calcAvgRiskAndProfile(clusObj, includeFixedEffects=F)
}
\arguments{
\item{clusObj}{Object of type clusObj.}
\item{includeFixedEffects}{By default this is set to 'FALSE'.}
}
\value{
Need to write this 
\item{none}{dont know}
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

