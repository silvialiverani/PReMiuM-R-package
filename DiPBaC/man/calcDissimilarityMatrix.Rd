\name{calcDissimilarityMatrix}
\alias{calcDissimilarityMatrix}
\title{Calculates the dissimilarity matrix}
\description{Calculates the dissimilarity matrix.}
\usage{
calcDissimilarityMatrix(runInfoObj)
}
\arguments{
\item{runInfoObj}{Object of type runInfoObj.}
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

