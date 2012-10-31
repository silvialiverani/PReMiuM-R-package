\name{calcPredictions}
\alias{calcPredictions}
\title{Calculates the predictions}
\description{Calculates the predictions.}
\usage{
calcPredictions(riskProfObj,predictResponseFileName=NULL,
    doRaoBlackwell=F,fullSweepPredictions=F,fullSweepLogOR=F)
}
\arguments{
\item{riskProfObj}{Object of type riskProfObj.}
\item{predictResponseFileName}{The file which contains the predictive scenarios.}
\item{doRaoBlackwell}{By default this is set to FALSE.}
\item{fullSweepPredictions}{By default this is set to FALSE.}
\item{fullSweepLogOR}{By default this is set to FALSE.}
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

