\name{compareClustering}
\alias{compareClustering}
\title{Function to compare different partitions}
\description{Function to compare different partitions.}
\usage{
compareClustering(riskProfileObjA,riskProfileObjB,clusterOrder=NULL)
}
\arguments{
\item{riskProfileObjA}{dont know}
\item{riskProfileObjB}{dont know}
\item{clusterOrder}{dont know}
}
\value{
Need to look into this, I don't know how it compares them.
}
\author{David Hastie and Silvia Liverani}
\examples{
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))

}


