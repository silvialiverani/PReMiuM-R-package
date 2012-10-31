\name{plotClustering}
\alias{plotClustering}
\title{Plot the clustering using principal components}
\description{Plot the clustering using principal components.}
\usage{
plotClustering(clusObj,outFile,clusterPlotOrder=NULL,whichCovariates=NULL)
}
\arguments{
\item{clusObj}{Object of type clusObj.}
\item{outFile}{The path and filename of the output, a plot saved in PNG format.}
\item{clusterPlotOrder}{don't know.}
\item{whichCovariates}{A vector of integers to select the covariates to print.}
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
\keyword{plots}
