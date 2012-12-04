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
\item{disSimRunInfoObj}{These are details regarding the run and in the same format as runInfoObj.}
\item{disSimMat}{The dissimilarity matrix, in vector format. Note that it is diagonal, so this contains the upper triangle diagonal entries.}
\item{disSimMatPred}{The dissimilarity matrix, again in vector format as above, for the predicted subjects.}
\item{lsOptSweep}{The optimal partition among those explored by the MCMC, as defined by the least squares method. See Dahl (2006).}
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

