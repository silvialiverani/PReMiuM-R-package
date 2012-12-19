\name{calcOptimalClustering}
\alias{calcOptimalClustering}
\title{Calculation of the optimal clustering}
\description{Calculates the optimal clustering.}
\usage{
calcOptimalClustering(disSimObj, maxNClusters=NULL, useLS=F)
}
\arguments{
\item{disSimObj}{A dissimilarity matrix (in vector format, as the output of the function calcDissimilarityMatrix(), and as described in ?calcDissimilarityMatrix) or a list of dissimilarity matrix, to combine the output of several runs of the MCMC.}
\item{maxNClusters}{Set the maximum number of clusters allowed. This is set to the maximum number explored.}
\item{useLS}{This is set to FALSE by default. If it is set to TRUE then the least-squares method is used for the calculation of the optimal clustering, as described in Molitor et al (2010).}
}
\value{
the output is a list with the following elements.
\item{clusObjRunInfoObj}{Details on this run.}
\item{clusterSizes}{Cluster sizes.}
\item{clusteringPred}{The predicted cluster memberships for the predicted scenarios.}
\item{clusObjDisSimMat}{Dissimilarity matrix.}
\item{clustering}{Cluster memberships.}
\item{nClusters}{Optimal number of clusters.}
\item{avgSilhouetteWidth}{Not sure}
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

}

