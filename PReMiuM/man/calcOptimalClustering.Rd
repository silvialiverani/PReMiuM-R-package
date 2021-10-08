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
\item{useLS}{This is set to FALSE by default. If it is set to TRUE then the least-squares method is used for the calculation of the optimal clustering, as described in Molitor et al (2010). Note that this is set to TRUE by default if disSimObj$onlyLS is set to TRUE.}
}
\value{
the output is a list with the following elements. This is an object of type clusObj.
\item{clusObjRunInfoObj}{Details on this run. An object of type runInfoObj.}
\item{clusterSizes}{Cluster sizes.}
\item{clusteringPred}{The predicted cluster memberships for the predicted scenarios.}
\item{clusObjDisSimMat}{Dissimilarity matrix.}
\item{clustering}{Cluster memberships.}
\item{nClusters}{Optimal number of clusters.}
\item{avgSilhouetteWidth}{Average silhouette width when using medoids method for clustering.}
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}
\examples{
\dontrun{
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=10, nBurn=20, data=inputs$inputData, output="output", 
    covNames=inputs$covNames, nClusInit=15)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
}
}
\keyword{postprocessing}
