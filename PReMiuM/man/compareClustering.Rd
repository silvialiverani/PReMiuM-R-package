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

# simulate a dataset
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)

# do clustering a first time
runInfoObjA<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)
dissimObjA<-calcDissimilarityMatrix(runInfoObjA)
clusObjA<-calcOptimalClustering(dissimObjA)
riskProfileObjA<-calcAvgRiskAndProfile(clusObjA)

# do clustering a second time
runInfoObjB<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)
dissimObjB<-calcDissimilarityMatrix(runInfoObjB)
clusObjB<-calcOptimalClustering(dissimObjB)
riskProfileObjB<-calcAvgRiskAndProfile(clusObjB)

# compare clustering
compareClustering(riskProfileObjA,riskProfileObjB)
}


