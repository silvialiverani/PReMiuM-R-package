\name{generateSampleDataFile}
\alias{generateSampleDataFile}
\title{Generate sample data files for profile regression}
\description{.}
\usage{
generateSampleDataFile(clusterSummary)
}
\arguments{
\item{clusterSummary}{A vector of strings of the covariate names as by the column names in the data argument.}
}
\value{
The output of this function is a list with the following elements
\item{yModel}{The outcome model according to which the data has been generated.}
\item{xModel}{The covariate model according to which the data has been generated.}
\item{inputData}{The data.frame that contains the data.}
\item{covNames}{The names of the covariates.}
\item{fixedEffectNames}{The names of the fixed effects.}

}
\author{David Hastie and Silvia Liverani}
\examples{
# generation of data for clustering

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)


}
\keyword{simulation}

