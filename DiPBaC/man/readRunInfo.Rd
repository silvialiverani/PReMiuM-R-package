\name{readRunInfo}
\alias{readRunInfo}
\title{Reads the information about a run of profile regression}
\description{Reads the information about a run of profile regression.}
\usage{
readRunInfo(directoryPath, fileStem='output')
}
\arguments{
\item{directoryPath}{The directory path of the profile regression run of interest.}
\item{fileStem}{The beginning of each file of the profile regression run of interest.}
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
}

