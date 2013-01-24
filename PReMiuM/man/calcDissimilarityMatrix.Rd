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
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Profile Regression Mixture Models using Dirichlet Processes. \emph{Submitted}.

}
\examples{
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=10, nBurn=20, data=inputs$inputData, output="output", 
    covNames=inputs$covNames,nClusInit=15)

dissimObj<-calcDissimilarityMatrix(runInfoObj)

}
\keyword{postprocessing}