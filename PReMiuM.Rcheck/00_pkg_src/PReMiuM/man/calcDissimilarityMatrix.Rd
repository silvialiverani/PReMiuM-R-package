\name{calcDissimilarityMatrix}
\alias{calcDissimilarityMatrix}
\title{Calculates the dissimilarity matrix}
\description{Calculates the dissimilarity matrix.}
\usage{
calcDissimilarityMatrix(runInfoObj, onlyLS=FALSE)
}
\arguments{
\item{runInfoObj}{Object of type runInfoObj.}
\item{onlyLS}{Logical. It is set to FALSE by default. When it is equal to TRUE the dissimilarity matrix is not returned and the only method available to identify the optimal partition using 'calcOptimalClustering' is least squares. This parameter is to be used for datasets with many subjects, as C++ can compute the dissimilarity matrix but it cannot pass it to R for usage in the function 'calcOptimalClustering'. As guidance, be aware that a dataset with 85,000 subjects will require a RAM of about 26Gb, even if onlyLS=TRUE.}
}
\value{
Need to write this 
\item{disSimRunInfoObj}{These are details regarding the run and in the same format as runInfoObj.}
\item{disSimMat}{The dissimilarity matrix, in vector format. Note that it is diagonal, so this contains the upper triangle diagonal entries.}
\item{disSimMatPred}{The dissimilarity matrix, again in vector format as above, for the predicted subjects.}
\item{lsOptSweep}{The optimal partition among those explored by the MCMC, as defined by the least squares method. See Dahl (2006).}
\item{onlyLS}{Logical. If it set to TRUE the only method available to identify the optimal partition using 'calcOptimalClustering' is least squares.}
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
    covNames=inputs$covNames,nClusInit=15)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
}
}
\keyword{postprocessing}
