\name{margModelPosterior}
\alias{margModelPosterior}
\title{Marginal Model Posterior}
\description{Compute the marginal model posterior.}
\usage{
margModelPosterior(runInfoObj,allocation)
}
\arguments{
\item{runInfoObj}{An object of type runInfoObj.}
\item{allocation}{By default, if allocation is not provided, the _z.txt file is read to compute the marginal model posterior for all the partitions available there. If allocation is equal to a vector that corresponds to a partition, the marginal model posterior is computed for that given partition.}
}
\value{
It returns a file in the output folder, with name ending in "_margModPost.txt", that contains the marginal model posterior. It also returns a list. The first argument is called margModPost and it is the mean of the values of the marginal model posterior as they appear in the file ending in "_margModPost.txt" in the output folder. The second argument is an updated runInfoObj which also include some hyperparameter values. 
}
\section{Authors}{
Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. URL http://www.jstatsoft.org/v64/i07/.

}
\examples{
\dontrun{
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=5, 
         nBurn=10, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, nClusInit=15,
         fixedEffectsNames = inputs$fixedEffectNames)

margModelPost<-margModelPosterior(runInfoObj)
}
}
\keyword{margModelPosterior}

