\name{margModelPosterior}
\alias{margModelPosterior}
\title{Marginal Model Posterior}
\description{Compute the marginal model posterior.}
\usage{
margModelPosterior(runInfoObj)
}
\arguments{
\item{runInfoObj}{An object of type runInfoObj.}
}
\value{
It returns a file in the output folder, with name ending in "_margModPost.txt", that contains the marginal model posterior. It also returns the mean of the values of the marginal model posterior as they appear in the file ending in "_margModPost.txt" in the output folder.
}
\section{Authors}{
Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Profile Regression Mixture Models using Dirichlet Processes. \emph{Submitted}.

}
\examples{
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=10, 
         nBurn=10, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, 
         fixedEffectsNames = inputs$fixedEffectNames,seed=12567)

margModelPosterior(runInfoObj)

}
\keyword{margModelPosterior}

