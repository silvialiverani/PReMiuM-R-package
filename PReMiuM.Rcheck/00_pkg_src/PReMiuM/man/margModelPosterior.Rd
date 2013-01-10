\name{margModelPosterior}
\alias{margModelPosterior}
\title{Marginal Model Posterior}
\description{Compute the marginal model posterior.}
\usage{
margModelPosterior(runInfoObj)
}
\arguments{
\item{runInfoObj}{boh.}
}
\value{
It returns a number of files in the output directory as well as a list with the following elements.
}
\author{Silvia Liverani}
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

