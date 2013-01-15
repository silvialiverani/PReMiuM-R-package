\name{summariseVarSelectRho}
\alias{summariseVarSelectRho}
\title{summariseVarSelectRho}
\description{This function summarises the posterior distribution of rho, a parameter for variable selection only.}
\usage{
summariseVarSelectRho(runInfoObj)
}
\arguments{
\item{runInfoObj}{Object of type runInfoObj}
}
\value{A list with the following elements.
\item{rho}{A matrix that has as many columns as the number of covariates and as many rows as the number of sweeps. This matrix records the samples from the posterior distribution of rho for each covariate at each sweep.}
\item{rhoMean}{Vector with the column means of the matrix rho above. Each value corresponds to the posterior mean of rho for each covariate.}
\item{rhoMedian}{Vector with the column medians of the matrix rho above. Each value corresponds to the posterior median of rho for each covariate.}
\item{rhoLowerCI}{Vector with the column lower confidence intervals of the matrix rho above. Each value corresponds to the lower confidence interval of the posterior distribution of rho for each covariate.}
\item{rhoUpperCI}{Vector with the column upper confidence intervals of the matrix rho above. Each value corresponds to the upper confidence interval of the posterior distribution of rho for each covariate.}
}
\author{David Hastie and Silvia Liverani}
\examples{

inputs <- generateSampleDataFile(clusSummaryVarSelectBernoulliDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=100, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, varSelect="Continuous")

rho<-summariseVarSelectRho(runInfoObj)
}
\keyword{variableSelection}

