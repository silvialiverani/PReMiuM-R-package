\name{globalParsTrace}
\alias{globalParsTrace}
\title{Plot of the trace of some of the global parameters}
\description{Function to lot the trace of some global parameters}
\usage{
globalParsTrace(runInfoObj, parameters = "nClusters",plotBurnIn=FALSE,whichBeta=1)
}
\arguments{
This function allows to visualise the trace of the global parameters. 

Note that this function has not been optimised for large datasets. 

\item{runInfoObj}{An object of class runInfoObj.}
\item{parameters}{The parameter whose trace will be plotted. This can be set equal to "nClusters" (default), "alpha", "mpp" and "beta", as by the model. As beta can be a vector, we advise to also set the option "whichBeta" below to select which fixed effect parameter to visualise in the plot. "mpp" stands for marginal partition posterior, also referred to as marginal model posterior. }
\item{plotBurnIn}{Set to FALSE (default) it does not plot the trace for the burn in period. Set to TRUE it plots the trace including the burn in period.}
\item{whichBeta}{Integer which selects which fixed effect parameter is plotted.}
}
\value{
Plot of trace of some global parameters.
}
\section{Authors}{

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}
\examples{
\dontrun{
# generate simulated dataset
generateDataList <- clusSummaryBernoulliDiscreteSmall()
inputs <- generateSampleDataFile(generateDataList)

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
 nSweeps=10, nBurn=20, data=inputs$inputData, output="output", nFilter=3,
 covNames=inputs$covNames,nClusInit=15,reportBurnIn=FALSE,
 fixedEffectsNames = inputs$fixedEffectNames)

# plot trace for alpha
globalParsTrace(runInfoObj,parameters="alpha",plotBurnIn=FALSE)
}
}


