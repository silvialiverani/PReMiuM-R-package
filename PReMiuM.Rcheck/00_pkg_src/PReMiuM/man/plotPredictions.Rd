\name{plotPredictions}
\alias{plotPredictions}
\title{Plot the conditional density using the predicted scenarios}
\description{Plots the conditional density for the predicted scenarios provided. It produces a pdf with a page for each predictive scenario provided. Each page has a plot of the predicted response, in the order as they were provided to the function. Note that fixed effects are not processed in this function. This function has been developed for Bernoulli, Normal and Survival response only. This function has been developed for Discrete and Normal covariates only.}
\usage{
plotPredictions(outfile, runInfoObj, predictions,
    logOR=FALSE)
}
\arguments{
\item{outfile}{String. The name of the output PDF file. The default is "condDensity.pdf".}
\item{runInfoObj}{An object of type runInfoObj which contains all the details about the run of profRegr.}
\item{predictions}{An object of type predictions which contains all the details about the run of calcPredictions.}
\item{logOR}{Whether to plot the response probability or log odds ratios. The default is FALSE and the response probability is plotted. }
}

\value{
The output is a plot in PDF format.
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
# example with Bernoulli outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
# prediction profiles
preds<-data.frame(matrix(c(
2, 2, 2, 2, 2,
0, 0, NA, 0, 0),ncol=5,byrow=TRUE))

colnames(preds)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
 nSweeps=10000, nBurn=10000, data=inputs$inputData, output="output", 
 covNames=inputs$covNames,predict=preds,
 fixedEffectsNames = inputs$fixedEffectNames)        
dissimObj <- calcDissimilarityMatrix(runInfoObj)
clusObj <- calcOptimalClustering(dissimObj)
riskProfileObj <- calcAvgRiskAndProfile(clusObj)
predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE,fullSweepLogOR=TRUE)

plotPredictions(outfile="predictiveDensity.pdf",runInfoObj=runInfoObj,
 predictions=predictions,logOR=TRUE)

}
}
\keyword{predictions}
