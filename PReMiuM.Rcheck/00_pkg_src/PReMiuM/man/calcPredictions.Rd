\name{calcPredictions}
\alias{calcPredictions}
\title{Calculates the predictions}
\description{Calculates the predictions.}
\usage{
calcPredictions(riskProfObj, predictResponseFileName=NULL,
    doRaoBlackwell=F, fullSweepPredictions=F, fullSweepLogOR=F)
}
\arguments{
\item{riskProfObj}{Object of type riskProfObj.}
\item{predictResponseFileName}{If this function is run after the function profRegr, and outcome (and possibly fixed effects) are known for the predicted profiles, then there is no need to set this, as the function profRegr will have produced a file ending in "_predictFull.txt". This file allows the computation of measures of fit for cross-validation. If the file has not been produced automatically, it can be produced manually and it can be provided here. We discourage this and we provide no documentation for doing so.}
\item{doRaoBlackwell}{By default this is set to FALSE. If it is set to TRUE then Rao-Blackwell predictions are computed.}
\item{fullSweepPredictions}{By default this is set to FALSE. If it is set to TRUE then a prediction is computed for each sweep.}
\item{fullSweepLogOR}{By default this is set to FALSE. If it is set to TRUE then a prediction log OR is computed for each sweep.}
}
\section{Details}{
This functions computes predicted responses, for various prediction scenarios. It is assumed that the predictive
allocations and Rao-Blackwell predictions have already been done in profRegr using the 'predict' input. 

The user can provide the function profRegr with a data.frame through the predict argument. This data.frame has a row for each subject, where each row contains values for the response, fixed effects and offset / number of trials (depending
on the response model) where available. Missing values in this data.frame are denoted by 'NA'. If the data.frame is not provided then the response, fixed effect and offset data is treated as missing for all subjects. If a subject is missing fixed
effect values, then the mean value or 0 category fixed effect is used in the
predictions (i.e. no fixed effect contribution to predicted response). If the
offset / number of trials is missing this value is taken to be 1 when making
predictions. If the response is provided for all subjects, the predicted
responses are compared with the observed responses and the bias and rmse are
computed.

The function can produce predicted values based on simple allocations
(the default), or a Rao-Blackwellised estimate of predictions, where
the probabilities of allocations are used instead of actually performing a
random allocation.
}

\value{
The output is a list with the following elements.

\item{bias}{The bias of the predicted values with respect to the observed outcome. If the response is not provided, this is set to NA.}
\item{rmse}{The root mean square error of the predicted values with respect to the observed outcome. If the response is not provided, this is set to NA.}
\item{mae}{The mean absolute error of the predicted values with respect to the observed outcome. If the response is not provided, this is set to NA.}
\item{observedY}{The values of the outcome provided by the user. This is in the case that predictions are run as a validation tool. If the response is not provided, this is set to NA.}
\item{predictedY}{This matrix has as many rows as predictions requested by the user. It is the mean of the predicted values over all the sweeps that have been run after the burn-in period.}
\item{doRaoBlackwell}{This is set to TRUE if it has done Rao-Blackwell predictions, and FALSE otherwise.}
\item{predictedYPerSweep}{This array has the first dimension equivalent to the number of sweeps and the second dimension as large as the number of predictions requested by the user. It contains the predicted values per sweep.}
\item{logORPerSweep}{This array has the first dimension equivalent to the number of sweeps and the second dimension as large as the number of predictions requested by the user. It contains the predicted log OR values per sweep (not available for Poisson and Normal outcome).}
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Liverani, S., Hastie, D. I., Papathomas, M. and Richardson, S. (2013) PReMiuM: An R package for Profile Regression Mixture Models using Dirichlet Processes. \emph{Submitted}. Available at http://uk.arxiv.org/abs/1303.2836

}
\examples{
\dontrun{
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
     
# prediction profiles
preds<-data.frame(matrix(c(0, 0, 1, 0, 0,
0, 0, 1, NA, 0),ncol=5,byrow=TRUE))
colnames(preds)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
     
# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=1000, data=inputs$inputData, output="output", 
    covNames=inputs$covNames,predict=preds)
     
# postprocessing
dissimObj <- calcDissimilarityMatrix(runInfoObj)
clusObj <- calcOptimalClustering(dissimObj)
riskProfileObj <- calcAvgRiskAndProfile(clusObj)
clusterOrderObj <- plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)

# example where the fixed effects can be provided for prediction 
# but the observed response is missing 
# (there are 2 fixed effects in this example). 
# in this example we also use the Rao Blackwellised predictions

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())

# prediction profiles
predsPoisson<- data.frame(matrix(c(7, 2.27, -0.66, 1.07, 9, 
     -0.01, -0.18, 0.91, 12, -0.09, -1.76, 1.04, 16, 1.55, 1.20, 0.89,
     10, -1.35, 0.79, 0.95),ncol=5,byrow=TRUE))
colnames(predsPoisson)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=100, 
         nBurn=100, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, outcomeT="outcomeT",
         fixedEffectsNames = inputs$fixedEffectNames,predict=predsPoisson)

# postprocessing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)


# example where both the observed response and fixed effects are present 
#(there are no fixed effects in this example, but 
# these would just be added as columns between the first and last columns). 

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())

# prediction profiles
predsPoisson<- data.frame(matrix(c(NA, 2.27, -0.66, 1.07, NA, 
     -0.01, -0.18, 0.91, NA, -0.09, -1.76, 1.04, NA, 1.55, 1.20, 0.89,
     NA, -1.35, 0.79, 0.95),ncol=5,byrow=TRUE))
colnames(predsPoisson)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=10, 
         nBurn=20, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, outcomeT="outcomeT",
         fixedEffectsNames = inputs$fixedEffectNames,
         nClusInit=15, predict=predsPoisson)

# postprocessing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)

}
}
\keyword{predictions}
