\name{calcPredictions}
\alias{calcPredictions}
\title{Calculates the predictions}
\description{Calculates the predictions.}
\usage{
calcPredictions(riskProfObj,
    predictResponseFileName=NULL,
    doRaoBlackwell=F,fullSweepPredictions=F,fullSweepLogOR=F)
}
\arguments{

NOTE: THIS DESCRIPTION HAS NOT BEEN CHANGED FOR R AND THE FUNCTIONALITY HAVE NOT BEEN ADAPTED YET

There is also now an additional function calcPredictions for computing predicted
responses, for various prediction scenarios. It is assumed that the predictive
allocations and Rao-Blackwell predictions have already been done in C++ (using
the --predict=<filename> run time option). The user can provide the function
with a file through the predictResponseFileName argument. This file has the
number of subjects, followed by a row for each subject, where each row contains
values for the response, fixed effects and offset / number of trials (depending
on the response model) where available. Missing values in this file are denoted
(-999). If the file is not provided then the response, fixed effect and offset
data is treated as missing for all subjects. If a subject is missing fixed
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

An example file where the fixed effects can be provided for prediction but the
observed response is missing is data/input/example_Poisson_Normal_predictW.txt
(there are 2 fixed effects in this example). An example of using the function,
which would do the Rao Blackwellised predictions, is given by

>calcPredictions(riskProfileObj,
  predictResponseFileName='../data/input/example_Poisson_Normal_predictW.txt',
  doRaoBlackwell=T)

An example file where both the observed response and fixed effects are present
is in data/input/example_Poisson_Normal_predictYW.txt (there are no fixed
effects in this example, but these would just be added as columns between the
first and last columns). An example of using the function, which would do the
simple predictions using the allocations produced by the C++, is given by

>calcPredictions(riskProfileObj,
  predictResponseFileName='../data/input/example_Poisson_Normal_predictYW.txt',
  doRaoBlackwell=F)


\item{riskProfObj}{Object of type riskProfObj.}
\item{predictResponseFileName}{The file which contains the predictive scenarios, as it is produced by the code and stored.}
\item{doRaoBlackwell}{By default this is set to FALSE.}
\item{fullSweepPredictions}{By default this is set to FALSE.}
\item{fullSweepLogOR}{By default this is set to FALSE.}
}
\value{
Need to write this 
\item{none}{it creates files outside}
}
\author{David Hastie and Silvia Liverani}
\examples{
# generate data
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)

# prediction profiles
preds<-data.frame(matrix(c(0, 0, 1, 0, 0,
0, 0, 1, NA, 0),ncol=5,byrow=TRUE))
colnames(preds)<-names(inputs$inputData)[-1]

# run profiel regression
profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
	nBurn=1000, data=inputs$inputData, output="output", 
	covNames=inputs$covNames,predict=preds)

# postprocessing
runInfoObj<-readRunInfo(fileStem="output",directoryPath=".")
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
	whichCovariates=c(1,2))
calcPredictions(riskProfileObj)
}

