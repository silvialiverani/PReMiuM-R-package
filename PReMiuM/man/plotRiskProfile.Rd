\name{plotRiskProfile}
\alias{plotRiskProfile}
\title{Plot the Risk Profiles}
\description{Plots the risk profiles for a profile regression model.}
\usage{
plotRiskProfile(riskProfObj, outFile, showRelativeRisk=F, 
    orderBy=NULL, whichClusters=NULL, 
    whichCovariates=NULL, useProfileStar=F,riskLim=NULL)
}
\arguments{
\item{riskProfObj}{An object of type riskProfObj.}
\item{outFile}{Path and file name to save the plot.}
\item{showRelativeRisk}{Whether to show the relative risk (with respect to the risk of the first cluster). This option is not available for Normal outcome. For Survival outcomes it computed proportional hazards, but only if the option proportionalHazards=T was used in the function calcAvgRiskAndProfile().} 
\item{orderBy}{Order by which the clusters are to be displayed. It can take values "Empirical", "ClusterSize" and "Risk" (the latter only if the outcome is provided). It can also take the name of a covariate to order the clusters, in which case the clusters are ordered.}
\item{whichClusters}{Either a vector of indeces that corresponds to the clusters that are to be displayed. The length of this vector must be greater than 1. The default is that all clusters are shown.}
\item{whichCovariates}{Either a vector of indeces or a vector of strings that corresponds to the covariates that are to be displayed. The length of this vector must be greater than 1. The default is that all covariates are shown.}
\item{useProfileStar}{To be set equal to TRUE only if a variable selection procedure has been run. The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.}
\item{riskLim}{Limits of the y-axis for the plot of the boxplots for the response variable. The default is NULL. If the riskLim are provided, they should be a vector of length 2. }
}
\value{
This function creates a png plot saved in the path given by outFile. All clusters are visually displayed together.

For discrete covariates, instead of plotting the probability that a phi is above
or below the mean value, we plot the actual phi values (and plot the mean value
across clusters as a horizontal line).

For normal covariates, for each covariate the upper plot is the posterior
distribution for the mean mu, and the lower plot is the posterior distribution
of sqrt(Sigma[j,j]) (i.e. the standard deviation for that covariate).

It also returns the following vector.
\item{meanSortIndex}{This vector is the index that represents the order that the clusters are represented. The default ordering is by empirical risk.}
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. URL http://www.jstatsoft.org/v64/i07/.

}
\examples{
\dontrun{
# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=10, nClusInit=15,
    nBurn=20, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")
}
}
\keyword{plots}
\keyword{postprocessing}
