\name{calcAvgRiskAndProfile}
\alias{calcAvgRiskAndProfile}
\title{Calculation of the average risks and profiles}
\description{Calculation of the average risks and profiles.}
\usage{
calcAvgRiskAndProfile(clusObj, includeFixedEffects=F,
    proportionalHazards=F)
}
\arguments{
\item{clusObj}{Object of type clusObj.}
\item{includeFixedEffects}{By default this is set to FALSE. If it is set to FALSE then the risk profile is computed with the parameters beta of the fixed effects assumed equal to zero. If it is set to TRUE, then risk profile at each sweep is computed adjusting for the sample of the beta parameter at that sweep.}
\item{proportionalHazards}{Whether the risk matrix should include lambda only for the yModel="Survival" case so that the proportional hazards can be computed in the plotting function. The default is the average survival time.}
}
\value{
A list with the following components. This is an object of type riskProfileObj.
\item{riskProfClusObj}{The object of type clusObj as given in the input of this function.}
\item{risk}{A matrix that has a column for each cluster and a row for each sweep. Each element of the matrix represents the estimated risk at each sweep for each cluster.}
\item{profile}{An array whose first dimension is the number of sweeps, the second is the number of clusters, the third is the number of discrete covariates and the fourth is the number of categories of each of the covariates. Each element of the array represents the covariate profile at each sweep for each cluster. The fourth dimension does not exists if the covariate type is Normal. If the covariate type is mixed, then instead of this element, the two elements below are defined, 'profilePhi' and 'profileMu'.}
\item{profileStar}{This is NULL if there has not been any variable selection. otherwise it contains the }
\item{empiricals}{A vector of length of the optimal number of clusters, where each value is the empirical mean of the outcome for each cluster.}
\item{profileStdDev}{An array whose first dimension is the number of sweeps, the second is the number of clusters, the third and the fourth are the number of continuous covariates. Each square matrix identified by the first and second dimension of the array represents the standard deviation at each sweep for each cluster. This element is only available if the covariate type is continuous or mixed.}
\item{profilePhi}{This array is the equivalent of the 'profile' above for discrete covariates in case of mixed covariates.}
\item{profileStarPhi}{This array is defined as profile and profilePhi, but the values are computed only if a variable selection procedure has been run. The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.}
\item{profileMu}{This array is the equivalent of the 'profile' above for Normal covariates in case of mixed covariates.}
\item{profileStarMu}{This array is defined as profile and profileMu, but the values are computed only if a variable selection procedure has been run. The definition of the star profile is given in Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Bayesian profile regression.}
\item{nuArray}{For yModel=Survival when weibullFixedShape=FALSE this array contains the sampled values of the shape parameter nu. The first dimension is the number of sweeps, the second is the number of clusters.}
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
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=10, 
    nBurn=20, data=inputs$inputData, output="output", nClusInit=15,
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
}
}
\keyword{postprocessing}
