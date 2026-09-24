\name{generateSampleDataFile}
\alias{generateSampleDataFile}
\title{Generate sample data files for profile regression}
\description{Generation of random sample datasets for profile regression.}
\usage{
generateSampleDataFile(clusterSummary, pQuantile=0.05)
}
\arguments{
\item{clusterSummary}{A vector of strings of the covariate names as by the column names in the data argument.}
\item{pQuantile}{pQuantile is the quantile parameter of the Asymmetric Laplace Distribution used to generate data to test the model for the quantiles.}
}
\value{
The output of this function is a list with the following elements
\item{yModel}{The outcome model according to which the data has been generated.}
\item{xModel}{The covariate model according to which the data has been generated.}
\item{inputData}{The data.frame that contains the data.}
\item{covNames}{The names of the covariates.}
\item{fixedEffectNames}{The names of the fixed effects.}
\item{uCAR}{The spatial gaussian effect. It is sample into the intrinsic autoregressive model with precision TauCAR under the constraint that the sum of term is null. Only used if includeCAR is TRUE.}
\item{TauCAR}{The precision of the spatial CAR effect. Only used if includeCAR is TRUE.}
\item{Permutation}{A vector of size nSubject given the cluster name of each subject. When spatial CAR is added to the model, for preventing potential identifiability problems, the clusters are randomly distributed within the all subjects. Only used if includeCAR is TRUE.}
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Aurore J. Lavigne, Department of Epidemiology and Biostatistics, Imperial College London, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}
\examples{
# generation of data for clustering

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)

}
\keyword{simulation}

