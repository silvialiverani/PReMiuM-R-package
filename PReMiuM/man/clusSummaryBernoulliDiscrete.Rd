\name{clusSummaryBernoulliDiscrete}
\alias{clusSummaryBernoulliDiscrete}
\alias{clusSummaryBernoulliNormal}
\alias{clusSummaryBernoulliDiscreteSmall}
\alias{clusSummaryBinomialNormal}
\alias{clusSummaryCategoricalDiscrete}
\alias{clusSummaryNormalDiscrete}
\alias{clusSummaryNormalNormal}
\alias{clusSummaryNormalNormalSpatial}
\alias{clusSummaryPoissonDiscrete}
\alias{clusSummaryPoissonNormal}
\alias{clusSummaryPoissonNormalSpatial}
\alias{clusSummaryVarSelectBernoulliDiscrete}
\alias{clusSummaryBernoulliMixed}
\alias{clusSummaryWeibullDiscrete}
\alias{clusSummaryQuantileNormal}
\alias{clusSummaryGammaNormal}
\title{Sample datasets for profile regression}
\description{Definition of skeleton of sample datasets for profile regression.}
\usage{
clusSummaryBernoulliDiscrete()
clusSummaryBernoulliNormal
clusSummaryBernoulliDiscreteSmall()
clusSummaryBinomialNormal()
clusSummaryCategoricalDiscrete()
clusSummaryNormalDiscrete()
clusSummaryNormalNormal()
clusSummaryNormalNormalSpatial()
clusSummaryPoissonDiscrete()
clusSummaryPoissonNormal()
clusSummaryPoissonNormalSpatial()
clusSummaryVarSelectBernoulliDiscrete()
clusSummaryBernoulliMixed()
clusSummaryWeibullDiscrete()
clusSummaryQuantileNormal()
clusSummaryGammaNormal()
}
\section{Details}{

clusSummaryBernoulliDiscrete generates a dataset with Bernoulli outcome and discrete covariates.

clusSummaryBernoulliNormal generates a dataset with Bernoulli outcome and Normal covariates.

clusSummaryBernoulliDiscreteSmall generates a dataset with Bernoulli outcome and discrete covariates (with smaller cluster sizes).

clusSummaryBinomialNormal generates a dataset with Binomial outcome and discrete covariates.

clusSummaryCategoricalDiscrete generates a dataset with categorical outcome and discrete covariates.

clusSummaryNormalDiscrete generates a dataset with Normal outcome and discrete covariates.

clusSummaryNormalNormal generates a dataset with Normal outcome and Normal covariates.

clusSummaryNormalNormalSpatial generates a dataset with Normal outcome, Normal covariates and a spatial conditional autoregressive term in the log relative risk.

clusSummaryPoissonDiscrete generates a dataset with Poisson outcome and discrete covariates.

clusSummaryPoissonNormal generates a dataset with Poisson outcome and Normal covariates.

clusSummaryPoissonNormalSpatial generates a dataset with Poisson outcome, Normal covariates and a spatial conditional autoregressive term in the log relative risk.

clusSummaryVarSelectBernoulliDiscrete generates a dataset with Bernoulli outcome and discrete covariates, suitable for variable selection as some covariates are not driving the clustering.

clusSummaryBernoulliMixed generates a dataset with Bernoulli outcome and mixed covariates.

clusSummaryWeibullDiscrete generates a dataset with a Weibull outcome and censored observations.

clusSummaryQuantileNormal generates a dataset with a Quantile outcome.

clusSummaryGammaNormal generates a dataset with a Gamma outcome with scale=1.
}

\value{
The output of these function is a list with the following components. These can be used as inputs for profile regression function profRegr().

\item{outcomeType}{The outcome type of the dataset.}
\item{covariateType}{The covariate type of the dataset.}
\item{nCovariates}{The number of covariates generated.}
\item{nCategories}{The number of categories of the covariates if the covariates are discrete or mixed.}
\item{nFixedEffects}{The number of fixed effects.}
\item{fixedEffectsCoeffs}{The names of the fixed effects.}
\item{missingDataProb}{The pobability of generating missing data.}
\item{nClusters}{The number of clusters.}
\item{clusterSizes}{The number of observations in each cluster.}
\item{clusterData}{The dataset, including the outcome, the covariates, the fixed effects, the number of trials (if Binomial outcome) and the offset (for Poisson outcome).}
\item{covNames}{The names of the covariates of the dataset.}
\item{nDiscreteCovs}{The number of discrete covariates, if the covariate type is mixed.}
\item{nContinuousCovs}{The number of continuous covariates, if the covariate type is mixed.}
\item{outcomeT}{The name of the column of the dataset containing the number of trials (if Binomial outcome) or the offset (for Poisson outcome).}
\item{includeCAR}{A boolean specifying wether a spatial CAR term is included.}
\item{TauCAR}{The precision for the spatial CAR term.}
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
names(clusSummaryBernoulliDiscrete())

}
\keyword{simulation}

