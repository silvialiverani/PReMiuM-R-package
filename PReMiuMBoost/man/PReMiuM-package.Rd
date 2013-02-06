\name{PReMiuM-package}
\alias{PReMiuM-package}
\alias{PReMiuM}
\docType{package}
\title{
Dirichlet Process Bayesian Clustering
}
\description{
Dirichlet process Bayesian clustering and functions for the post-processing of its output. 
}
\section{Details}{
  \pkg{PReMiuM} provides the following:
  \itemize{
\item Implements an infinite Dirichlet process model
\item Can do dependent or independent slice sampling (Kalli et al., 2011) 
     or truncated Dirichlet process model (Ishwaran and James, 2001)
\item Handles categorical or Normal covariates, or a mixture of them
\item Handles Bernoulli, Binomial, Categorical, Poisson or Normal responses
\item Handles inclusion of fixed effects in the response model
\item Handles Extra Variation in the response (for Bernoulli, Binomial and Poisson response only)
\item Handles variable selection (tested in Discrete covariate case only)
\item Includes label switching moves for better mixing
\item Allows user to exclude the response from the model
\item Allows user to compute the entropy of the allocation
\item Allows user to run with a fixed alpha or update alpha (default)
\item Allows users to run predictive scenarios (at C++ run time)
\item Basic or Rao-Blackwellised predictions can be produced
\item Handling of missing data
\item C++ for model fitting
\item Uses Eigen Linear Algebra Library and Boost C++
\item Completely self contained (all library code in included in distribution)
\item Adaptive MCMC where appropriate
\item R package for generating simulation data and post processing
\item R plotting functions allow user choice of what to order clusters by
}
}
\details{
\tabular{ll}{
Package: \tab PReMiuM\cr
Type: \tab Package\cr
Version: \tab 3.0.12\cr
Date: \tab 2013-01-29\cr
License: \tab GPL3\cr
LazyLoad: \tab yes\cr
}

Program to implement Dirichlet Process Bayesian Clustering as described in
Liverani et al. 2013. Previously this project was called profile regression.

}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\section{Acknowledgements}{
Silvia Liverani thanks The Leverhulme Trust for financial support. 
}
\references{
Molitor J, Papathomas M, Jerrett M and Richardson S. (2010) Bayesian Profile Regression with an Application to the National Survey of Children's Health, Biostatistics 11: 484-498.

Papathomas M, Molitor J, Richardson S. et al (2011) Examining the joint effect of multiple risk factors using exposure risk profiles: lung cancer in non smokers. Environmental Health Perspectives 119: 84-91.

Hastie, D. I., Liverani, S., Richardson, S. and Stucker I. (2013) A semi-parametric approach to estimate risk functions associated with multi-dimensional exposure profiles: application to smoking and lung cancer. \emph{Submitted}.

Molitor, J., Brown, I. J., Papathomas, M., Molitor, N., Liverani, S., Chan, Q., Richardson, S., Van Horn, L., Daviglus, M. L., Stamler, J. and Elliott, P. (2013) Blood pressure differences associated with DASH-like lower sodium compared with typical American higher sodium nutrient profile: INTERMAP USA. \emph{Submitted}.

Hastie, D. I., Liverani, S. and Richardson, S. (2013) Sampling from Dirichlet process mixture models with unknown concentration parameter. \emph{Submitted}.

Liverani, S., Hastie, D. I. and Richardson, S. (2013) PReMiuM: An R package for Profile Regression Mixture Models using Dirichlet Processes. \emph{Submitted}.

}

\examples{
\dontrun{
# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=10, nClusInit=20,
    nBurn=20, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")
}
}
