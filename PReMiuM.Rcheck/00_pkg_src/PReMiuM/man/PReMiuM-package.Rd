\name{PReMiuM-package}
\alias{PReMiuMpackage}
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
\item Handles Bernoulli, Binomial, Categorical, Poisson, survival or Normal responses
\item Handles inclusion of fixed effects in the response model, including a spatial CAR (conditional autoregressive) term
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
Version: \tab 3.0.33\cr
Date: \tab 2015-03-13\cr
License: \tab GPL2\cr
LazyLoad: \tab yes\cr
}

Program to implement Dirichlet Process Bayesian Clustering as described in
Liverani et al. 2014. This is a package for Bayesian clustering using a Dirichlet process mixture model. This model is an alternative to regression models, non-parametrically linking a response vector to covariate data through cluster membership. The package allows Bernoulli, Binomial, Poisson, Normal, survival and categorical response, as well as Normal and discrete covariates. It also allows for fixed effects in the response model, where a spatial CAR (conditional autoregressive) term can be also included. Additionally, predictions may be made for the response, and missing values for the covariates are handled. Several samplers and label switching moves are implemented along with diagnostic tools to assess convergence. A number of R functions for post-processing of the output are also provided. In addition to fitting mixtures, it may additionally be of interest to determine which covariates actively drive the mixture components. This is implemented in the package as variable selection. 

The R package PReMiuM is supported through research grants. One key requirement of such funding applications is the ability to demonstrate the impact of the work we seek funding for can. Whatever you are using PReMiuM for, it would be very helpful for us to learn about our users, to tailor our future methodological developments to your needs. Please email us at liveranis@gmail.com or visit http://www.silvialiverani.com/support-premium/. 

}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Aurore J. Lavigne, Department of Epidemiology and Biostatistics, Imperial College London, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\section{Acknowledgements}{
Silvia Liverani thanks The Leverhulme Trust for financial support. 

The R package PReMiuM is supported through research grants. One key requirement of such funding applications is the ability to demonstrate the impact of the work we seek funding for can. Whatever you are using PReMiuM for, it would be very helpful for us to learn about our users, to tailor our future methodological developments to your needs. Please email us at liveranis@gmail.com or visit http://www.silvialiverani.com/support-premium/. 

}
\references{
Molitor J, Papathomas M, Jerrett M and Richardson S. (2010) Bayesian Profile Regression with an Application to the National Survey of Children's Health, Biostatistics 11: 484-498.

Papathomas M, Molitor J, Richardson S. et al (2011) Examining the joint effect of multiple risk factors using exposure risk profiles: lung cancer in non smokers. Environmental Health Perspectives 119: 84-91.

Hastie, D. I., Liverani, S., Azizi, L., Richardson, S. and Stucker I. (2013) A semi-parametric approach to estimate risk functions associated with multi-dimensional exposure profiles: application to smoking and lung cancer. \emph{BMC Medical Research Methodology}. 13 (1), 129.

Molitor, J., Brown, I. J., Papathomas, M., Molitor, N., Liverani, S., Chan, Q., Richardson, S., Van Horn, L., Daviglus, M. L., Stamler, J. and Elliott, P. (2014) Blood pressure differences associated with DASH-like lower sodium compared with typical American higher sodium nutrient profile: INTERMAP USA.
Hypertension 64 (6), 1198-1204. Available at http://www.ncbi.nlm.nih.gov/pubmed/25201893

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. URL http://www.jstatsoft.org/v64/i07/.

Hastie, D. I., Liverani, S. and Richardson, S. (2014) Sampling from Dirichlet process mixture models with unknown concentration parameter: Mixing issues in large data implementations. \emph{Forthcoming in the Statistics \& Computing}. Available at http://link.springer.com/article/10.1007%2Fs11222-014-9471-3

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
