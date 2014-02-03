\name{setHyperparams}
\alias{setHyperparams}
\title{Definition of characteristics of sample datasets for profile regression}
\description{
Hyperparameters for the priors can be specified here and passed as an argument to profRegr. 

The user can specify some or all hyperparameters.
Those hyperparameters not specified will take their default values. Where the
file is not provided, all hyperparameters will take their default values.
}
\usage{
setHyperparams(shapeAlpha=NULL,rateAlpha=NULL,
     useReciprocalNCatsPhi=NULL,aPhi=NULL,mu0=NULL,Tau0=NULL,R0=NULL,
     kapp0=NULL,muTheta=NULL,sigmaTheta=NULL,dofTheta=NULL,muBeta=NULL,
     sigmaBeta=NULL,dofBeta=NULL,shapeTauEpsilon=NULL,
     rateTauEpsilon=NULL,aRho=NULL,bRho=NULL,shapeSigmaSqY=NULL,
     scaleSigmaSqY=NULL,rSlice=NULL,truncationEps=NULL)
}
\arguments{

\item{shapeAlpha}{The shape parameter for Gamma prior on alpha (default=2)}
\item{rateAlpha}{The inverse-scale (rate) parameter for the Gamma prior on alpha (default=1)}
\item{useReciprocalNCatsPhi}{Boolean denoting whether the vector phi_j (for covariate j) have all elements equal (only used in the discrete covariate case, default=true)}
\item{aPhi}{The vector of parameters for the Dirichlet prior on phi_j. Element j corresponds to covariate j which then has a prior Dirichlet(aPhi[j],aPhi[j],....,aPhi[j]). (Only used in discrete case if useReciprocalNCatsPhi is false, default=(1 1 1 ... 1))}
\item{mu0}{The mean vector for mu_c in the Normal covariate case (only used in Normal covariate case, default=empirical covariate means)}
\item{Tau0}{The precision matrix for mu_c in the Normal covariate case (only used in Normal covariate case, default=inverse of diagonal matrix with elements equal to squareof empirical range for each covariate)}
\item{R0}{The matrix parameter for the Wishart distribution for Tau_c (only used in Normal covariate case, default=1/nCovariates * inverse of empirical covariance matrix)}
\item{kapp0}{The degrees of freedom parameter for the Wishart distribution for Tau_c (only used in Normal covariate case, default=nCovariates).}
\item{muTheta}{The location parameter for the t-Distribution for theta_c (only used if response included in model, default=0)}
\item{sigmaTheta}{The scale parameter for the t-Distribution for theta_c (only used if response included in model, default=2.5)}
\item{dofTheta}{The degrees of freedom parameter for the t-Distribution for theta_c (only used if response included in model, default=7)}
\item{muBeta}{The location parameter for the t-Distribution for beta (only used when fixed effects present, default=0)}
\item{sigmaBeta}{The scale parameter for the t-Distribution for beta (only used when fixed effects present, default=2.5)}
\item{dofBeta}{The dof parameter for the t-Distribution for beta (only used when fixed effects present, default=7)}
\item{shapeTauEpsilon}{Shape parameter for gamma distribution for prior for precision tau of extra variation errors epsilon (only used if extra variation is used i.e. extraYVar argument is included, default=5.0)}
\item{rateTauEpsilon}{Inverse-scale (rate) parameter for gamma distribution for prior for precision tau of extra variation errors epsilon (only used if extra variation is used i.e. extraYVar argument is used, default=0.5)}
\item{aRho}{Parameter for beta distribution for prior on rho in variable selection (default=0.5)}
\item{bRho}{Parameter for beta distribution for prior on rho in variable selection (default=0.5)}
\item{shapeSigmaSqY}{Shape parameter of inverse-gamma prior for sigma_Y^2 (only used in the Normal response model, default =2.5)}
\item{scaleSigmaSqY}{Scale parameter of inverse-gamma prior for sigma_Y^2 (only used in the Normal response model, default =2.5)}
\item{rSlice}{Slice parameter for independent slice sampler such that xi_c = (1-rSlice)*rSlice^c for c=0,1,2,... (only used for slice independent sampler i.e. sampler=SliceIndependent, default 0.75).}
\item{truncationEps}{Parameter for determining the truncation level of the finite Dirichlet process (only used for truncated sampler i.e. sampler=Truncated}
}

\value{
The output of this function is a list with the components defined as above.
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Liverani, S., Hastie, D. I., Azizi, L., Papathomas, M. and Richardson, S. (2013) PReMiuM: An R package for Profile Regression Mixture Models using Dirichlet Processes. \emph{Submitted}. Available at http://uk.arxiv.org/abs/1303.2836

}
\examples{

hyp <- setHyperparams(shapeAlpha=3,rateAlpha=2,mu0=c(30,13),R0=3.2*diag(2))

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=2, nClusInit=15,
    nBurn=2, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames,
    hyper=hyp)

}
\keyword{hyperparameters}
