\name{profRegr}
\alias{profRegr}
\title{Profile Regression}
\description{Fit a profile regression model.}
\usage{
profRegr(covNames, fixedEffectsNames, outcome="outcome", 
  outcomeT=NA, data, output="output", hyper, predict, 
	predictType="RaoBlackwell",
  nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, 
  nClusInit, seed, yModel="Bernoulli", xModel="Discrete", 
  sampler="SliceDependent", alpha=-2, dPitmanYor = 0, excludeY=FALSE, 
  extraYVar=FALSE, varSelectType="None", entropy,reportBurnIn=FALSE,
  run=TRUE, discreteCovs, continuousCovs, whichLabelSwitch="123",
  includeCAR=FALSE, neighboursFile="Neighbours.txt", uCARinit=FALSE,
	PoissonCARadaptive=FALSE,	weibullFixedShape=TRUE, 
	useNormInvWishPrior=FALSE, useHyperpriorR1=FALSE)
}
\arguments{
\item{covNames}{A vector of strings of the covariate names as by the column names in the data argument. The names of the covariates cannot include space characters.}
\item{fixedEffectsNames}{A vector of strings of the fixed effect names as by the column names in the data argument. Each fixed effect must be of class 'numeric'. If a fixed effect is of class 'character', an error message will appear and the fixed effect will need to be recoded as numeric. The names of the fixed effects cannot include space characters.}
\item{outcome}{A string of column of the data argument that contains the outcome. The outcome cannot have missing values - you could consider predicting the value of the outcome for those subjects for which it has not been observed. The name cannot include space characters.}
\item{outcomeT}{A string of column of the data argument that contains the offset (for Poisson outcome) or the number of trials (for Binomial outcome) or censoring for Survival reponse (coded as 0 or 1). The name cannot include space characters.}
\item{data}{A data frame which has as columns the outcome, the covariates, the fixed effects if any and the offset (for Poisson outcome) or the number of trials (for Binomial outcome) or censoring (for Survival outcome). The outcome cannot have missing values - you could consider predicting the value of the outcome for those subjects for which it has not been observed. For Survival response censoring must be coded as 0 if the event has not occurred (ie, there has been censoring) and 1 if the event has occurred (no censoring has taken place). The names of the columns cannot include space characters.}
\item{output}{Path to folder to save all output files. The covariates can have missing values, which must be coded as 'NA'. There cannot be missing values in the fixed effects - if there are, use an imputation method before using profile regression.}
\item{hyper}{Object of type setHyperparams with hyperparameters specifications. This is optional, default values are provided for all hyperparameters. See ?setHyperparams for details.}
\item{predict}{Data frame containing the predictive scenarios. This is only required if predictions are requested. 

At each iteration the predictive subjects are assigned to one of the current
clusters according to their covariate profiles (but ignoring missing values), or
their Rao Blackwellised estimate of theta is recorded (a weighted average of all
theta, weighted by the probability of allocation into each cluster. For Normal and Quantile response they can also be randomly allocated. See also the option predictType below.

The predictive subjects have no impact on the likelihood and so do not determine
the clustering or parameters at each iteration. The predictive allocations are
then recorded as extra entries in each row of the output_z.txt file. This can
then be processed in the post processing to create a dissimilarity matrix with
the fitting subjects. The post procesing function calcPredictions will create
predicted response values for these subjects.

See ?calcPredictions for more details and examples. 
}
\item{predictType}{This can be set equal to "RaoBlackwell" and "random". The default is RaoBlackwell. The random option can only be used for Normal and Quantile response, where the estimated variance of the clusters is considered and the predictive subjects are randomly assigned to a mixture component and then are also randomly sampled within that component.}
\item{nSweeps}{Number of iterations of the MCMC after the burn-in period. By default this is 1000.}
\item{nBurn}{Number of initial iterations of the MCMC to be discarded. By default this is 1000.}
\item{reportBurnIn}{If TRUE then the burn in iterations are reported in the output files, if set to FALSE they are not. It is set to FALSE by default.}
\item{nProgress}{The number of sweeps at which to print a progress update. By default this is 500.}
\item{nFilter}{The frequency (in sweeps) with which to write the output to file. The default value is 1.}
\item{nClusInit}{The number of clusters individuals should be initially randomly assigned to (Unif[50,60]).}
\item{seed}{The value for the seed for the random number generator. The default value is the current time.}
\item{yModel}{The model type for the outcome variable. The options currently available are "Bernoulli", "Poisson", "Binomial", "Categorical", "Normal", "Quantile" and "Survival". The default value is Bernoulli.}
\item{xModel}{The model type for the covariates. The options currently available are "Discrete", "Normal" and "Mixed". The default value is "Discrete".}
\item{sampler}{The sampler type to be used. Options are "SliceDependent", "SliceIndependent" and "Truncated".  The default value is "SliceDependent".}
\item{alpha}{The value to be used if alpha is fixed. If a value smaller than or equal to -1 is used then alpha is random, if dPitmanYor is equal to zero (the random alpha option is available for Dirichlet process prior only). The default value is -2 (random alpha). For fixed alpha, if dPitmanYor is in the interval (0,1) then a Pitman-Yor process prior is used instead of a Dirichlet process prior.}
\item{dPitmanYor}{The discount parameter for the Pitman-Yor process prior. The default value is 0, which is equivalent to a Dirichlet process prior. This parameter must belong to the interval [0,1) and it must be provided together with a non-negative value for alpha. The Pitman-Yor process prior is only available for non-random parameters. Note that the third label switching move is only available for Dirichlet process priors, so it will not be run if dPitmanYor>0. Therefore setting dPitmanYor to a value greater than zero will forse whichLabelSwitch=12.}
\item{excludeY}{If TRUE only the covariate data X is modelled. By default this is set to FALSE.}
\item{extraYVar}{If set equal to TRUE extra Gaussian variance is included in the response model. This option is available only for Bernoulli, Binomial and Poisson response. By default the extra Gaussian variance is not included, so extraYVar=FALSE.}
\item{varSelectType}{The type of variable selection to be used "None", "BinaryCluster" or "Continuous".  The "Continuous" variable selection is the implementation of the novel variable selection formulation proposed by Papathomas, Molitor, Hoggart, Hastie, Richardson (2012) "Exploring data from genetic association studies using Bayesian variable selection and the Dirichlet process: application to searching for gene x gene patterns" in Genetic Epidemiology. The "BinaryCluster" variable selection is based on the method proposed by Chung and Dunson (2009) "Nonparametric Bayes conditional distribution modelling with variable selection" in the Journal of the American Statistical Association. Both types of variable selection can be used with discrete, continuous or mixed covariates. The default value is "None".}
\item{entropy}{If included then we compute allocation entropy. By default the allocation entropy is not included.}
\item{run}{Logical. If TRUE then the MCMC is run. Set run=FALSE if the MCMC has been run already and it is only required to collect information about the run.}
\item{discreteCovs}{The names of the discrete covariates among the covariate names, if xModel="Mixed". This and continuousCovs must be defined if xModel="Mixed", while covNames is ignored.}
\item{continuousCovs}{The names of the discrete covariates among the covariate names, if xModel="Mixed". This and continuousCovs must be defined if xModel="Mixed", while covNames is ignored.}
\item{whichLabelSwitch}{The label switching moves to run. The options available are moves 1, 2 and 3 ("123"), moves 1 and 2 ("12") and move 3 only ("3"). The moves are described in Hastie et al. (2013). Note that the third label switching move is only available for Dirichlet process priors, so it will not be run if dPitmanYor>0. Therefore setting dPitmanYor to a value greater than zero will forse whichLabelSwitch=12.}
\item{includeCAR}{A boolean specifying wether a conditional autoregressive term should be introduced within the model, to take into account possible spatial correlation within residuals. Only for Poisson and Normal response models.}
\item{neighboursFile}{The file name of the file specifying neighbourhood graph. It should have the same structure than neighbourhood graph files used in the "INLA" package, and can be produced from a nb object of package "spdep", by the function "nb2INLA" of package "spdep". See ?nb2INLA for details. Each file must have at least one neighbour.}
\item{uCARinit}{This parameter gives the possibility of giving initialisation values for the spatial residuals u of the spatial CAR. It is set to FALSE by default (meaning that the spatial residuals are initialised randomly). It can be set alternatively to a vector of values, one for each of the observations available.}
\item{PoissonCARadaptive}{This parameter controls which sampler is used for the parameters of the spatial random effect when the outcome is Poisson. When it is set to TRUE, the adaptive rejection sampler is used. When it is set to FALSE (default) a random walk Metropolis is used. }
\item{weibullFixedShape}{This parameter controls whether the shape parameter of the Weibull distribution (for yModel=Survival only) is a global parameter (fixed) or cluster specific. It is equal to TRUE by default.}
\item{useNormInvWishPrior}{By default this variable equals FALSE. When this variable equals TRUE, the conjugate Normal-inverse-Wishart prior is used rather
than the independant normal and inverse Wishart priors. If this prior is used, variable selection cannot be used as it has not been implemented.}
\item{useHyperpriorR1}{Adds a hyperprior for R1, which is the (global) scale parameter for Tau, which is the precision matrix for xModel=Normal or Mixed. The default for this option is TRUE.}
}

\value{
Once the C++ has completed the output from fitting the regression is stored in a number of text files in the directory specified. Files are produced
containing the MCMC traces for all of the values of interest, along with a log file and files for monitoring the acceptance rates of the adaptive Metropolis Hastings moves.

It returns a number of files in the output directory as well as a list with the following elements. This an object of type runInfoObj. The files that are produced in the output directory are described below. 
\item{directoryPath}{String. Directory path of the output files.}
\item{fileStem}{String. The }
\item{inputFileName}{String. Location and file name of input dataset as created by this function for the C++ routines}
\item{nSweeps}{Integer. The number of sweeps of the MCMC after the burn-in.}
\item{nBurn}{Integer. The number of iterations in the burn-in period of the MCMC.}
\item{reportBurnIn}{Logical. Whether the output of the burn-in report should be included.}
\item{nFilter}{Integer. The frequency (in sweeps) with which to write the output to file.}
\item{nProgress}{The number of sweeps at which to print a progress update.}
\item{nSubjects}{Integer. The number of subjects.}
\item{nPredictSubjects}{Integer. The number of subjects for which to run predictions.}
\item{fullPredictFile}{Logical. It is FALSE by default. It is equal to TRUE if the outcome or the outcome and the fixed effects were included in the dataframe provided in the input predict. If TRUE, the function will have a produced a file ending in "_predictFull.txt" which contains the values of the outcome and fixed effects for the computation of measures of fit in the function calcPredictions.}
\item{covNames}{A vector of strings with the names of the covariates.}
\item{xModel}{String. The model type for the covariates.}
\item{includeResponse}{Logical. If FALSE only the covariate data X is modelled.}
\item{yModel}{String. The model type for the outcome.}
\item{varSelect}{Logical. If FALSE no variable selection is performed.}
\item{varSelectType}{String. It specifies what type of variable selection has been performed, if any.} 
\item{nCovariates}{Integer. The number of covariates.}
\item{nFixedEffects}{Integer. The number of fixed effects.}
\item{nCategoriesY}{Integer. The number of categories of the outcome, if yModel = "Categorical". It is 1 otherwise.}
\item{nCategories}{Vector of integers. The number of categories of each covariate, if xModel = "Discrete". It is 1 otherwise.}
\item{extraYVar}{TRUE if extra Gaussian variance is included in the response model.}
\item{xMat}{A matrix of the covariate data.}
\item{yMat}{A matrix of the outcome data, including the offset if the outcome is Poisson, the number of trials if the outcome is Binomial and 0 or 1 for Survival outcome (1 for censored individuals, 0 otherwise).}
\item{wMat}{A matrix of the fixed effect data.}
\item{whichLabelSwitch}{The label switching moves that have been run. The options available are moves 1, 2 and 3 ("123"), moves 1 and 2 ("12") and move 3 only ("3"). The moves are described in Hastie et al. (2013).}
\item{includeCAR}{Logical. Whether a spatial CAR term is included.}
\item{predictType}{String. Whether a RaoBlackwell or random predictions have been computed.}
\item{weibullFixedShape}{Logical. Whether the shape parameter of the Weibull distribution for the survival response is fixed or cluster specific.}

These are the files produced in the output directory. We refer to Liverani et al. (2015) 

\item{_alpha.txt}{If alpha is random, each row is a draw from a posterior distribution of alpha (including burn in if reportBurnIn=TRUE).}

\item{_beta.txt}{If fixed effects are included, this file provides the draws from the posterior distribution of the beta parameters at each sweep. Each row represents the vector of beta's at each sweep (including burn in if reportBurnIn=TRUE).}

\item{_hyper.txt}{Internal file to communicate between R and C++ the values of the hyperparamters.}

\item{_input_txt}{Internal file to communicate the data between R and C++.}

\item{_log.txt}{This file logs some information about the run, such as what variables were included, which hyperparameters were used, the seed of the random numbers, the acceptance rates of the MCMC moves that were included in the run.}

\item{_logPost.txt}{This file report the logPosterior, the logLikelihood and logPrior for the model fit at each sweep (including burn in if reportBurnIn=TRUE).}

\item{_nClusters.txt}{This file includes the number of clusters at each sweep. Each row represents a sweep (including burn in if reportBurnIn=TRUE) and each element in the rows is the number of clusters per sweep. This includes the number of empty clusters, if any.}

\item{_nMembers.txt}{This file includes the number of observations in each cluster at each sweep. Each row represents a sweep (including burn in if reportBurnIn=TRUE) and each element in the rows is the number of observations in each cluster per sweep. The last number in each row is the total number of observations, computed as the sum of the elements in the row as a check that all observations have been assigned to a cluster.}

\item{_theta.xt}{This file includes the value of theta (cluster specific parameter for the response variable) for each cluster at each sweep. Each row represents a sweep (including burn in if reportBurnIn=TRUE) and each element in the rows is the value of theta for each cluster at that sweep. The thetas provided her are in the same order as the clusters in _nMembers.txt and they are drawn from the prior when they correspond to empty clusters.}

\item{_z.txt}{This file includes the cluster membership for each observation at each sweep. Each row represents a sweep (including burn in if reportBurnIn=TRUE) and each element in the rows is the cluster membership for each of the observations, ordered as they are provided to profRegr in the dataframe.}

There are more files that can be in the output, depending on which options are used in profRegr. The file _mu.txt for example reports the mean for xModel=Normal, _phi.txt reports the multinomial probabilities for xModel=Discrete, _rho.txt reports the paramters for variable selection, etc. The files usually report one line for each sweep (including burn in if reportBurnIn=TRUE). See Liverani et al. (2015) for more details of the parameters. 

Note that for the _gamma.txt for variable selection the results are reported per sweep (each line is a sweep) and within each line by cluster (so for each covariate the switches per cluster are reported in order, before the second covariate is reported for each cluster, etc).

}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Aurore J. Lavigne, Department of Epidemiology and Biostatistics, Imperial College London, UK

Lamiae Azizi, MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>

The R package PReMiuM is supported through research grants. One key requirement of such funding applications is the ability to demonstrate the impact of the work we seek funding for can. Whatever you are using PReMiuM for, it would be very helpful for us to learn about our users, to tailor our future methodological developments to your needs. Please email us at liveranis@gmail.com or visit http://www.silvialiverani.com/support-premium/. 

}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. URL http://www.jstatsoft.org/v64/i07/.

Hastie, D. I., Liverani, S. and Richardson, S. (2014) Sampling from Dirichlet process mixture models with unknown concentration parameter: Mixing issues in large data implementations. \emph{Forthcoming in the Statistics \& Computing}. Available at http://link.springer.com/article/10.1007%2Fs11222-014-9471-3

}
\examples{
# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=10, nClusInit=20,
    nBurn=20, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)


# example with Bernoulli outcome and Mixed covariates
inputs <- generateSampleDataFile(clusSummaryBernoulliMixed())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=10, nClusInit=15,
    nBurn=20, data=inputs$inputData, output="output", 
    discreteCovs = inputs$discreteCovs,
    continuousCovs = inputs$continuousCovs)

}
\keyword{profileRegression}

