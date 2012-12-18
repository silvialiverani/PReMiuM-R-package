\name{profRegr}
\alias{profRegr}
\title{Profile Regression}
\description{Fit a profile regression model.}
\usage{
profRegr(covNames, fixedEffectsNames, outcome="outcome", 
        outcomeT=NA, data, output="output", hyper, predict, 
        nSweeps=1000, nBurn=1000, nProgress=500, nFilter=1, 
        nClusInit, seed, yModel="Bernoulli", xModel="Discrete", 
        sampler="SliceDependent", alpha=-1, excludeY, 
        extraYVar, varSelectType, entropy,reportBurnIn=FALSE,
        run=TRUE, discreteCovs, continuousCovs)
}
\arguments{
\item{covNames}{A vector of strings of the covariate names as by the column names in the data argument.}
\item{fixedEffectsNames}{A vector of strings of the fixed effect names as by the column names in the data argument.}
\item{outcome}{A string of column of the data argument that contains the outcome. The outcome cannot have missing values - you could consider predicting the value of the outcome for those subjects for which it has not been observed.}
\item{outcomeT}{A string of column of the data argument that contains the offset (for Poisson outcome) or the number of trials (for Binomial outcome).}
\item{data}{A data frame which has as columns the outcome, the covariates, the fixed effects if any and the offset (for Poisson outcome) or the number of trials (for Binomial outcome). The outcome cannot have missing values - you could consider predicting the value of the outcome for those subjects for which it has not been observed.}
\item{output}{Path to folder to save all output files. The covariates can have missing values, which must be coded as 'NA'. There cannot be missing values in the fixed effects - if there are, use an imputation method before using profile regression.}
\item{hyper}{Path to file with hyperparameters specifications. This is optional: default values are provided for all hyperparameters.}
\item{predict}{Data frame containing the predictive scenarios. See ?calcPredictions for more details. This is only required if predictions are requested.}
\item{nSweeps}{Number of iterations of the MCMC after the burn-in period. By default this is 1000.}
\item{nBurn}{Number of initial iterations of the MCMC to be discarded. By default this is 1000.}
\item{reportBurnIn}{If TRUE then the burn in iterations are reported in the output files, if set to FALSE they are not. It is set to FALSE by default.}
\item{nProgress}{The number of sweeps at which to print a progress update. By default this is 500.}
\item{nFilter}{The frequency (in sweeps) with which to write the output to file. The default value is 1.}
\item{nClusInit}{The number of clusters individuals should be initially randomly assigned to (Unif[50,60]).}
\item{seed}{The value for the seed for the random number generator. The default value is the current time.}
\item{yModel}{The model type for the outcome variable. The options currently available are "Bernoulli", "Poisson", "Binomial", "Categorical" and "Normal". The default value is Bernoulli.}
\item{xModel}{The model type for the covariates. The options currently available are "Discrete", "Normal" and "Mixed". The default value is "Discrete".}
\item{sampler}{The sampler type to be used. Options are "SliceDependent", "SliceIndependent" and "Truncated".  The default value is "SliceDependent".}
\item{alpha}{The value to be used if alpha is to remain fixed. If a negative value is used then alpha is updated. The default value is -1.}
\item{excludeY}{If included only the covariate data X is modelled. By default this is not included.}
\item{extraYVar}{If included extra Gaussian variance is included in the response model. By default the extra Gaussian variance is not included.}
\item{varSelectType}{The type of variable selection to be used "None", "BinaryCluster" or "Continuous".  The default value is "None".}
\item{entropy}{If included then we compute allocation entropy. By default the allocation entropy is not included.}
\item{run}{Logical. If TRUE then the MCMC is run. Set run=FALSE if the MCMC has been run already and it is only required to collect information about the run.}
\item{discreteCovs}{The names of the discrete covariates among the covariate names, if xModel="Mixed". This and continuousCovs must be defined if xModel="Mixed", while covNames is ignored.}
\item{continuousCovs}{The names of the discrete covariates among the covariate names, if xModel="Mixed". This and continuousCovs must be defined if xModel="Mixed", while covNames is ignored.}
}
\value{
It returns a number of files in the output directory as well as a list with the following elements.
\item{directoryPath}{String. Directory path of the output files.}
\item{fileStem}{String. The }
\item{inputFileName}{String. Location and file name of input dataset as created by this function for the C++ routines}
\item{nSweeps}{Integer. The number of sweeps of the MCMC after the burn-in.}
\item{nBurn}{Integer. The number of iterations in the burn-in period of the MCMC.}
\item{reportBurnIn}{Logical. Whether the output of the burn-in report should be included.}
\item{nFilter}{Integer. The frequency (in sweeps) with which to write the output to file. The default value is 1.}
\item{nSubjects}{Integer. The number of subjects.}
\item{nPredictSubjects}{Integer. The number of subjects for which to run predictions.}
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
\item{xMat}{A matrix of the covariate data.}
\item{yMat}{A matrix of the outcome data, including the offset if the outcome is Poisson and the number of trials if the outcome is Binomial.}
\item{wMat}{A matrix of the fixed effect data.}
}
\author{David Hastie and Silvia Liverani}
\examples{
# generation of data for clustering
## generation of fixed effects
fe1<-rnorm(200,0,1)
fe2<-runif(200,0,1)
## generation of the outcome 
beta<-c(2,3)
W <- cbind(fe1,fe2) %*% beta
theta<- c(-7,0,3)
clusterIndex<-c(rep(1,80),rep(2,60),rep(3,60))
mu<-theta[clusterIndex]+W
p<-1/(1+exp(-mu))
outcome<-vector()
for (i in 1:200){
    if(runif(1)<p[i]){
        outcome[i]<-1
    }else{
        outcome[i]<-0
    }
}
## generation of the covariates
covariateProbs<-list(list(c(0.8,0.1,0.1),
    c(0.8,0.1,0.1),
    c(0.8,0.1,0.1),
    c(0.8,0.1,0.1),
    c(0.8,0.1,0.1)),
    list(c(0.1,0.8,0.1),
    c(0.1,0.8,0.1),
    c(0.1,0.8,0.1),
    c(0.1,0.8,0.1),
    c(0.8,0.1,0.1)),
    list(c(0.8,0.1,0.1),
    c(0.1,0.1,0.8),
    c(0.1,0.1,0.8),
    c(0.1,0.1,0.8),
    c(0.1,0.1,0.8)))
X<-data.frame(Var1=rep(NA,200),Var2=rep(NA,200),
    Var3=rep(NA,200),Var4=rep(NA,200),Var5=rep(NA,200))
for (i in 1:200){
    for (j in 1:5){
        u<-runif(1)
        for(kk in 1:3){
            if(u<cumsum(covariateProbs[[clusterIndex[i]]][[j]])[kk]){
                X[i,j]<-kk-1
                break
            }
        }
    }	
}

inputData<-data.frame(cbind(outcome,X,fe1,fe2))

runInfoObj<-profRegr(yModel="Bernoulli", xModel="Discrete", 
    nSweeps=100, nBurn=100, data=inputData, output="output", 
    covNames=c("Var1","Var2","Var3","Var4","Var5"),
    fixedEffectsNames=c("fe1","fe2"))

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))

}
\keyword{profileRegression}

