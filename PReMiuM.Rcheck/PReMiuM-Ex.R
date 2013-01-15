pkgname <- "PReMiuM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PReMiuM')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PReMiuM-package")
### * PReMiuM-package

flush(stderr()); flush(stdout())

### Name: PReMiuM-package
### Title: Dirichlet Process Bayesian Clustering
### Aliases: PReMiuM-package PReMiuM

### ** Examples


# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=1000, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")




cleanEx()
nameEx("calcAvgRiskAndProfile")
### * calcAvgRiskAndProfile

flush(stderr()); flush(stdout())

### Name: calcAvgRiskAndProfile
### Title: Calculation of the average risks and profiles
### Aliases: calcAvgRiskAndProfile
### Keywords: postprocessing

### ** Examples


generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)



cleanEx()
nameEx("calcDissimilarityMatrix")
### * calcDissimilarityMatrix

flush(stderr()); flush(stdout())

### Name: calcDissimilarityMatrix
### Title: Calculates the dissimilarity matrix
### Aliases: calcDissimilarityMatrix
### Keywords: postprocessing

### ** Examples

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))




cleanEx()
nameEx("calcOptimalClustering")
### * calcOptimalClustering

flush(stderr()); flush(stdout())

### Name: calcOptimalClustering
### Title: Calculation of the optimal clustering
### Aliases: calcOptimalClustering
### Keywords: postprocessing

### ** Examples

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)




cleanEx()
nameEx("calcPredictions")
### * calcPredictions

flush(stderr()); flush(stdout())

### Name: calcPredictions
### Title: Calculates the predictions
### Aliases: calcPredictions
### Keywords: predictions

### ** Examples

inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
     
# prediction profiles
preds<-data.frame(matrix(c(0, 0, 1, 0, 0,
0, 0, 1, NA, 0),ncol=5,byrow=TRUE))
colnames(preds)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
     
# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
    nSweeps=100, nBurn=1000, data=inputs$inputData, output="output", 
    covNames=inputs$covNames,predict=preds)
     
# postprocessing
dissimObj <- calcDissimilarityMatrix(runInfoObj)
clusObj <- calcOptimalClustering(dissimObj)
riskProfileObj <- calcAvgRiskAndProfile(clusObj)
clusterOrderObj <- plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)

# example where the fixed effects can be provided for prediction 
# but the observed response is missing 
# (there are 2 fixed effects in this example). 
# in this example we also use the Rao Blackwellised predictions

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())

# prediction profiles
predsPoisson<- data.frame(matrix(c(7, 2.27, -0.66, 1.07, 9, 
     -0.01, -0.18, 0.91, 12, -0.09, -1.76, 1.04, 16, 1.55, 1.20, 0.89,
     10, -1.35, 0.79, 0.95),ncol=5,byrow=TRUE))
colnames(predsPoisson)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=100, 
         nBurn=100, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, outcomeT="outcomeT",
         fixedEffectsNames = inputs$fixedEffectNames,predict=predsPoisson)

# postprocessing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)


# example where both the observed response and fixed effects are present 
#(there are no fixed effects in this example, but 
# these would just be added as columns between the first and last columns). 

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())

# prediction profiles
predsPoisson<- data.frame(matrix(c(NA, 2.27, -0.66, 1.07, NA, 
     -0.01, -0.18, 0.91, NA, -0.09, -1.76, 1.04, NA, 1.55, 1.20, 0.89,
     NA, -1.35, 0.79, 0.95),ncol=5,byrow=TRUE))
colnames(predsPoisson)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=100, 
         nBurn=100, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, outcomeT="outcomeT",
         fixedEffectsNames = inputs$fixedEffectNames,predict=predsPoisson)

# postprocessing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)





cleanEx()
nameEx("clusSummaryBernoulliDiscrete")
### * clusSummaryBernoulliDiscrete

flush(stderr()); flush(stdout())

### Name: clusSummaryBernoulliDiscrete
### Title: Sample datasets for profile regression
### Aliases: clusSummaryBernoulliDiscrete clusSummaryBinomialNormal
###   clusSummaryCategoricalDiscrete clusSummaryNormalDiscrete
###   clusSummaryNormalNormal clusSummaryPoissonDiscrete
###   clusSummaryPoissonNormal clusSummaryVarSelectBernoulliDiscrete
###   clusSummaryBernoulliMixed
### Keywords: simulation

### ** Examples


clusSummaryBernoulliDiscrete()




cleanEx()
nameEx("generateSampleDataFile")
### * generateSampleDataFile

flush(stderr()); flush(stdout())

### Name: generateSampleDataFile
### Title: Generate sample data files for profile regression
### Aliases: generateSampleDataFile
### Keywords: simulation

### ** Examples

# generation of data for clustering

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)





cleanEx()
nameEx("is.wholenumber")
### * is.wholenumber

flush(stderr()); flush(stdout())

### Name: is.wholenumber
### Title: Function to check if a number is a whole number
### Aliases: is.wholenumber

### ** Examples

is.wholenumber(4) # TRUE
is.wholenumber(3.4) # FALSE



cleanEx()
nameEx("margModelPosterior")
### * margModelPosterior

flush(stderr()); flush(stdout())

### Name: margModelPosterior
### Title: Marginal Model Posterior
### Aliases: margModelPosterior
### Keywords: margModelPosterior

### ** Examples

inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=10, 
         nBurn=10, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, 
         fixedEffectsNames = inputs$fixedEffectNames,seed=12567)

margModelPosterior(runInfoObj)




cleanEx()
nameEx("plotRiskProfile")
### * plotRiskProfile

flush(stderr()); flush(stdout())

### Name: plotRiskProfile
### Title: Plot the Risk Profiles
### Aliases: plotRiskProfile
### Keywords: plots postprocessing

### ** Examples

# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=1000, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")



cleanEx()
nameEx("profRegr")
### * profRegr

flush(stderr()); flush(stdout())

### Name: profRegr
### Title: Profile Regression
### Aliases: profRegr
### Keywords: profileRegression

### ** Examples

# example for Poisson outcome and Discrete covariates
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=1000, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")

# example with Bernoulli outcome and Mixed covariates
inputs <- generateSampleDataFile(clusSummaryBernoulliMixed())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=1000, 
    nBurn=1000, data=inputs$inputData, output="output", 
    discreteCovs = inputs$discreteCovs,
    continuousCovs = inputs$continuousCovs, nClusInit=10)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2,4,5))




cleanEx()
nameEx("setHyperparams")
### * setHyperparams

flush(stderr()); flush(stdout())

### Name: setHyperparams
### Title: Definition of characteristics of sample datasets for profile
###   regression
### Aliases: setHyperparams
### Keywords: hyperparameters

### ** Examples


hyp <- setHyperparams(shapeAlpha=3,rateAlpha=2,mu0=c(30,13),R0=3.2*diag(2))

inputs <- generateSampleDataFile(clusSummaryPoissonNormal())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=2, 
    nBurn=2, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames,
    hyper=hyp)




cleanEx()
nameEx("summariseVarSelectRho")
### * summariseVarSelectRho

flush(stderr()); flush(stdout())

### Name: summariseVarSelectRho
### Title: summariseVarSelectRho
### Aliases: summariseVarSelectRho
### Keywords: variableSelection

### ** Examples


inputs <- generateSampleDataFile(clusSummaryVarSelectBernoulliDiscrete())

runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=100, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, varSelect="Continuous")

rho<-summariseVarSelectRho(runInfoObj)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
