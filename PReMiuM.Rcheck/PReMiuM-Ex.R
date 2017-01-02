pkgname <- "PReMiuM"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PReMiuM')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PReMiuM-package")
### * PReMiuM-package

flush(stderr()); flush(stdout())

### Name: PReMiuM-package
### Title: Dirichlet Process Bayesian Clustering
### Aliases: PReMiuMpackage PReMiuM

### ** Examples

## Not run: 
##D # example for Poisson outcome and Discrete covariates
##D inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
##D runInfoObj<-profRegr(yModel=inputs$yModel, 
##D     xModel=inputs$xModel, nSweeps=10, nClusInit=20,
##D     nBurn=20, data=inputs$inputData, output="output", 
##D     covNames = inputs$covNames, outcomeT = inputs$outcomeT,
##D     fixedEffectsNames = inputs$fixedEffectNames)
##D 
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D clusObj<-calcOptimalClustering(dissimObj)
##D riskProfileObj<-calcAvgRiskAndProfile(clusObj)
##D clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")
## End(Not run)



cleanEx()
nameEx("calcAvgRiskAndProfile")
### * calcAvgRiskAndProfile

flush(stderr()); flush(stdout())

### Name: calcAvgRiskAndProfile
### Title: Calculation of the average risks and profiles
### Aliases: calcAvgRiskAndProfile
### Keywords: postprocessing

### ** Examples

## Not run: 
##D generateDataList <- clusSummaryBernoulliDiscrete()
##D inputs <- generateSampleDataFile(generateDataList)
##D runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=10, 
##D     nBurn=20, data=inputs$inputData, output="output", nClusInit=15,
##D     covNames=inputs$covNames)
##D 
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D clusObj<-calcOptimalClustering(dissimObj)
##D riskProfileObj<-calcAvgRiskAndProfile(clusObj)
## End(Not run)



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
    nSweeps=10, nBurn=20, data=inputs$inputData, output="output", 
    covNames=inputs$covNames,nClusInit=15)

dissimObj<-calcDissimilarityMatrix(runInfoObj)




cleanEx()
nameEx("calcOptimalClustering")
### * calcOptimalClustering

flush(stderr()); flush(stdout())

### Name: calcOptimalClustering
### Title: Calculation of the optimal clustering
### Aliases: calcOptimalClustering
### Keywords: postprocessing

### ** Examples

## Not run: 
##D generateDataList <- clusSummaryBernoulliDiscrete()
##D inputs <- generateSampleDataFile(generateDataList)
##D runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
##D     nSweeps=10, nBurn=20, data=inputs$inputData, output="output", 
##D     covNames=inputs$covNames, nClusInit=15)
##D 
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D clusObj<-calcOptimalClustering(dissimObj)
## End(Not run)



cleanEx()
nameEx("calcPredictions")
### * calcPredictions

flush(stderr()); flush(stdout())

### Name: calcPredictions
### Title: Calculates the predictions
### Aliases: calcPredictions
### Keywords: predictions

### ** Examples

## Not run: 
##D inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
##D      
##D # prediction profiles
##D preds<-data.frame(matrix(c(0, 0, 1, 0, 0,
##D 0, 0, 1, NA, 0),ncol=5,byrow=TRUE))
##D colnames(preds)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
##D      
##D # run profile regression
##D runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
##D     nSweeps=100, nBurn=1000, data=inputs$inputData, output="output", 
##D     covNames=inputs$covNames,predict=preds)
##D      
##D # postprocessing
##D dissimObj <- calcDissimilarityMatrix(runInfoObj)
##D clusObj <- calcOptimalClustering(dissimObj)
##D riskProfileObj <- calcAvgRiskAndProfile(clusObj)
##D clusterOrderObj <- plotRiskProfile(riskProfileObj,"summary.png",
##D     whichCovariates=c(1,2))
##D output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)
##D 
##D # example where the fixed effects can be provided for prediction 
##D # but the observed response is missing 
##D # (there are 2 fixed effects in this example). 
##D # in this example we also use the Rao Blackwellised predictions
##D 
##D inputs <- generateSampleDataFile(clusSummaryPoissonNormal())
##D 
##D # prediction profiles
##D predsPoisson<- data.frame(matrix(c(7, 2.27, -0.66, 1.07, 9, 
##D      -0.01, -0.18, 0.91, 12, -0.09, -1.76, 1.04, 16, 1.55, 1.20, 0.89,
##D      10, -1.35, 0.79, 0.95),ncol=5,byrow=TRUE))
##D colnames(predsPoisson)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
##D 
##D # run profile regression
##D runInfoObj<-profRegr(yModel=inputs$yModel, 
##D          xModel=inputs$xModel, nSweeps=100, 
##D          nBurn=100, data=inputs$inputData, output="output", 
##D          covNames = inputs$covNames, outcomeT="outcomeT",
##D          fixedEffectsNames = inputs$fixedEffectNames,predict=predsPoisson)
##D 
##D # postprocessing
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D clusObj<-calcOptimalClustering(dissimObj)
##D riskProfileObj<-calcAvgRiskAndProfile(clusObj)
##D output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)
##D 
##D 
##D # example where both the observed response and fixed effects are present 
##D #(there are no fixed effects in this example, but 
##D # these would just be added as columns between the first and last columns). 
##D 
##D inputs <- generateSampleDataFile(clusSummaryPoissonNormal())
##D 
##D # prediction profiles
##D predsPoisson<- data.frame(matrix(c(NA, 2.27, -0.66, 1.07, NA, 
##D      -0.01, -0.18, 0.91, NA, -0.09, -1.76, 1.04, NA, 1.55, 1.20, 0.89,
##D      NA, -1.35, 0.79, 0.95),ncol=5,byrow=TRUE))
##D colnames(predsPoisson)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
##D 
##D # run profile regression
##D runInfoObj<-profRegr(yModel=inputs$yModel, 
##D          xModel=inputs$xModel, nSweeps=10, 
##D          nBurn=20, data=inputs$inputData, output="output", 
##D          covNames = inputs$covNames, outcomeT="outcomeT",
##D          fixedEffectsNames = inputs$fixedEffectNames,
##D          nClusInit=15, predict=predsPoisson)
##D 
##D # postprocessing
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D clusObj<-calcOptimalClustering(dissimObj)
##D riskProfileObj<-calcAvgRiskAndProfile(clusObj)
##D output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)
##D 
## End(Not run)



cleanEx()
nameEx("clusSummaryBernoulliDiscrete")
### * clusSummaryBernoulliDiscrete

flush(stderr()); flush(stdout())

### Name: clusSummaryBernoulliDiscrete
### Title: Sample datasets for profile regression
### Aliases: clusSummaryBernoulliDiscrete clusSummaryBernoulliNormal
###   clusSummaryBernoulliDiscreteSmall clusSummaryBinomialNormal
###   clusSummaryCategoricalDiscrete clusSummaryNormalDiscrete
###   clusSummaryNormalNormal clusSummaryNormalNormalSpatial
###   clusSummaryPoissonDiscrete clusSummaryPoissonNormal
###   clusSummaryPoissonNormalSpatial clusSummaryVarSelectBernoulliDiscrete
###   clusSummaryBernoulliMixed clusSummaryWeibullDiscrete
###   clusSummaryQuantileNormal
### Keywords: simulation

### ** Examples

names(clusSummaryBernoulliDiscrete())




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




cleanEx()
nameEx("globalParsTrace")
### * globalParsTrace

flush(stderr()); flush(stdout())

### Name: globalParsTrace
### Title: Plot of the trace of some of the global parameters
### Aliases: globalParsTrace

### ** Examples


# generate simulated dataset
generateDataList <- clusSummaryBernoulliDiscreteSmall()
inputs <- generateSampleDataFile(generateDataList)

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
 nSweeps=10, nBurn=20, data=inputs$inputData, output="output", nFilter=3,
 covNames=inputs$covNames,nClusInit=15,reportBurnIn=FALSE,
 fixedEffectsNames = inputs$fixedEffectNames)

# plot trace for alpha
globalParsTrace(runInfoObj,parameters="alpha",plotBurnIn=FALSE)




cleanEx()
nameEx("heatDissMat")
### * heatDissMat

flush(stderr()); flush(stdout())

### Name: heatDissMat
### Title: Plot the heatmap of the dissimilarity matrix
### Aliases: heatDissMat

### ** Examples

## Not run: 
##D # generate simulated dataset
##D generateDataList <- clusSummaryBernoulliDiscreteSmall()
##D inputs <- generateSampleDataFile(generateDataList)
##D 
##D # run profile regression
##D runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
##D  nSweeps=10, nBurn=2000, data=inputs$inputData, output="output", 
##D  covNames=inputs$covNames,nClusInit=15)
##D 
##D # compute dissimilarity matrix     
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D 
##D # plot heatmap
##D heatDissMat(dissimObj)
## End(Not run)



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
nameEx("mapforGeneratedData")
### * mapforGeneratedData

flush(stderr()); flush(stdout())

### Name: mapforGeneratedData
### Title: Map generated data
### Aliases: mapforGeneratedData

### ** Examples

inputs=generateSampleDataFile(clusSummaryPoissonNormalSpatial())
mapforGeneratedData(inputs$uCAR)



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
         xModel=inputs$xModel, nSweeps=5, 
         nBurn=10, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, nClusInit=15,
         fixedEffectsNames = inputs$fixedEffectNames)

margModelPosterior(runInfoObj)




cleanEx()
nameEx("plotPredictions")
### * plotPredictions

flush(stderr()); flush(stdout())

### Name: plotPredictions
### Title: Plot the conditional density using the predicted scenarios
### Aliases: plotPredictions
### Keywords: predictions, plots

### ** Examples

## Not run: 
##D # example with Bernoulli outcome and Discrete covariates
##D inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
##D # prediction profiles
##D preds<-data.frame(matrix(c(
##D 2, 2, 2, 2, 2,
##D 0, 0, NA, 0, 0),ncol=5,byrow=TRUE))
##D 
##D colnames(preds)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
##D # run profile regression
##D runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
##D  nSweeps=10000, nBurn=10000, data=inputs$inputData, output="output", 
##D  covNames=inputs$covNames,predict=preds,
##D  fixedEffectsNames = inputs$fixedEffectNames)        
##D dissimObj <- calcDissimilarityMatrix(runInfoObj)
##D clusObj <- calcOptimalClustering(dissimObj)
##D riskProfileObj <- calcAvgRiskAndProfile(clusObj)
##D predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE,fullSweepLogOR=TRUE)
##D 
##D plotPredictions(outfile="predictiveDensity.pdf",runInfoObj=runInfoObj,
##D  predictions=predictions,logOR=TRUE)
##D 
## End(Not run)



cleanEx()
nameEx("plotRiskProfile")
### * plotRiskProfile

flush(stderr()); flush(stdout())

### Name: plotRiskProfile
### Title: Plot the Risk Profiles
### Aliases: plotRiskProfile
### Keywords: plots postprocessing

### ** Examples

## Not run: 
##D # example for Poisson outcome and Discrete covariates
##D inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
##D runInfoObj<-profRegr(yModel=inputs$yModel, 
##D     xModel=inputs$xModel, nSweeps=10, nClusInit=15,
##D     nBurn=20, data=inputs$inputData, output="output", 
##D     covNames = inputs$covNames, outcomeT = inputs$outcomeT,
##D     fixedEffectsNames = inputs$fixedEffectNames)
##D 
##D dissimObj<-calcDissimilarityMatrix(runInfoObj)
##D clusObj<-calcOptimalClustering(dissimObj)
##D riskProfileObj<-calcAvgRiskAndProfile(clusObj)
##D clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")
## End(Not run)



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
    xModel=inputs$xModel, nSweeps=2, nClusInit=15,
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
    xModel=inputs$xModel, nSweeps=10, nClusInit=15, 
    nBurn=20, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, varSelect="Continuous")

rho<-summariseVarSelectRho(runInfoObj)



cleanEx()
nameEx("vec2mat")
### * vec2mat

flush(stderr()); flush(stdout())

### Name: vec2mat
### Title: Vector to upper triangular matrix
### Aliases: vec2mat

### ** Examples

vec2mat(data=c(1,2,3),nrow=3)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
