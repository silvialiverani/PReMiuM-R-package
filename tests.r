# test all functions

library(PReMiuM)

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryBernoulliMixed()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryBinomialNormal()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryCategoricalDiscrete()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryNormalDiscrete()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryNormalNormal()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryPoissonDiscrete()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryPoissonNormal()
inputs <- generateSampleDataFile(generateDataList)
generateDataList <- clusSummaryVarSelectBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList) 

generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
# prediction profiles
preds<-data.frame(matrix(c(0, 0, 1, 0, 0,0, 0, 1, NA, 0),ncol=5,byrow=TRUE))
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
         nBurn=1000, data=inputs$inputData, output="output", 
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
         nBurn=1000, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, outcomeT="outcomeT",
         fixedEffectsNames = inputs$fixedEffectNames,
         nClusInit=15, predict=predsPoisson)
# postprocessing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
output_predictions <- calcPredictions(riskProfileObj,fullSweepPredictions=TRUE)
generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=1000, data=inputs$inputData, output="output", nClusInit=15,
    covNames=inputs$covNames,extraYVar=TRUE)
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
computeRatioOfVariance(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj,useLS=T)
is.wholenumber(4) # TRUE
is.wholenumber(3.4) # FALSE
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
         xModel=inputs$xModel, nSweeps=100, 
         nBurn=1000, data=inputs$inputData, output="output", 
         covNames = inputs$covNames, nClusInit=15,
         fixedEffectsNames = inputs$fixedEffectNames)
margModelPosterior(runInfoObj)
inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=100, nClusInit=20,
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames)
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")
hyp <- setHyperparams(shapeAlpha=3,rateAlpha=2,mu0=c(30,13),R0=3.2*diag(2))
inputs <- generateSampleDataFile(clusSummaryPoissonNormal())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=100, nClusInit=15,
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, outcomeT = inputs$outcomeT,
    fixedEffectsNames = inputs$fixedEffectNames,
    hyper=hyp)
inputs <- generateSampleDataFile(clusSummaryVarSelectBernoulliDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
    xModel=inputs$xModel, nSweeps=100, nClusInit=15, 
    nBurn=1000, data=inputs$inputData, output="output", 
    covNames = inputs$covNames, varSelect="Continuous")
rho<-summariseVarSelectRho(runInfoObj)


