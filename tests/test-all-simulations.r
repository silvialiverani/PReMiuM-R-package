#######################################
# test devised to check if all simulations were producing the same result for old and new code
#######################################

## test code for 
whichPackage<-"new" # or "old" 

library(PReMiuM)
set.seed(123)

# example for Poisson outcome and Discrete covariates
simulations<-c("clusSummaryBernoulliDiscrete",
  "clusSummaryBernoulliDiscrete",
  "clusSummaryBernoulliNormal",
  "clusSummaryBernoulliDiscreteSmall",
  "clusSummaryBinomialNormal",
  "clusSummaryCategoricalDiscrete",
  "clusSummaryNormalDiscrete",
  "clusSummaryNormalNormal",
  "clusSummaryNormalNormalSpatial",
  "clusSummaryPoissonDiscrete",
  "clusSummaryPoissonNormal",
  "clusSummaryPoissonNormalSpatial",
  "clusSummaryVarSelectBernoulliDiscrete",
  "clusSummaryBernoulliMixed",
  "clusSummaryWeibullDiscrete",
  "clusSummaryQuantileNormal")


for (i in 1:length(simulations)){
  inputs <- generateSampleDataFile(get(simulations[i])())
  if (inputs$yModel=="Poisson"||inputs$yModel=="Binomial"||inputs$yModel=="Survival") {
      runInfoObj<-profRegr(yModel=inputs$yModel, 
                     xModel=inputs$xModel, nSweeps=100, nClusInit=15,
                     nBurn=100, data=inputs$inputData, output=paste("output",i,whichPackage,sep=""), 
                     covNames = inputs$covNames, outcomeT = inputs$outcomeT,
                     fixedEffectsNames = inputs$fixedEffectNames,seed=123)
  } else {
    if (inputs$xModel=="Mixed"){
      runInfoObj<-profRegr(yModel=inputs$yModel, 
                         xModel=inputs$xModel, nSweeps=100, nClusInit=15,
                         nBurn=100, data=inputs$inputData, output=paste("output",i,whichPackage,sep=""),
                         discreteCovs = inputs$discreteCovs, continuousCovs=inputs$continuousCovs, 
                         fixedEffectsNames = inputs$fixedEffectNames,seed=123)
    } else {
      runInfoObj<-profRegr(yModel=inputs$yModel, 
                           xModel=inputs$xModel, nSweeps=100, nClusInit=15,
                           nBurn=100, data=inputs$inputData, output=paste("output",i,whichPackage,sep=""),
                           covNames = inputs$covNames, 
                           fixedEffectsNames = inputs$fixedEffectNames,seed=123)
    }
  }
  dissimObj<-calcDissimilarityMatrix(runInfoObj)
  clusObj<-calcOptimalClustering(dissimObj,maxNClusters = 8)
  riskProfileObj<-calcAvgRiskAndProfile(clusObj)
  clusterOrderObj<-plotRiskProfile(riskProfileObj,paste("summary",i,whichPackage,".png",sep=""))
}

#testing excludeY=TRUE
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
runInfoObj<-profRegr(yModel=inputs$yModel, 
                     xModel=inputs$xModel, nSweeps=100, nClusInit=15,
                     nBurn=100, data=inputs$inputData, output=paste("output",i+1,"new",sep=""),
                     covNames = inputs$covNames, #excludeY=TRUE, 
                     fixedEffectsNames = inputs$fixedEffectNames,seed=123)
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj,maxNClusters = 8)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,paste("summary",i+1,whichPackage,".png",sep=""))


