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
nameEx("calcAvgRiskAndProfile")
### * calcAvgRiskAndProfile

flush(stderr()); flush(stdout())

### Name: calcAvgRiskAndProfile
### Title: Calculation of the average risks and profiles
### Aliases: calcAvgRiskAndProfile

### ** Examples


generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))



cleanEx()
nameEx("calcDissimilarityMatrix")
### * calcDissimilarityMatrix

flush(stderr()); flush(stdout())

### Name: calcDissimilarityMatrix
### Title: Calculates the dissimilarity matrix
### Aliases: calcDissimilarityMatrix

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

### ** Examples

inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete())
     
# prediction profiles
preds<-data.frame(matrix(c(0, 0, 1, 0, 0,
0, 0, 1, NA, 0),ncol=5,byrow=TRUE))
colnames(preds)<-names(inputs$inputData)[2:(inputs$nCovariates+1)]
     
# run profiel regression
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



cleanEx()
nameEx("clusSummaryBernoulliDiscrete")
### * clusSummaryBernoulliDiscrete

flush(stderr()); flush(stdout())

### Name: clusSummaryBernoulliDiscrete
### Title: Definition of characteristics of sample datasets for profile
###   regression
### Aliases: clusSummaryBernoulliDiscrete clusSummaryBinomialNormal
###   clusSummaryCategoricalDiscrete clusSummaryNormalDiscrete
###   clusSummaryNormalNormal clusSummaryPoissonDiscrete
###   clusSummaryPoissonNormal clusSummaryVarSelectBernoulliDiscrete
###   clusSummaryBernoulliMixed
### Keywords: simulation

### ** Examples


clusSummaryBernoulliDiscrete()




cleanEx()
nameEx("compareClustering")
### * compareClustering

flush(stderr()); flush(stdout())

### Name: compareClustering
### Title: Function to compare different partitions
### Aliases: compareClustering

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
nameEx("plotClustering")
### * plotClustering

flush(stderr()); flush(stdout())

### Name: plotClustering
### Title: Plot the clustering using principal components
### Aliases: plotClustering
### Keywords: plots

### ** Examples


generateDataList <- clusSummaryBernoulliDiscrete()
inputs <- generateSampleDataFile(generateDataList)
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, 
    nBurn=100, data=inputs$inputData, output="output", 
    covNames=inputs$covNames)

dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png",
    whichCovariates=c(1,2))




cleanEx()
nameEx("plotRiskProfile")
### * plotRiskProfile

flush(stderr()); flush(stdout())

### Name: plotRiskProfile
### Title: Plot the Risk Profiles
### Aliases: plotRiskProfile
### Keywords: plot

### ** Examples

# generation of data for clustering
## generation of fixed effects
fe1<-rnorm(200,0,1)
fe2<-runif(200,0,1)
## generation of the outcome 
beta<-c(2,3)
W <- cbind(fe1,fe2) 
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
