prem.bench <- function(classNum = 3){
  
  #create list of function names for each distribution. Will later be called to generate the sample data#
  #
  function.list = c("clusSummaryBernoulliDiscrete", "clusSummaryBernoulliNormal", "clusSummaryBernoulliDiscreteSmall","clusSummaryCategoricalDiscrete","clusSummaryNormalDiscrete","clusSummaryNormalNormal","clusSummaryNormalNormalSpatial","clusSummaryVarSelectBernoulliDiscrete","clusSummaryBernoulliMixed")
  
  #set up the final output variable#
  
  final.out = NULL
  #
  #for loop runs the premium functions on each distribution and outputs the optimal cluster, Y value, and known truth cluster for all elements in each distribution's sample data.
  
  for (i in function.list){
    
    #set seed for consistent results among all distributions#
    #
    set.seed(0001)
    
    #run the function to generate the sample data using the name of the distribution from the function list.
    #
    inputs = do.call("generateSampleDataFile",list(do.call(i,list())))
    
    #directly pull the population size from each generated sample sicne they are different for each distribution
    #
    popSize = dim(inputs$inputData)[1]
    
    #use the population size/number of clusters to create randomly sized clusters that add up to the population total
    #
    clusSizes = c(sample(1:classNum,popSize,replace = TRUE,prob = c(sample(20:80, classNum,replace = TRUE)/100)))
    table(clusSizes)
    
    #create a vector 'known' that repeats "Known 1", "Known 2", etc for the total amount in known cluster 1, known cluster 2, etc.
    #
    known = NULL
    for (j in 1:classNum){
      newdat = rep(paste("Known",j),table(clusSizes)[j])
      known = c(known, newdat)
    }
    
    #Adds a vector of the current distribution to identify the clusters in the final output
    #
    function.vector = c(rep(i,length(known)))
    
    #multiple 'if' statements to determine which parameters need to be filled in. Some distributions have more parameters than others.
    #
    if (exists("fixedEffectsNames",where = inputs)){
      runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, nClusInit=15,nBurn=300, data=inputs$inputData, output="output",covNames = inputs$covNames,fixedEffectsNames = inputs$fixedEffectNames, seed=12345)
    } else {
      if (exists("discreteCovs",where = inputs)){
        runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, nClusInit=15,nBurn=300, data=inputs$inputData, output="output",covNames = inputs$covNames,discreteCovs = inputs$discreteCovs, continuousCovs = inputs$continuousCovs, seed=12345)
      } else {
        runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, nSweeps=100, nClusInit=15,nBurn=300, data=inputs$inputData, output="output",covNames = inputs$covNames, seed=12345)
      }
    }
    
    #the rest of the PReMiuM steps to get the optimal clustering
    #
    dissimObj<-calcDissimilarityMatrix(runInfoObj)
    clusObj<-calcOptimalClustering(dissimObj,maxNClusters =7)
    riskProfileObj<-calcAvgRiskAndProfile(clusObj)
    clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary2.png")
    optAlloc<-clusObj$clustering
    
    #creates a data frame with the optimal clusters, Y values for the simulated data, and the known truth clusters
    #
    tmp_boxplot<-data.frame(opt=as.factor(optAlloc),outcome=inputs$inputData$outcome,known=as.factor(known))
    
    #combine the output of the current distribution with the previous distribution, so the final output is one data frame
    #
    final.out = rbind(final.out,cbind(function.vector,tmp_boxplot))
  }
  return(final.out)
}

#run the function, save output as 'tester'
#
tester = prem.bench()

#Seperate each distribution's output by the name of distribution. Removes the first column (the name of the distributon) to avoid redundancy.
#
BernoulliDiscrete = tester[which(tester == "clusSummaryBernoulliDiscrete"),-1]
BernoulliNormal = tester[which(tester == "clusSummaryBernoulliNormal"),-1]
BernoulliDiscreteSmall = tester[which(tester == "clusSummaryBernoulliDiscreteSmall"),-1]
CategoricalDiscrete = tester[which(tester == "clusSummaryCategoricalDiscrete"),-1]
NormalDiscrete = tester[which(tester == "clusSummaryNormalDiscrete"),-1]
NormalNormal = tester[which(tester == "clusSummaryNormalNormal"),-1]
NormalNormalSpatial = tester[which(tester == "clusSummaryNormalNormalSpatial"),-1]
VarSelectBernoulliDiscrete = tester[which(tester == "clusSummaryVarSelectBernoulliDiscrete"),-1]
BernoulliMixed = tester[which(tester == "clusSummaryBernoulliMixed"),-1]
