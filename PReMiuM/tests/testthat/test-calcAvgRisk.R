context("Main profile regression function")

test_that("MCMC output of main profile regression function runs correctly", {

  set.seed(1234)
  generateDataList <- clusSummaryBernoulliDiscrete()
  inputs <- generateSampleDataFile(generateDataList)
  runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
                       nSweeps=10, nBurn=20, data=inputs$inputData, output="output", 
                       covNames=inputs$covNames,nClusInit=15, seed=12345)
  
  dissimObj<-calcDissimilarityMatrix(runInfoObj)
    
  expect_equal(dissimObj$disSimMat[3], 0.8)
})
