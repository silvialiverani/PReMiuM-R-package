context("Dissimilarity function")

test_that("Dissimilarity function is running correctly", {
  library(PReMiuM)
  library(testthat)
  set.seed(1234)
  generateDataList <- clusSummaryBernoulliDiscrete()
  inputs <- generateSampleDataFile(generateDataList)
  runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
                       nSweeps=1, nBurn=0, data=inputs$inputData, output="output", 
                       covNames=inputs$covNames,nClusInit=15, seed=12345)
  dissimObj<-calcDissimilarityMatrix(runInfoObj)
  expect_equal(dissimObj$disSimMat[3], 0)
})

