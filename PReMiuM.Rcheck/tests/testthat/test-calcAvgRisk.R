test_that("Dissimilarity function is running correctly", {
  library(PReMiuM)
  library(testthat)
  set.seed(1234)
  generateDataList <- clusSummaryBernoulliDiscrete()
  inputs <- generateSampleDataFile(generateDataList)
  runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
                       nSweeps=1, nBurn=0, data=inputs$inputData, 
                       output=paste(tempdir(),"/output",sep=""), 
                       covNames=inputs$covNames,nClusInit=15, seed=12345)
  dissimObj<-calcDissimilarityMatrix(runInfoObj)
  expect_equal(head(dissimObj$disSimMat,10),c(1,0, 1, 0, 1, 0, 1, 1, 1, 1))
  expect_equal(tail(dissimObj$disSimMat,10),c(1, 1, 1, 1, 1, 1, 1, 0, 1, 1))
  expect_equal(dissimObj$disSimMat[20304:20306],c(1,1,0))
  expect_equal(length(dissimObj$disSimMat),499500)
  clusObj<-calcOptimalClustering(dissimObj)
  expect_equal(clusObj$clusterSizes[1:3],c(99,100,64))
})

