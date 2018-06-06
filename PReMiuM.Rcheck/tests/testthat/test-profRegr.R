context("Main profile regression function")

test_that("MCMC output of main profile regression function runs correctly", {
  library(PReMiuM)
  library(testthat)
  set.seed(1234)
  inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
  runInfoObj<-profRegr(yModel=inputs$yModel, 
                       xModel=inputs$xModel, nSweeps=1, nClusInit=20,
                       nBurn=0, data=inputs$inputData, output="output", 
                       covNames = inputs$covNames, outcomeT = inputs$outcomeT,
                       fixedEffectsNames = inputs$fixedEffectNames,seed=12345)
  testthis<-read.table("output_nClusters.txt")[1]
  expect_equal(testthis[1,1], 21)
})
