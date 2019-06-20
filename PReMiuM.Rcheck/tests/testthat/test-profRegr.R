test_that("MCMC output of main profile regression function runs correctly", {
  library(PReMiuM)
  library(testthat)
  set.seed(1234)
  inputs <- generateSampleDataFile(clusSummaryPoissonDiscrete())
  runInfoObj<-profRegr(yModel=inputs$yModel, 
                       xModel=inputs$xModel, nSweeps=1, nClusInit=20,
                       nBurn=0, data=inputs$inputData, 
                       output=paste(tempdir(),"/output",sep=""), 
                       covNames = inputs$covNames, outcomeT = inputs$outcomeT,
                       fixedEffectsNames = inputs$fixedEffectNames,seed=12345)
  testthis<-read.table(paste(tempdir(),"/output_nClusters.txt",sep=""))[1]
  expect_equal(testthis[1,1], 21)
})
