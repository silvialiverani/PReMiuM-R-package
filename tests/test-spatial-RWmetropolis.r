library(PReMiuM)
set.seed(123)
# generate data
inputs<-generateSampleDataFile(clusSummaryPoissonNormalSpatial())
# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, 
                     xModel=inputs$xModel, nSweeps=2000, nClusInit=15,
                     nBurn=2000, data=inputs$inputData, output="output",
                     covNames = inputs$covNames, outcomeT = "outcomeT", 
                     includeCAR = TRUE, 
                     fixedEffectsNames = inputs$fixedEffectNames,seed=123,
                     PoissonCARadaptive=FALSE)

uCAR<-read.table("output_uCAR.txt")
# check if true value of uCAR random effects are in the credible intervals
pdf(onefile = TRUE)
for (i in 1:100){
  plot(uCAR[,i],type="l")
  abline(h=inputs$uCAR[i],col="red",lwd=2)
  abline(h=mean(uCAR[,i]),col="blue",lwd=2)
  abline(h=quantile(uCAR[,i],0.05),col="blue",lwd=2,lty=2)
  abline(h=quantile(uCAR[,i],0.95),col="blue",lwd=2,lty=2)
}
dev.off()

# what % of credible intervals includes the true value of the uCAR coefficients?
correct<-0
for (i in 1:100){
  if (inputs$uCAR[i]<quantile(uCAR[,i],0.95)&&inputs$uCAR[i]>quantile(uCAR[,i],0.05))
    correct<-correct+1
}
correct/100

# run postprocessing
dissimObj<-calcDissimilarityMatrix(runInfoObj)
clusObj<-calcOptimalClustering(dissimObj,maxNClusters = 8)
riskProfileObj<-calcAvgRiskAndProfile(clusObj)
clusterOrderObj<-plotRiskProfile(riskProfileObj,"summary.png")

