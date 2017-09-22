##install packages##

install.packages("plyr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("PReMiuM")

##load packages
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(PReMiuM)

## read.csv() will read in your data if it is a csv, if it is not a .csv use read.table() or the correct variation##
subin=read.csv("Subset_of_Final_Input.csv")

##puts quotes around brand names to make them strings##
subin$brand_hybrid<-as.character(subin$brand_hybrid)

##pulls the variables of class 'numeric' and 'categorical' from the data, excluding yield because its the Y variable##
numericVars <- which(sapply(subin, class)=='numeric' & names(subin) != 'Yield')
categoricalVars <- which(sapply(subin, class)=='character' & names(subin) != 'Yield')

##OPTIONAL: create a data frame with your predicted values if you have values for what you expect for what your covariates are##
preds <- data.frame(matrix(c("Beck-XL 5939AMXT",290,0,98000,78,1,1,"Beck-XL 5939AMXT",290,0,98000,78,1,1), ncol = 7, byrow = TRUE))
colnames(preds) <- names(c(categoricalVars,numericVars))

##defines covariate names for the profile regression function##
covName = names(c(categoricalVars,numericVars))

##run profile regression, use '? profRegr' to see which inputs to use##
runInfoObj <- profRegr(covName, outcome = 'Yield', 
                  yModel = 'Normal', xModel = "Mixed",
                  discreteCovs = c(names(subin[categoricalVars])),
                  continuousCovs = c(names(subin[numericVars])),
                  data = subin, predict = preds)

##calculate dissimilarity matrix##
calcDists <- calcDissimilarityMatrix(runInfoObj)

##find the optimal clustering#
clusts <- calcOptimalClustering(calcDists)


#The values of the boxplots are computed running this line:
#riskProfileObj<-calcAvgRiskAndProfile(clusObj)
#riskProfileObj$risk gives you the values used for the boxplots for the outcome
#riskProfileObj$profile gives you the values used for the boxplots of the covariates

riskProfileObj <- calcAvgRiskAndProfile(clusts)

clusterOrderObj<-plotRiskProfile(riskProfileObj,"boxplots.png")

optAlloc <- clusts$clustering

##create data frame for optimal allocation with outcomes and brand##
tmp_boxplot<-data.frame(opt=as.factor(optAlloc),outcome=subin$Yield,known=as.factor(subin$brand_hybrid))

##open *.jpeg device to put plots in##
jpeg('violinplot.jpeg')

##create violin plot##
p <- ggplot(tmp_boxplot, aes(x=known, y=outcome, fill=opt)) + geom_violin()+ labs(title="",x="Brand Hybrids", y = "Yield") + facet_grid(~known,scales='free',space='free') + guides(fill=guide_legend(title="Clusters")) + theme(strip.text.x = element_blank(), strip.background = element_blank(),axis.text.x = element_text(angle=60,hjust=1)) 

##set color palette##
p+scale_fill_brewer(palette="Set1")

##dev.off() or it will corrupt the files##
dev.off()

##calculate the predictions)
predictions <- calcPredictions(riskProfileObj, fullSweepPredictions = TRUE, fullSweepLogOR = TRUE)

##sets up function to plot prediction density curves##
plotPredictions= function (outfile, runInfoObj, predictions, logOR = FALSE) 
{
  nPredictedSubjects = NULL
  directoryPath = NULL
  fileStem = NULL
  yModel = NULL
  xModel = NULL
  logOddsRatio = NULL
  weibullFixedShape = NULL
  for (i in 1:length(runInfoObj)) assign(names(runInfoObj)[i], 
                                         runInfoObj[[i]])
  if (yModel != "Bernoulli" && yModel != "Normal" && yModel != 
      "Survival" && yModel != "Quantile") 
    stop("This function has been developed for Bernoulli, Normal, Quantile and Survival response only.")
   if (yModel == "Normal" || yModel == "Quantile") 
    logOR <- FALSE
  predictResponseFileName = file.path(runInfoObj$directoryPath, 
                                      paste(runInfoObj$fileStem, "_predict.txt", sep = ""))
  relScenarios <- read.table(predictResponseFileName, header = FALSE, 
                             skip = 1)
  if (logOR == FALSE) {
    preds <- predictions$predictedYPerSweep[, , 1]
  }
  else {
    if (!is.null(predictions$logORPerSweep)) {
      preds <- predictions$logORPerSweep
    }
    else {
      stop("Log OR (odds ratios) cannot be plotted because they have not been computed by calcPredictions. Re-run calcPredictions with option fullSweepLogOR=TRUE.")
    }
  }
  pdf(outfile, onefile = TRUE)
  nPredictSubjects <- runInfoObj$nPredictSubjects
  denObj <- vector(mode = "list")
  for (i in 1:nPredictSubjects) {
    denObj[[i]] <- density(na.omit(preds[, i]))
  }
  for (k in 1:nPredictSubjects) {
    plotDF <- data.frame(logOddsRatio = denObj[[k]]$x, density = denObj[[k]]$y)
    plotObj <- ggplot(plotDF)
    plotObj <- plotObj + geom_line(aes(x = logOddsRatio, 
                                       y = density), size = 0.2)
    plotObj <- plotObj + theme(legend.position = "none")
    plotObj <- plotObj + labs(x = ifelse(logOR == TRUE, "Log OR of response", 
                                         "Response")) + theme(axis.title.x = element_text(size = 15)) + 
      labs(y = "Density") + theme(axis.title.y = element_text(size = 15, 
                                                              angle = 90))
    plotObj <- plotObj + theme(axis.text.x = element_text(size = 15)) + 
      theme(axis.text.y = element_text(size = 15))
    print(plotObj)
  }
  dev.off()
}

##outputs file as a pdf after plotting density curve##
plotPredictions(outfile = "predictiveDensity.pdf", runInfoObj = runInfoObj,predictions = predictions)

##output files will be where##