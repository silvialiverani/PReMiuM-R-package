\name{heatDissMat}
\alias{heatDissMat}
\title{Plot the heatmap of the dissimilarity matrix}
\description{Function to plot the heatmap of the dissimilarity matrix}
\usage{
heatDissMat(dissimObj, main=NULL, xlab=NULL, ylab=NULL)
}
\arguments{
\item{dissimObj}{An object of class dissimObj.}
\item{main}{The usual plot option, to be passed to the heatmap function.}
\item{ylab}{The usual plot option, to be passed to the heatmap function.}
\item{xlab}{The usual plot option, to be passed to the heatmap function.}
}
\value{
Plot of the heatmap of the dissimilary matrix. This functions uses the function 'heatmap' of package 'stats'. Note that this function has not been optimised for large datasets. 
}
\section{Authors}{

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}
\examples{
\dontrun{
# generate simulated dataset
generateDataList <- clusSummaryBernoulliDiscreteSmall()
inputs <- generateSampleDataFile(generateDataList)

# run profile regression
runInfoObj<-profRegr(yModel=inputs$yModel, xModel=inputs$xModel, 
 nSweeps=10, nBurn=2000, data=inputs$inputData, output="output", 
 covNames=inputs$covNames,nClusInit=15)

# compute dissimilarity matrix     
dissimObj<-calcDissimilarityMatrix(runInfoObj)

# plot heatmap
heatDissMat(dissimObj)
}
}


