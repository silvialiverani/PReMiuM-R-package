\name{simBenchmark}
\alias{simBenchmark}
\title{Benchmark for simulated examples}
\description{This function checks the cluster allocation of profile regression against the generating clusters for a selection of the simulated dataset provided within the package. This can be used to compute confusion matrices for simulated examples, as shown in the example below.}
\usage{
simBenchmark(whichModel = "clusSummaryBernoulliDiscrete",
  nSweeps = 1000, nBurn = 1000, seedProfRegr = 123)
}
\arguments{
\item{whichModel}{Which simulated dataset this benchmark should be carried out for. At the moment this function works only for these datasets structures provided in the package: "clusSummaryBernoulliNormal","clusSummaryBernoulliDiscreteSmall","clusSummaryCategoricalDiscrete", "clusSummaryNormalDiscrete","clusSummaryNormalNormal", "clusSummaryNormalNormalSpatial",
"clusSummaryVarSelectBernoulliDiscrete", "clusSummaryBernoulliMixed". These dataset structures can be used by the function generateSampleDataFile to create simulated datasets.}
\item{nSweeps}{The number of sweeps of the profile regression algorithm for this benchmarking.}
\item{nBurn}{The number of sweeps in the burn in of the profile regression algorithm for this benchmarking.} 
\item{seedProfRegr}{Sets the seed for the random number generation in profile regression (ie. sets the seed for the portion of the MCMC code in C++). Note that setting this seed does not mean that the function simBenchmark will give the same answer. This is because the first step of this function generates a random sample, which will vary in each run unless a global seed is set in R using the function set.seed.}
}
\value{
This function creates a data.frame. Each row corresponds to each observation in the generated dataset. The columns are: 
\item{clusterAllocation}{Cluster allocation carried out by profile regression. These values are integers, corresponding to cluster numbers.}
\item{outcome}{Value of the outcome (y) in the dataset.}
\item{generatingCluster}{Cluster allocation in the data generating mechanism. These values are characters which include the word 'Known' and then the original numbering of the cluster. The word 'Known' is included to avoid confusion with the cluster allocations identified by profile regression. }
}
\section{Authors}{
Silvia Liverani, Queen Mary University of London, UK

Austin Gratton, University of North Carolina Wilmington, USA

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. URL http://www.jstatsoft.org/v64/i07/.

}
\examples{
\dontrun{
# vector of all test datasets allowed by this bencmarking function
testDatasets<-c("clusSummaryBernoulliNormal",
  "clusSummaryBernoulliDiscreteSmall","clusSummaryCategoricalDiscrete",
  "clusSummaryNormalDiscrete","clusSummaryNormalNormal", 
  "clusSummaryNormalNormalSpatial","clusSummaryVarSelectBernoulliDiscrete", 
  "clusSummaryBernoulliMixed")

# runs profile regression on all datasets and 
# computes confusion matrix for each one
for (i in 1:length(testDatasets)){
  tester<-simBenchmark(testDatasets[i])
  print(table(tester[,c(1,3)]))
}


}
}
\keyword{benchmark}
\keyword{simulations}
\keyword{testing}
