\name{computeRatioOfVariance}
\alias{computeRatioOfVariance}
\title{computeRatioOfVariance}
\description{Computes of the ratio between the variance of the extra variation and the total variance.}
\usage{
computeRatioOfVariance(runInfoObj)
}
\arguments{
This function can only be used when the extra variation is included in the response model.
\item{runInfoObj}{Object of type runInfoObj}
}
\value{
For each sweep this function outputs the ratio between the variance of the thetas' and the sum of the variances of the thetas' and the extra variation epsilon as described in Liverani et al. (2013). 
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}



