\name{vec2mat}
\alias{vec2mat}
\title{Vector to upper triangular matrix}
\description{Function to convert a vector to an upper triangular matrix. The vector does not include the diagonal values, which are then set equal to 1 in the matrix. The matrix is filled by row. }
\usage{
vec2mat(data = NA, nrow = 1) 
}
\arguments{
\item{data}{The vector to be converted, excluding the diagonal which is set equal to 1.}
\item{nrow}{The number of rows (and columns) of the resulting matrix.}
}
\value{
The symmetric matrix. The matrix is filled by column.
}
\section{Authors}{

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}
\examples{
vec2mat(data=c(1,2,3),nrow=3)
}


