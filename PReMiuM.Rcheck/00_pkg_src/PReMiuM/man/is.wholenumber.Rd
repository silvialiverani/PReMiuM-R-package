\name{is.wholenumber}
\alias{is.wholenumber}
\title{Function to check if a number is a whole number}
\description{Function to check if a number is whole, accounting for a rounding error.}
\usage{
is.wholenumber(x, tol = .Machine$double.eps^0.5) 
}
\arguments{
\item{x}{The number to be checked.}
\item{tol}{Tolerance level.}
}
\value{
The default method for 'is.wholenumber' returns 'TRUE' if the number provided is a whole number.
}
\section{Authors}{
David Hastie, Department of Epidemiology and Biostatistics, Imperial College London, UK

Silvia Liverani, Department of Epidemiology and Biostatistics, Imperial College London and MRC Biostatistics Unit, Cambridge, UK

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Silvia Liverani, David I. Hastie, Lamiae Azizi, Michail Papathomas, Sylvia Richardson (2015). PReMiuM: An R Package for Profile Regression Mixture Models Using Dirichlet Processes. Journal of Statistical Software, 64(7), 1-30. \doi{10.18637/jss.v064.i07}.

}
\examples{
is.wholenumber(4) # TRUE
is.wholenumber(3.4) # FALSE
}


