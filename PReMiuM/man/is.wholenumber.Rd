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

Liverani, S., Hastie, D. I., Azizi, L., Papathomas, M. and Richardson, S. (2013) PReMiuM: An R package for Profile Regression Mixture Models using Dirichlet Processes. \emph{Submitted}. Available at http://uk.arxiv.org/abs/1303.2836

}
\examples{
is.wholenumber(4) # TRUE
is.wholenumber(3.4) # FALSE
}


