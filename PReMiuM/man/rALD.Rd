\name{rALD}
\alias{rALD}
\alias{qALD}
\title{Asymmetric Laplace Distribution}
\description{Random generation and quantile function for a Three Parameter Asymmetric Laplace Distribution as defined in Koenker and Machado (1999) for quantile regression with location parameter equal to mu, scale parameter sigma and skewness parameter p..}
\usage{
rALD(n, mu = 0, sigma = 1, p = 0.5)
qALD(prob, mu = 0, sigma = 1, p = 0.5, lower.tail = TRUE)
}
\arguments{
\item{prob}{Vector of probabilities.}
\item{n}{Number of observations.}
\item{mu}{Location parameter.}
\item{sigma}{Scale parameter.}
\item{p}{Skewness parameter.}
\item{lower.tail}{Logical; if TRUE (default), probabilities are P[X strictly smaller than x] otherwise, P[X > x].}
}
\details{
If mu, sigma or p are not specified they assume the default values of 0, 1 and 0.5, respectively, belonging to the Symmetric Standard Laplace Distribution denoted by ALD(0,1,0.5).

As discussed in Koenker and Machado (1999) and Yu and Moyeed (2001) we say that a random variable Y is distributed as an ALD with location parameter mu, scale parameter sigma>0 and skewness parameter p in (0,1), if its probability density function (pdf) is given by

f(y|mu,sigma,p)=p(1-p)/sigma * e^{-p_{p}(y-mu)/sigma)}

where p_p(.) is the so called check (or loss) function defined by

p_p(u)=u(p - I_{u<0})

with I_{.} denoting the usual indicator function. This distribution is denoted by ALD(mu,sigma,p) and it's p-th quantile is equal to mu.

The scale parameter sigma must be positive and non zero. The skew parameter p must be between zero and one (0<p<1).
}
\value{
The length of the result is determined by n for rALD, and is the maximum of the lengths of the numerical arguments for the other functions dALD, pALD and qALD.

}
\section{Authors}{
Silvia Liverani, Queen Mary University of London, UK. 

Maintainer: Silvia Liverani <liveranis@gmail.com>
}
\references{

Galarza Morales, C., Lachos Davila, V., Barbosa Cabral, C., and Castro Cepero, L. (2017) Robust quantile regression using a generalized class of skewed distributions. Stat,6: 113-130 \doi{10.1002/sta4.140}.

}
\examples{
is.wholenumber(4) # TRUE
is.wholenumber(3.4) # FALSE
}


