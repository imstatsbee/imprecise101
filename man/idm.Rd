% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/walley1996.R
\name{idm}
\alias{idm}
\alias{hpd}
\title{Imprecise Dirichlet Model}
\usage{
idm(nj, s = 1, N, tj = NA_real_, k, cA = 1)

hpd(
  alpha = 3,
  beta = 5,
  p = 0.95,
  tolerance = 1e-04,
  maxiter = 100,
  verbose = FALSE
)
}
\arguments{
\item{nj}{number of observations in the j th category}

\item{s}{learning parameter}

\item{N}{total number of drawings}

\item{tj}{mean probability associated with the j th category}

\item{k}{number of elements in the sample space}

\item{cA}{the number of elements in the event A}

\item{alpha}{shape1 parameter of beta distribution}

\item{beta}{shape2 parameter of beta distribution}

\item{p}{level of credible interval}

\item{tolerance}{level of error allowed}

\item{maxiter}{maximum number of iterations}

\item{verbose}{logical option suppressing messages}
}
\value{
\code{idm} returns a list of lower and upper probabilities.
\item{p.lower}{Minimum of imprecise probabilities}
\item{p.upper}{Maximum of imprecise probabilities}
\item{v.lower}{Variance of lower bound}
\item{v.upper}{Variance of upper bound}
\item{s.lower}{Standard deviation of lower bound}
\item{s.upper}{Standard deviation of upper bound}
\item{p}{Precise probabilty}
\item{p.delta}{Degree of imprecision}

\code{hpd} gives a list of scalar values corresponding to the lower and upper bounds of highest posterior probability density region.
}
\description{
This function computes lower and upper posterior probabilities under an imprecise Dirichlet model when prior information is not available.

This function searches for the lower and upper bounds of a given level of the highest posterior density interval under the imprecise Dirichlet prior.
}
\examples{
idm(nj=1, N=6, s=2, k=4)
x <- hpd(alpha=3, beta=5, p=0.95) # c(0.0031, 0.6587) when s=2
# round(x,4); x*(1-x)^5
}
\references{
Walley, P. (1996), Inferences from Multinomial Data: Learning About a Bag of Marbles. Journal of the Royal Statistical Society: Series B (Methodological), 58: 3-34. https://doi.org/10.1111/j.2517-6161.1996.tb02065.x
}
